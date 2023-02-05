#!/usr/bin/env python3
import requests
import json
import codecs
import argparse
import pandas as pd
import numpy as np

from pathlib import Path
from tqdm.auto import tqdm
from io import StringIO
from scipy import stats
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.MolStandardize import rdMolStandardize


def listkey_request(aid, sids_or_cids):
    # Note: the code doesn't handle errors or bad requests, etc.
    # We use the PubChem REST API to generate a ListKey request.
    # This is required in order to obtain long lists of data.
    # Official documentation:
    # https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial#section=Storing-Lists-on-the-Server

    # This function will help us retrieve rows of the requested
    # data for a specific PubChem AID.

    assert sids_or_cids in ['sids', 'cids']
    listkey_template = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/{ids}/JSON?list_return=listkey'
    listkey_response = requests.get(listkey_template.format(aid=aid, ids=sids_or_cids)).json()

    listkey = listkey_response['IdentifierList']['ListKey']
    size = listkey_response['IdentifierList']['Size']

    return listkey, size


def retrieve_csv_from_listkey(template, size, cid_col='PUBCHEM_CID'):
    # We will use the PubChem REST API to get 'chunks' (rows) of the requested
    # data as CSV files.

    retrieved = []

    # We download a maximum of 5000 rows at once. Larger values can cause errors.
    for i in tqdm(range(0, size + 1, 5000)):
        retrieved.append(requests.get(template.format(iter_start=i)).content)

    # Each chunk is loaded from the string representation into a dataframe.
    processed = []
    for csv in retrieved:
        processed_csv = pd.read_csv(StringIO(str(csv, 'utf-8')))
        processed_csv = processed_csv.rename(columns=lambda c: c.strip())
        processed_csv = processed_csv.dropna(subset=[cid_col])
        processed.append(processed_csv)

    # All chunks are joined into a single dataframe.
    df = pd.concat(processed).reset_index(drop=True)
    return df


def retrieve_aid_as_csv(aid):
    # This functions leverages PubChem's REST API to download the full data
    # table corresponding to an AID.

    print('Retrieving AID data from PubChem...')
    listkey, size = listkey_request(aid, 'sids')
    retrieve_template = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/CSV?sid=listkey&listkey={listkey_id}&listkey_start={{iter_start}}&listkey_count=5000'
    df = retrieve_csv_from_listkey(retrieve_template.format(aid=aid, listkey_id=listkey), size, cid_col='PUBCHEM_CID')

    df['CID'] = df['PUBCHEM_CID'].astype(int)
    df['Activity'] = df['PUBCHEM_ACTIVITY_OUTCOME']

    del df['PUBCHEM_CID']
    del df['PUBCHEM_ACTIVITY_OUTCOME']

    return df


def get_smiles_from_pubchem_by_cid(aid):
    # The CSV files downloaded from PubChem using retrieve_aid_as_csv() do not
    # contain SMILES strings. We use different REST calls to get a SMILES
    # string for each CID.

    print('Retrieving SMILES from PubChem...')

    listkey, size = listkey_request(aid, 'cids')
    retrieve_template = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{listkey_id}/property/CanonicalSMILES/CSV?&listkey_start={{iter_start}}&listkey_count=5000'
    smiles_df = retrieve_csv_from_listkey(retrieve_template.format(listkey_id=listkey), size, cid_col='CID')
    smiles_df['CID'] = smiles_df['CID'].astype(int)

    return smiles_df


def aggregate_sd(df, list_of_replicate_column_names, agg_fn='median'):
    # This function aggregates replicate measurements from the primary screen
    # (given by the different column names, if more than one). The functions
    # also computed a Z-Score based on the activity values (not currently
    # used in modelling).

    assert agg_fn in ['median', 'mean']
    df['SD'] = df[list_of_replicate_column_names].apply(pd.to_numeric, errors='coerce', axis=1) \
                                                 .agg(func=agg_fn, axis=1, skipna=True).astype(float)
    df['SD Z-score'] = stats.zscore(df['SD'].values, nan_policy='omit')
    return df


def aggregate_dr(df, list_of_replicate_column_names, transform_dr='no', agg_fn='median'):
    # Aggregate multiple measurements (columns) for confirmatory data, if they
    # exist. The function also takes care of converting a generic 'XC50' score
    # into the corresponding pXC50.

    assert transform_dr in ['no', 'pXC50', 'minus']
    assert agg_fn in ['median', 'mean']

    processed_dr = df[list_of_replicate_column_names].apply(pd.to_numeric, errors='coerce', axis=1) \
                                                     .agg(func=agg_fn, axis=1, skipna=True)

    if transform_dr == 'pXC50':
        df['DR'] = -np.log10(processed_dr * 1e-6)
        df['XC50'] = processed_dr
        dr_names = ['DR', 'XC50']
    elif transform_dr == 'minus':
        df['DR'] = -processed_dr
        df['Log XC50'] = processed_dr
        dr_names = ['DR', 'Log XC50']
    elif transform_dr == 'no':
        df['DR'] = processed_dr
        dr_names = ['DR']

    return df, dr_names


def reduce_df(df, column_names):
    # Reformat the dataframe to keep only the CID and the new columns.

    col_elems = dict(zip(column_names, ['mean' for _ in range(len(column_names))]))
    return df.groupby('CID').agg({'CID': 'first', **col_elems, 'Activity': 'first'}).reset_index(drop=True)


def add_mols_from_smiles(df):
    # This function generates an RDKit mol object for each molecule.

    print('Reading SMILES into RDKit...')
    df['mol'] = [Chem.MolFromSmiles(s) for s in tqdm(df['CanonicalSMILES'])]
    df['rdkit-smiles'] = [Chem.MolToSmiles(m) for m in tqdm(df['mol'])]

    return df


def largest_fragment_filtering(df, remove_if_multiple_fragments=True, fragment_list=None):
    # Removes or keeps the largest fragment for each compound, according
    # to the provided options.

    print('Largest fragment filtering...')

    # Save largest fragments
    largest_Fragment = rdMolStandardize.LargestFragmentChooser()
    df['l-fragsmiles'] = [Chem.MolToSmiles(largest_Fragment.choose(m)) for m in tqdm(df['mol'])]

    if fragment_list is not None:
        print('Comparing smaller fragments with provided list...')

        # Compute all fragments
        df['other-fragsmiles'] = [
            [Chem.MolToSmiles(fr) for fr in Chem.rdmolops.GetMolFrags(m, asMols=True)] for m in tqdm(df['mol'])
            ]
        
        # Get all fragments besides the largest one
        smaller_fragments = df.apply(
            lambda row: [s for s in row['other-fragsmiles'] if s != row['l-fragsmiles']], axis=1
            )
        
        # Log instances of fragments not in the provided fragment list
        if fragment_list is not None:
            for small_fr_list in smaller_fragments:
                for fr in small_fr_list:
                    if fr not in fragment_list:
                        print(f'Fragment {fr} not in provided list!')

    # In this case we remove all compounds where multiple fragments where
    # encountered
    if remove_if_multiple_fragments:
        print('Removing compounds with multiple fragments...')
        df = df.loc[df['rdkit-smiles'] == df['l-fragsmiles']]
        return df
    
    # Otherwise, keep only the largest fragment and recompute related data
    print('Keeping compounds with multiple fragments.')
    df.loc[:, 'rdkit-smiles'] = df['l-fragsmiles']
    df.loc[:, 'mol'] = [Chem.MolFromSmiles(s) for s in tqdm(df['rdkit-smiles'])]
    df = df.dropna(subset=['mol'])

    return df


def remove_stereoisomers(df, sort_col='SD'):
    # Removes the stereochemical information from all molecules in a dataframe
    # and keeps only the unique resulting SMILES.

    print('Removing stereoisomers...')

    def smiles_no_stereo(mol):
        rdmolops.RemoveStereochemistry(mol)
        return Chem.MolToSmiles(mol)
    
    df['sf-smiles'] = df['mol'].apply(func=smiles_no_stereo)
    df = df.sort_values(sort_col, ascending=False).drop_duplicates('sf-smiles', keep='first')

    return df


def neutralise_df(df):
    # Retain only dataframe compounds that can be recorded with a formal
    # charge of 0.

    print('Neutralising molecules...')

    def neutralise_smiles(mol):
        try:
            pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
            at_matches = mol.GetSubstructMatches(pattern)
            at_matches_list = [y[0] for y in at_matches]
            if len(at_matches_list) > 0:
                for at_idx in at_matches_list:
                    atom = mol.GetAtomWithIdx(at_idx)
                    chg = atom.GetFormalCharge()
                    hcount = atom.GetTotalNumHs()
                    atom.SetFormalCharge(0)
                    atom.SetNumExplicitHs(hcount - chg)
                    atom.UpdatePropertyCache()
            return Chem.MolToSmiles(mol)
        except Chem.rdchem.AtomValenceException:
            return None
        
    df.loc[:, 'neut-smiles'] = df['mol'].apply(func=neutralise_smiles)
    df = df.dropna(subset=['neut-smiles'])
    df = df.copy()
    df.loc[:, 'neut-mol'] = [Chem.MolFromSmiles(s) for s in tqdm(df['neut-smiles'])]
    df = df.dropna(subset=['neut-mol'])
    df.loc[:, 'neut-chr'] = [rdmolops.GetFormalCharge(m) for m in tqdm(df['neut-mol'])]
    df = df.loc[df['neut-chr'] == 0]

    return df


def filter(aid, list_of_sd_cols, list_of_dr_cols, transform_dr, aid_dr=None,
           remove_if_multiple_fragments=False, fragment_smiles=None):
    # Main processing function, taking care of downloading and filtering AID
    # data. Keeps and returns the associated metadata.

    metadata = {
        'AID': aid,
        'SD columns': list_of_sd_cols,
        'DR columns': list_of_dr_cols,
        'Separate DR': 'Yes' if aid_dr else 'No'
        }

    # Retrieve the AID data table.
    main_df = retrieve_aid_as_csv(aid)
    metadata['Unfiltered SD size'] = main_df.shape[0]

    # If the confirmatory (DR) data is in a separate AID, retrieve it.
    if aid_dr:
        dr_df = retrieve_aid_as_csv(aid_dr)
        metadata['Unfiltered DR size'] = dr_df.shape[0]

    # Clean the resulting dataframe(s) and keep only the CID and activities.
    # Again, if there is a different AID for DR, this is handled separately.
    if aid_dr:
        main_df = reduce_df(aggregate_sd(main_df, list_of_sd_cols), column_names=['SD', 'SD Z-score'])
        dr_df, dr_names = aggregate_dr(dr_df, list_of_dr_cols, transform_dr=transform_dr)
        dr_df = reduce_df(dr_df, column_names=dr_names)
    else:
        main_df = aggregate_sd(main_df, list_of_sd_cols)
        main_df, dr_names = aggregate_dr(main_df, list_of_dr_cols, transform_dr=transform_dr)
        main_df = reduce_df(main_df, column_names=['SD', 'SD Z-score', *dr_names])
        metadata['Unfiltered DR size'] = main_df[main_df['DR'].notna()].shape[0]

    # Add SMILES to the dataframe(s), based on the CIDs.
    main_smiles_df = get_smiles_from_pubchem_by_cid(aid)
    main_df = main_df.merge(right=main_smiles_df, on='CID')

    if aid_dr:
        dr_smiles_df = get_smiles_from_pubchem_by_cid(aid_dr)
        dr_df = dr_df.merge(right=dr_smiles_df, on='CID')

    # Add RDKit mol objects based on the SMILES.
    main_df = add_mols_from_smiles(main_df)
    metadata['After RDKit SD size'] = main_df.shape[0]
    if aid_dr:
        dr_df = add_mols_from_smiles(dr_df)
        metadata['After RDKit DR size'] = dr_df.shape[0]
    else:
        metadata['After RDKit DR size'] = main_df[main_df['DR'].notna()].shape[0]

    # Remove CIDs that denoted a combination of compounds.
    main_df = largest_fragment_filtering(main_df, remove_if_multiple_fragments, fragment_smiles)
    metadata['After largest fragment SD size'] = main_df.shape[0]
    if aid_dr:
        dr_df = largest_fragment_filtering(dr_df, remove_if_multiple_fragments, fragment_smiles)
        metadata['After largest fragment DR size'] = dr_df.shape[0]
    else:
        metadata['After largest fragment DR size'] = main_df[main_df['DR'].notna()].shape[0]

    # Remove stereochemistry
    main_df = remove_stereoisomers(main_df, 'SD')
    metadata['After stereoisomers SD size'] = main_df.shape[0]
    if aid_dr:
        dr_df = remove_stereoisomers(dr_df, 'DR')
        metadata['After stereoisomers DR size'] = dr_df.shape[0]
    else:
        metadata['After stereoisomers DR size'] = main_df[main_df['DR'].notna()].shape[0]

    # Keep only compounds that can be recorded with a formal charge of 0.
    main_df = neutralise_df(main_df)
    metadata['After neutralisation SD size'] = main_df.shape[0]
    if aid_dr:
        dr_df = neutralise_df(dr_df)
        metadata['After neutralisation DR size'] = dr_df.shape[0]
    else:
        metadata['After neutralisation DR size'] = main_df[main_df['DR'].notna()].shape[0]

    # Check all CIDs are unique
    print('Checking resulting CIDs are unique.')
    assert len(main_df['CID'].values) == len(set(main_df['CID'].values.tolist()))
    if aid_dr:
        assert len(dr_df['CID'].values) == len(set(dr_df['CID'].values.tolist()))

    # Make sure to only keep columns of interest.
    if aid_dr:
        main_df = main_df[['CID', 'SD', 'SD Z-score', 'Activity', 'neut-smiles']]
        dr_df = dr_df[['CID', *dr_names, 'Activity', 'neut-smiles']]
    else:
        main_df = main_df[['CID', 'SD', 'SD Z-score', *dr_names, 'Activity', 'neut-smiles']]

    merged_df = None
    if aid_dr:
        merged_df = main_df.merge(right=dr_df, on='CID')
        merged_df = merged_df[['CID', 'SD', 'SD Z-score', *dr_names, 'Activity_y', 'neut-smiles_y']]
        merged_df = merged_df.rename(columns={'Activity_y': 'Activity', 'neut-smiles_y': 'neut-smiles'})


    # Separate compounds with a DR value but without an SD value.
    dr_not_in_sd_df = None
    if not aid_dr:
        main_dr_not_sd = main_df[(main_df['SD'].isna()) & (main_df['DR'].notna())]
        metadata['# of DR compounds not in primary screen'] = main_dr_not_sd.shape[0]
        if not main_dr_not_sd.empty:
            print('There are DR compounds not in SD.')
            dr_not_in_sd_df = main_dr_not_sd

    else:
        dr_cids = set(dr_df['CID'].values.tolist())
        dr_cids_not_in_sd = dr_cids.difference(set(main_df['CID'].values.tolist()))
        metadata['# of DR compounds not in primary screen'] = len(dr_cids_not_in_sd)
        if len(dr_cids_not_in_sd) > 0:
            print('There are DR compounds not in SD.')
            dr_not_in_sd_df = dr_df[dr_df['CID'].isin(dr_cids_not_in_sd)]


    # Update metadata with additional statistics.
    metadata['Filtered SD size'] = main_df.shape[0]
    pd.set_option('use_inf_as_na', True)
    if aid_dr:
        metadata['Filtered DR size'] = merged_df.shape[0]
        merged_df = merged_df.replace([np.inf, -np.inf], np.nan, inplace=False)
        r, pvalue = stats.pearsonr(merged_df[merged_df['DR'].notna()]['SD'].values, 
                                   merged_df[merged_df['DR'].notna()]['DR'].values)
    else:
        main_df_dr = main_df[main_df['DR'].notna()]
        main_df_dr = main_df_dr.replace([np.inf, -np.inf], np.nan, inplace=False)
        metadata['Filtered DR size'] = main_df_dr.shape[0]
        r, pvalue = stats.pearsonr(main_df_dr[main_df_dr['DR'].notna()]['SD'].values,
                                   main_df_dr[main_df_dr['DR'].notna()]['DR'].values)

    metadata['Pearson r'] = float(r)
    metadata['Pearson correlation p-value'] = float(pvalue)

    return main_df, merged_df, dr_not_in_sd_df, metadata


def unescaped_str(arg_str):
    return codecs.decode(str(arg_str), 'unicode_escape')


def main():
    parser = argparse.ArgumentParser(description='Retrieve multi-fidelity HTS data sets from PubChem.')
    parser.add_argument('--AID', type=str, required=True)
    parser.add_argument('--list_of_sd_cols', type=unescaped_str, nargs='+', required=True)
    parser.add_argument('--list_of_dr_cols', type=unescaped_str, nargs='+', required=True)
    parser.add_argument('--transform_dr', type=str, required=True)
    parser.add_argument('--AID_DR', type=str, required=False)
    parser.add_argument('--save_dir', type=str, required=True)
    parser.add_argument('--remove_if_multiple_fragments', action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument('--compare_fragments_with_list', action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()
    argsdict = vars(args)

    fragment_smiles = None
    if argsdict['compare_fragments_with_list']:
        fragment_smiles = pd.read_csv('small_fragments/fragment_counts.csv')['SMILES'].values

    # Get SD and DR dataframes (if separate), a dataframe of DR compounds with
    # no SD values (if it exists), and the computed metadata.
    main_df, dr_df, dr_not_in_sd_df, metadata = filter(
        argsdict['AID'],
        list_of_sd_cols=argsdict['list_of_sd_cols'],
        list_of_dr_cols=argsdict['list_of_dr_cols'],
        transform_dr=argsdict['transform_dr'],
        aid_dr=argsdict['AID_DR'],
        remove_if_multiple_fragments=argsdict['remove_if_multiple_fragments'],
        fragment_smiles=fragment_smiles
    )

    # Generate a save directory at the indicated path.
    if argsdict['AID_DR']:
        save_dir = '{save_dir}/AID{id2}-{id}/'.format(
            save_dir=argsdict['save_dir'], id=argsdict['AID'], id2=argsdict['AID_DR']
            )
    else:
        save_dir = '{save_dir}/AID{id}/'.format(save_dir=argsdict['save_dir'], id=argsdict['AID'])

    Path(save_dir).mkdir(exist_ok=True)

    # Save dataframes to disk.
    main_df.reset_index(drop=True).to_csv(f'{save_dir}/SD.csv', index=False)
    if argsdict['AID_DR']:
        dr_df.reset_index(drop=True).to_csv(f'{save_dir}/DR.csv', index=False)
    if dr_not_in_sd_df is not None and not dr_not_in_sd_df.empty:
        dr_not_in_sd_df.reset_index(drop=True).to_csv(f'{save_dir}/DR_not_in_SD.csv', index=False)

    # Save metadata to disk.
    with open(f'{save_dir}/metadata.json', 'w') as file_out:
        json.dump(metadata, file_out)


if __name__ == "__main__":
    main()