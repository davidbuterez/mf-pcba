
# MF-PCBA

This repository contains the code required to download, filter, and assemble multi-fidelity datasets from PubChem.

The main functionality is provided in the script `pubchem_retrieve.py`. The script takes the following arguments:

-  `--AID`, this corresponds to the AID of the SD dataset on PubChem

-  `--list_of_sd_cols`, this corresponds to all of the columns that contain SD activity values (e.g. activity, inhibition, etc.). Should include columns for replicates

-  `--list_of_dr_cols`, this corresponds to all of the columns that contain DR activity values (e.g. activity, inhibition, etc.)

-  `--transform_dr`, this allows the conversion of DR values from XC50 to pXC50 (`transform_dr "pXC50"`), from a Log XC50 value to pIC50 (`transform_dr "minus"`), or no transformation (`transform_dr "no"`)

-  `--AID_DR`, optionally, if the DR dataset is reported in a different PubChem assay, this corresponds to the AID of the DR dataset

-  `--save_dir`, path to a directory where the resulting datasets will be stored; will be created if it does not exist


Example invocations:

1. SD and DR data are stored in a single PubChem assay:

```python pubchem_retrieve.py --AID "1445" --list_of_sd_cols "Primary Inhibition" "Primary Inhibition Rep 2" "Primary Inhibition Rep 3" --list_of_dr_cols "IC50" --transform_dr "pXC50" --save_dir retrieved```

2. SD and DR data are stored in separate PubChem assays:

```python pubchem_retrieve.py --AID "873" --list_of_sd_cols "Percent inhibition" --list_of_dr_cols "IC50 #1" --transform_dr "pXC50" --AID_DR "1431" --save_dir retrieved```

The notebook `add_default_pXC50.ipynb` can be used to add the default pXC50 values mentioned in the paper and other useful information for machine learning.

### Provided scripts to download the MF-PCBA datasets
Scripts to download, filter, and assemble all 60 multi-fidelity datasets (individually) are provided in the directory `retrieve-scripts` (make sure to update with a valid save directory).  The run time for a dataset with 335,445 SD molecules (AID 504329) is **554.52 seconds** (9 minutes and 15 seconds) on a high-end workstation with a fast internet connection. The scripts can easily be called in parallel, for example using SLURM with multiple nodes.

### Random seeds
The five random seeds for each of the 60 multi-fidelity datasets are provided in the file `MF_PCBA_random_seeds.json`, with an example of their usage in `split_DR_with_random_seeds.ipynb`.

### Requirements
The code requires the following Python libraries: `pandas`, `numpy`, `rdkit`, `tqdm`, `scipy`.

## Example workflow
1. Download one dataset or a selection of datasets. For example, the following command downloads the AID 1445 dataset to a `save_dir` directory:

```
python pubchem_retrieve.py --AID "1445"
--list_of_sd_cols "Primary Inhibition" "Primary Inhibition Rep 2" "Primary Inhibition Rep 3" 
--list_of_dr_cols "IC50" --transform_dr "pXC50" --save_dir <save_dir>
```

2. The step above downloaded and filtered the data corresponding to AID 1445. To obtain train, validation, and test sets, the `split_DR_with_random_seeds.ipynb` notebook can be used. The same 5 random split seeds as used in the paper are provided in the MF-PCBA repository and are used by default. After this step, the DR data is split into train, validation, and test sets 5 different times, with the resulting `.csv` files being saved in different directories:

```
parent_dir/
├── 0/
│   ├── train.csv
│   ├── validate.csv
│   └── test.csv
├── 1/
│   ├── train.csv
│   ├── validate.csv
│   └── test.csv
| ...
└──
...
```
## Multi-fidelity machine learning worflows
In-depth examples of how to use the MF-PCBA datasets are provided in the repository https://github.com/davidbuterez/multi-fidelity-gnns-for-drug-discovery, corresponding to the *Improving molecular property prediction with multi-fidelity graph representation learning on millions of experimental endpoints* paper.
