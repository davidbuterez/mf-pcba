
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

## MF-PCBA datasets

Scripts to download, filter, and assemble all 60 multi-fidelity datasets (individually) are provided in the directory `retrieve-scripts`.  The run time for a dataset with 335,445 SD molecules (AID 504329) is **554.52 seconds** (9 minutes and 15 seconds) on a high-end workstation with a fast internet connection. The scripts can easily be called in parallel, for example using SLURM with multiple nodes.


## Requirements

The code requires the following Python libraries: `pandas`, `numpy`, `rdkit`, `tqdm`, `scipy`.