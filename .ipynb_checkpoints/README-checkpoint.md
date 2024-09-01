# Ancestral sequence reconstruction for biocatalytic stereoselective synthesis of azaphilones

This repo contains all the required data to replicate the molecular modelling of the acyltranferase ancestors.

## Install

Create the conda environment called anc_acyl

`conda env create -f environment.yml`

`conda activate anc_acyl`

Install [MMTSB Toolset](https://github.com/mmtsb/toolset)

Install [pyCHARMM](https://github.com/BrooksResearchGroup-UM/pyCHARMM-Workshop/tree/main/0Install_Tools)

In most scripts, the CHARMM_LIB_DIR is a required variable and is the lib directory generated when bulding pyCHARMM.
An example user CHARMM_LIB_DIR is `/home/azamh/charmm/pycharmm_2024/pycharmm_install/lib` which appears in the example scripts.

## Use

Each folder contains the data generated from modelling/docking and all the data necessary to reproduce our results.

Each folder (except [AlphaFold](https://github.com/google-deepmind/alphafold)), contains an **example** directory. The example directory contains an example script or notebook which should be run with the anc_acyl enviornment and the example directory as the current working directory. Modify the example scripts with your CHARMM_LIB_DIR. All generated output is placed in the **output** folder. See this folder for example output and logs before running. Each example can be run independently with the provided data in the repository.

An important note is that scripts use absolute paths and paths for CHARMM should be lowercase and less than 128 characters.

## Order

To reproduce the results from the paper, the following modelling order should be followed:

0. alphafold_runs
1. vacuum_min
2. add_5cthioester
3. rcdocker
4. fcdocker
5. fcdocker_flexible_5cthioester
6. add_substrate
7. dyn
8. analysis
