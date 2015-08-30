# NAMD Analysis job directory notes.

**v0.2**
**January 2015**
**Author: MKuiper (VLSCI) and Shane Gordon (La Trobe University)**

# Disclaimer!

You are free to use it, but it may not be entirely suitable for what you are
trying to achieve. Please email feedback, bugs or suggestions to:
<a href="mailto:mkuiper@unimelb.edu.au">mkuiper@unimelb.edu.au</a> or
<a href="mailto:se2gordon@students.latrobe.edu.au">se2gordon@students.latrobe.edu.au</a>

# How it all works

The `/Analysis` directory is designed to be used at the completion of a run as a
place to consolidate data and process the trajectory data.

There are a number of sub_directories and scripts here to help you do this:

```sh
/Data
/temp_data_dir
/Scripts
/protein_cluster_data
/ligand_cluster_data

- a1_extract_all_my_data.sh       - script to extract paths to dcd data in
/MainJob_dir

- a2_create_no_H_no_H2O_dcd.sh    - script to create a greatly reduced data
file containing not water or hydrogen.

- a3_protein_backbone_cluster_analysis.sh


- a4_ligand_cluster_analysis.sh

- clustering_configuration.tcl   - configuration script for analysis.
```

The scripts above are designed to be run in order,  a1_, a2_, ..etc but certain
analysis can be omitted if necessary.
