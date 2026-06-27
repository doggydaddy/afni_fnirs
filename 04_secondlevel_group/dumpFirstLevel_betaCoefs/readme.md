# DumpFirstLevel_allChannels

Dump the first-level beta coefficients onto CSV files for own 2nd level analysis without using AFNI 3dttests:

main program:

        ./dumpBetaCoefs_allChannels.sh <group1 sublist> <group1 label> <group2 sublist> <group2 label> <volume idx> <output filename>

the following support programs are needed for the main program to run:

* grab_betas.sh
  * Grabs 1st level beta coefs from 1 channel

* add_covariates.py
  * Adds covariates in covariates file (path hard-coded in the support program) to the csv

* merge_channel_dumps.py
  * Takes temporary csv single channel dumps from the output of grab_betas.sh
    and merge them into the final output csv

The post-processing program relabel_groups.sh can be called on the final output csv to relabel the csv.
