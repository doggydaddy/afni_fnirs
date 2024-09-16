#!/bin/bash

# usage/call:
#
# ./preprocessing_nii <data_directory>
#
# pre-processes ALL *.nii files in <data_directory>

data_folder=$1                                # first argument is <data_directory>
all_niis=`find ${data_folder} -name '*.nii'`  # grab ALL files with .nii extension and puts them in a list
echo $all_niis

for nii in $all_niis; 
do
    # for each .nii file, perform the following:

    # first, we need to get the filename and file path
    filepath=$(dirname "$nii")
    filename=$(basename "$nii")
    
 
    # performing demean/scaling
    # masking to only the 16 channels to avoid processing empty voxels
    #
    # We will also have to scale the data to ensure that the mean BOLD signal
    # across the run is 100, so that any deflections can be expressed in
    # percent signal change.
    #
    3dTstat -prefix ${filepath}/mean_${filename} $nii
    3dcalc -a $nii -b ${filepath}/mean_${filename} -c fNIRS_template_mask.nii -expr 'c * min(200, a/b*100) * step(b)' -prefix ${filepath}/scl_${filename}  


    # performing prewhitening/temporal filtering 
    # there are many ways to skin a cat, 
    # but here is one way that is qutie inuitive programmatically
    # additional things we want this program can do:
    # 1. despike prior to band-pass to get rid of massive spikes
    # 2. normalize all time-series to L2 norm = 1
    # masking to only the 16 channels to avoid processing empty voxels
    # 
    3dBandpass -despike -norm -mask fNIRS_template_mask.nii -band 0.08 0.12 -input ${filepath}/scl_${filename} -prefix ${filepath}/prp_${filename}
done
