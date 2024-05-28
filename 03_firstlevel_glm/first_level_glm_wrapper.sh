#!/bin/bash

inputdir=$1
mask=$2

fmri_files=`find $inputdir -name 'prp_*.mem.fnirs.nii'`

for i in $fmri_files; 
do
    ./first_level_glm.sh $i $mask
done
