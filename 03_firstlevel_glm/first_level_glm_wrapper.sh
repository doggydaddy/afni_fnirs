#!/bin/bash

inputdir=$1
mask=$2

fmri_files=`find $inputdir -name 'prp_*.mem.fnirs.nii'`

for i in $fmri_files; 
do
    echo "processing: $i"
    ./first_level_glm.sh $i $mask
done

# cleanup 
rm .REML_cmd
rm 3dDeconvolve.err
rm X.nocensor.xmat.1D
