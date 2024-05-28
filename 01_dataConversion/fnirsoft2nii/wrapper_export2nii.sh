#!/bin/bash
# use >> log.txt 2>&1 to output to log file.

data_folder=$1
subjects=`find $data_folder -mindepth 2 -maxdepth 2 -type d`
echo $subjects

touch toNii.log
for s in $subjects; do
    echo "processing: $s" 2>&1 | tee toNii.log

    IN=$s
    arrIN=(${IN//// })
    prefix=${arrIN[2]} 
    
    ./export2nii.py -f $s -tag hbt -prefix $prefix
done
