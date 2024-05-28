#!/bin/bash

datadir=$1
onset_durations=`find ${datadir} -name '*.onsetdur.txt'`

for f in $onset_durations;
do
    i_path=$(dirname "${f}")
    i_name=$(basename "${f}")
    base_data_name=`echo ${i_name} | cut -d'_' -f2`
    subj_id=`echo ${base_data_name} | cut -d'.' -f1`
    type_id=`echo ${base_data_name} | cut -d'.' -f2`
    task_id=`echo ${base_data_name} | cut -d'.' -f3`
    o_path=${i_path}/${subj_id}.${type_id}.${task_id}.timing.1D

    echo "i_path: " $i_path
    echo "i_name: " $i_name
    echo "o_path: " $o_path

	timing_tool.py -fsl_timing_files $f -write_timing ${o_path}
done
