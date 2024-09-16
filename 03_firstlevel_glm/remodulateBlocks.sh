#!/bin/bash

# remodulate blocks

# parsing arguments
datadir=$1          # 1. <data directory>
csv_file=$2         # 2. <corsi results csv>

filelist=`find $datadir -name '*.timing.1D'`

for f in $filelist
do
    input=$f
    echo "processing: $input"

    i_path=$(dirname "${input}")
    i_name=$(basename "${input}")

    subj_id=`echo ${i_name} | cut -d'.' -f1`
    type_id=`echo ${i_name} | cut -d'.' -f2`
    task_id=`echo ${i_name} | cut -d'.' -f3`

    blocks=`cat $input`
    c=1
    for b in $blocks
    do
        col_nr=$(($c))
        amp=`awk -F, '$1 == "'$subj_id'" {print $'$col_nr'}' $csv_file`
        if [ -z "$amp" ]; then
            modded_block=`echo $b`
        else
            modded_block=`echo $b | sed 's/\*\([0-9]\+\):/*'$amp':/'`
        fi

        if [ "$c" -eq 1 ]; then
            echo "skip first block"
        elif [ "$c" -eq 2 ]; then
            echo $modded_block > ${i_path}/${subj_id}.${type_id}.${task_id}.model.mD
        elif [ "$c" -eq 5 ]; then
            echo $modded_block >> ${i_path}/${subj_id}.${type_id}.${task_id}.model.mD 
        else
            echo $modded_block >> ${i_path}/${subj_id}.${type_id}.${task_id}.model.mD 
        fi
        c=$(($c+1))
    done
    awk '{printf "%s ", $0}' ${i_path}/${subj_id}.${type_id}.${task_id}.model.mD > ${i_path}/${subj_id}.${type_id}.${task_id}.model.1D
    rm ${i_path}/${subj_id}.${type_id}.${task_id}.model.mD
done
