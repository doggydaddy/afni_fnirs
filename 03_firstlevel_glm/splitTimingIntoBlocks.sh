#!/bin/bash
datadir=$1
filelist=`find $datadir -name '*.timing.1D'`

for f in $filelist; do
    input=$f

    i_path=$(dirname "${input}")
    i_name=$(basename "${input}")

    echo $i_name
    subj_id=`echo ${i_name} | cut -d'.' -f1`
    type_id=`echo ${i_name} | cut -d'.' -f2`
    task_id=`echo ${i_name} | cut -d'.' -f3`

    blocks=`cat $input`
    c=1
    for b in $blocks
    do
        echo $b "into: " > ${i_path}/${subj_id}.${type_id}.${task_id}.block${c}.1D
        c=$(($c+1))
    done
done
