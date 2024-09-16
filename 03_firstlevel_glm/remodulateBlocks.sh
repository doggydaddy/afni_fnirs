#!/bin/bash

# remodulate blocks

# parsing arguments
datadir=$1          # 1. <data directory>
csv_file=$2         # 2. <corsi results csv>

# filenames for the models generated
output_affix="model.b2-5.1D"        # this is used for the full model
output_affix1="model.b2-3.1D"       # used for first split model
output_affix2="model.b4-5.1D"       # used for the second split model

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
            #echo $b > ${i_path}/${subj_id}.${type_id}.${task_id}.model.mD
        elif [ "$c" -eq 2 ]; then
            echo $modded_block > ${i_path}/${subj_id}.${type_id}.${task_id}.model.mD
        elif [ "$c" -eq 5 ]; then
            echo $modded_block >> ${i_path}/${subj_id}.${type_id}.${task_id}.model.mD 
        else
            echo $modded_block >> ${i_path}/${subj_id}.${type_id}.${task_id}.model.mD 
        fi
        c=$(($c+1))
    done

    # save one file
    single_line=`awk '{printf "%s ", $0}' ${i_path}/${subj_id}.${type_id}.${task_id}.model.mD`
    echo $single_line > ${i_path}/${subj_id}.${type_id}.${task_id}.${output_affix}

    # split and save multiple files, 
    # uncomment these if you want to split the single model into two
    # (currently using the second space as delimiter)

    # Split the input at the second space and save the two parts into two variables.
    line1=$(echo "$single_line" | sed 's/ /§/2; s/§.*//')    # Extract text before the second space.
    line2=$(echo "$single_line" | sed 's/^[^ ]* [^ ]* //')   # Extract text after the second space.
    # Save each line into a separate file.
    echo "$line1" > ${i_path}/${subj_id}.${type_id}.${task_id}.${output_affix1}
    echo "$line2" > ${i_path}/${subj_id}.${type_id}.${task_id}.${output_affix2}
    
    # cleanup
    rm ${i_path}/${subj_id}.${type_id}.${task_id}.model.mD
done
