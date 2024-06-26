#!/bin/bash

input_nii=$1

# note: ABSOLUTE PATH to mask file
# not problematic due to the mask is assumed to be group common.
mask_nii=$2

i_path=$(dirname "${input_nii}")
i_name=$(basename "${input_nii}")

base_data_name=`echo ${i_name} | cut -d'_' -f2`
subj_id=`echo ${base_data_name} | cut -d'.' -f1`
type_id=`echo ${base_data_name} | cut -d'.' -f2`
task_id=`echo ${base_data_name} | cut -d'.' -f3`


# print-outs
echo "analysis path: " $i_path
echo "filename: " $i_name

echo "subj_id = " $subj_id
echo "type_id = " $type_id
echo "task_id = " $task_id

# note, add confounding variables:
# -stim_file 10 ${i_path}/${subj}_${ses}_${task}_white_matter_noHead.txt'[0]' -stim_base 10 -stim_label 10 wm \

3dDeconvolve -input ${input_nii} \
             -mask ${mask_nii} \
             -num_stimts 1 \
             -stim_times_AM1 1 ${i_path}/${subj_id}.${type_id}.${task_id}.timing.1D 'dmBLOCK' -stim_label 1 mem \
             -jobs 16 \
             -gltsym 'SYM: mem' -glt_label 1 mem \
             -fout -tout -x1D ${i_path}/X.xmat.1D -xjpeg ${i_path}/X.jpg \
             -x1D_uncensored X.nocensor.xmat.1D \
             -fitts ${i_path}/fitts.${subj_id}.${type_id}.${task_id}.nii.gz \
             -errts ${i_path}/errts.${subj_id}.${type_id}.${task_id}.nii.gz \
             -bucket ${i_path}/stats.${subj_id}.${type_id}.${task_id}.nii.gz

# this is an example if you wish to perform per-block contrasts

#3dDeconvolve -input ${input_nii} \
#    -mask ${mask_nii} \
#    -num_stimts 5 \
#    -stim_times_AM1 1 ${i_path}/${subj_id}.${type_id}.${task_id}.block1.1D 'dmBLOCK' -stim_label 1 b1 \
#    -stim_times_AM1 2 ${i_path}/${subj_id}.${type_id}.${task_id}.block2.1D 'dmBLOCK' -stim_label 2 b2 \
#    -stim_times_AM1 3 ${i_path}/${subj_id}.${type_id}.${task_id}.block3.1D 'dmBLOCK' -stim_label 3 b3 \
#    -stim_times_AM1 4 ${i_path}/${subj_id}.${type_id}.${task_id}.block4.1D 'dmBLOCK' -stim_label 2 b4 \
#    -stim_times_AM1 5 ${i_path}/${subj_id}.${type_id}.${task_id}.block5.1D 'dmBLOCK' -stim_label 3 b5 \
#    -jobs 16 \
#    -gltsym 'SYM: +0.2*b1 +0.2*b2 +0.2*b3 +0.2*b4 +0.2*b5'  -glt_label 1 mem \
#    -gltsym 'SYM: +b5 -b1'                                  -glt_label 2 b5-b1 \
#    -gltsym 'SYM: b1 -0.25*b2 -0.25*b3 -0.25*b4 -0.25*b5'   -glt_label 3 b1-b25) \
#    -fout -tout -x1D ${i_path}/X.xmat.1D -xjpeg ${i_path}/X.jpg \
#    -x1D_uncensored X.nocensor.xmat.1D \
#    -fitts ${i_path}/fitts.${pset}_${subj}_${task}.nii.gz \
#    -errts ${i_path}/errts.${pset}_${subj}_${task}.nii.gz \
#    -bucket ${i_path}/stats.${pset}_${subj}_${task}.nii.gz 
