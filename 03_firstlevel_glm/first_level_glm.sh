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

#3dDeconvolve -input ${input_nii} \
#             -mask ${mask_nii} \
#             -num_stimts 1 \
#             -stim_times_AM1 1 ${i_path}/${subj_id}.${type_id}.${task_id}.timing.1D 'dmBLOCK' -stim_label 1 mem \
#             -jobs 16 \
#             -gltsym 'SYM: mem' -glt_label 1 mem \
#             -fout -tout -x1D ${i_path}/X.xmat.1D -xjpeg ${i_path}/X.jpg \
#             -x1D_uncensored X.nocensor.xmat.1D \
#             -fitts ${i_path}/fitts.${subj_id}.${type_id}.${task_id}.nii.gz \
#             -errts ${i_path}/errts.${subj_id}.${type_id}.${task_id}.nii.gz \
#             -bucket ${i_path}/stats.${subj_id}.${type_id}.${task_id}.nii.gz


# this is an example if you wish to perform per-block contrasts
# use this you wish to perform per-block corsi score 
# modulated amplitude block in your contrasts

#3dDeconvolve -input ${input_nii} \
#    -polort 5 -GOFORIT 6 -noFDR \
#    -mask ${mask_nii} \
#    -num_stimts 5 \
#    -stim_times_AM1 1 ${i_path}/${subj_id}.${type_id}.${task_id}.mblock1.1D 'dmBLOCK(1)' -stim_label 1 'b1' -global_times \
#    -stim_times_AM1 2 ${i_path}/${subj_id}.${type_id}.${task_id}.mblock2.1D 'dmBLOCK(1)' -stim_label 2 'b2' -global_times \
#    -stim_times_AM1 3 ${i_path}/${subj_id}.${type_id}.${task_id}.mblock3.1D 'dmBLOCK(1)' -stim_label 3 'b3' -global_times \
#    -stim_times_AM1 4 ${i_path}/${subj_id}.${type_id}.${task_id}.mblock4.1D 'dmBLOCK(1)' -stim_label 4 'b4' -global_times \
#    -stim_times_AM1 5 ${i_path}/${subj_id}.${type_id}.${task_id}.mblock5.1D 'dmBLOCK(1)' -stim_label 5 'b5' -global_times \
#    -jobs 16 \
#    -gltsym 'SYM: +0.25*b2 +0.25*b3 +0.25*b4 +0.25*b5' -glt_label 1 b2-5 \
#    -fout -tout -x1D ${i_path}/X.xmat.1D -xjpeg ${i_path}/X.jpg \
#    -x1D_uncensored X.nocensor.xmat.1D \
#    -bucket ${i_path}/stats.${subj_id}_${type_id}_${task_id}.nii.gz 


#3dDeconvolve -input ${input_nii} \
#             -polort 6 -GOFORIT 6 -noFDR \
#             -mask ${mask_nii} \
#             -num_stimts 1 \
#             -stim_times_AM1 1 ${i_path}/${subj_id}.${type_id}.${task_id}.model.1D 'dmBLOCK' -stim_label 1 corsi -local_times \
#             -jobs 16 \
#             -gltsym 'SYM: corsi' -glt_label 1 corsi \
#             -tout \
#             -bucket ${i_path}/stats.${subj_id}.${type_id}.${task_id}.nii.gz \


3dDeconvolve -input ${input_nii} \
             -polort 6 -GOFORIT 6 -noFDR \
             -mask ${mask_nii} \
             -num_stimts 2 \
             -stim_times_AM1 1 ${i_path}/${subj_id}.${type_id}.${task_id}.model.b2-3.1D 'dmBLOCK' -stim_label 1 corsi_easy -local_times \
             -stim_times_AM1 2 ${i_path}/${subj_id}.${type_id}.${task_id}.model.b4-5.1D 'dmBLOCK' -stim_label 2 corsi_hard -local_times \
             -jobs 16 \
             -gltsym 'SYM: -corsi_easy +corsi_hard' -glt_label 1 hard-easy \
             -gltsym 'SYM: +0.5*corsi_easy +0.5*corsi_hard' -glt_label 2 corsi \
             -tout \
             -bucket ${i_path}/stats.${subj_id}.${type_id}.${task_id}.nii.gz \


# cleanup
rm Decon*
