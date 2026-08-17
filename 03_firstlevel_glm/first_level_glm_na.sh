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

# non-amplitude-modulated model: 4 separate block regressors (b1-b4),
# each a dmBLOCK with a single event pinned to amplitude 1 (requires
# splitTimingIntoStaticBlocks.sh to have generated the na.block1-4.1D
# files beforehand, in 'onset*1:duration' AM1 syntax - plain -stim_times
# doesn't support per-event durations). b1-b2 = easy blocks, b3-b4 = hard.

# polort A: auto-select detrending order per run (pnum = 1+int(D/150), D =
# run duration in seconds). A fixed polort 6 sets the drift-removal cutoff
# at ~140-165s for these run lengths - close enough to individual block
# durations (up to ~130s) to risk eating into the task response. Auto
# roughly doubles that cutoff while still adapting to each subject's
# (varying) run length.
3dDeconvolve -input ${input_nii} \
             -polort A -GOFORIT 6 -noFDR \
             -mask ${mask_nii} \
             -num_stimts 4 \
             -stim_times_AM1 1 ${i_path}/${subj_id}.${type_id}.${task_id}.na.block1.1D 'dmBLOCK' -stim_label 1 b1 -local_times \
             -stim_times_AM1 2 ${i_path}/${subj_id}.${type_id}.${task_id}.na.block2.1D 'dmBLOCK' -stim_label 2 b2 -local_times \
             -stim_times_AM1 3 ${i_path}/${subj_id}.${type_id}.${task_id}.na.block3.1D 'dmBLOCK' -stim_label 3 b3 -local_times \
             -stim_times_AM1 4 ${i_path}/${subj_id}.${type_id}.${task_id}.na.block4.1D 'dmBLOCK' -stim_label 4 b4 -local_times \
             -jobs 16 \
             -gltsym 'SYM: -0.5*b1 -0.5*b2 +0.5*b3 +0.5*b4' -glt_label 1 na-hard-easy \
             -gltsym 'SYM: +0.25*b1 +0.25*b2 +0.25*b3 +0.25*b4' -glt_label 2 na-corsi \
             -tout \
             -bucket ${i_path}/stats.${subj_id}.${type_id}.${task_id}.na.nii.gz \


# cleanup
rm Decon*
