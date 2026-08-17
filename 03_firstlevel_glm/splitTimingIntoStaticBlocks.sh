#!/bin/bash

# splitTimingIntoStaticBlocks.sh
#
# split a subject's *.timing.1D file into 4 individual, non-amplitude-
# modulated block regressors (static amplitude = 1).
#
# note: the output keeps the 'onset*1:duration' AM1 syntax (amplitude
# pinned to 1) rather than bare 'onset:duration', because 3dDeconvolve's
# plain -stim_times does not accept a per-event duration field - only
# -stim_times_AM1 does. Since each block regressor here has exactly one
# event, a fixed amplitude of 1 is equivalent to no modulation.
#
# The first entry in *.timing.1D is a non-task block and is dropped, same
# convention as remodulateBlocks.sh. The remaining 4 blocks are renumbered
# 1-4 (old block2 -> na.block1, old block3 -> na.block2, old block4 ->
# na.block3, old block5 -> na.block4), so that:
#   - na.block1 + na.block2 = easy blocks
#   - na.block3 + na.block4 = hard blocks
#
# Recordings are frequently shorter than the full 5-block corsi timing -
# most severely for the last block (na.block4): as of 2026-08-17, only 4 of
# 439 subjects have it fully captured, and 70 have essentially none of it
# (see ../notes.txt). A block's full intended duration would produce a
# wildly unstable dmBLOCK regressor when most of that duration is actually
# missing from the recording (the column ends up scaled to model a ~100s
# response off a handful of captured seconds). So each block's modeled
# duration is clipped here to whatever was actually captured
# (min(onset+duration, scan_duration) - onset). Blocks with less than
# min_capture_frac of their intended duration captured are dropped
# entirely (no na.blockN.1D written) - first_level_glm_na.sh's
# 3dDeconvolve call then fails cleanly for that subject/type from the
# missing stim_times file, which is how those runs get excluded from
# downstream analysis (no stats.*.na.nii.gz produced).
#
# Use this **instead** of splitTimingIntoBlocks.sh / remodulateBlocks.sh
# when you want plain dmBLOCK regressors with no amplitude modulation.

datadir=$1          # 1. <data directory>
min_capture_frac=0.10

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

    prp_nii="${i_path}/prp_${subj_id}.${type_id}.${task_id}.fnirs.nii"
    if [ ! -f "$prp_nii" ]; then
        echo "  [WARN] no matching fnirs nii (${prp_nii}), skipping"
        continue
    fi
    nt=`3dinfo -nt "$prp_nii"`
    tr=`3dinfo -tr "$prp_nii"`
    scan_dur=`python3 -c "print(${nt}*${tr})"`

    blocks=`cat $input`
    c=1
    nb=1
    for b in $blocks
    do
        if [ "$c" -eq 1 ]; then
            echo "  skip first block (non-task block)"
            c=$(($c+1))
            continue
        fi

        onset=`echo $b | sed -E 's/([0-9.]+)\*.*/\1/'`
        dur=`echo $b | sed -E 's/.*:([0-9.]+)/\1/'`

        read captured_dur frac <<< `python3 -c "
onset=${onset}; dur=${dur}; scan_dur=${scan_dur}
captured = max(0.0, min(onset+dur, scan_dur) - onset)
frac = captured/dur if dur > 0 else 0.0
print(f'{captured:.3f} {frac:.4f}')
"`

        below=`python3 -c "print(1 if ${frac} < ${min_capture_frac} else 0)"`
        outfile="${i_path}/${subj_id}.${type_id}.${task_id}.na.block${nb}.1D"
        rm -f "$outfile"
        if [ "$below" -eq 1 ]; then
            echo "  na.block${nb}: only ${frac} of intended duration captured - dropping (no file written)"
        else
            echo "${onset}*1:${captured_dur}" > "$outfile"
        fi

        nb=$(($nb+1))
        c=$(($c+1))
    done
done
