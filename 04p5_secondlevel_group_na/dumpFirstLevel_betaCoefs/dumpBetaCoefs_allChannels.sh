#!/bin/bash

# non-amplitude-modulated (na) model bucket sub-brick indices (0-based,
# confirmed via `3dinfo -label`):
#   9  = na-hard-easy_GLT#0_Coef
#   11 = na-corsi_GLT#0_Coef

group1_subjectlist=$1
group1_groupname=$2
group2_subjectlist=$3
group2_groupname=$4

volume_index=$5

output_filename=$6

template_mask='../../templates/biopac16ch_template_mask.nii'
# local, regenerated covariate file (includes wais_matrix/sud) - see
# prep_select_cov_na.py; the canonical data_aux/select.cov.1D on this mount
# is stale (missing those columns).
covariate_file='select.cov.na.1D'


for((channel_nr=1; channel_nr<=16; channel_nr++))
do
    ./grab_betas.sh ${group1_subjectlist} ${group1_groupname} ${group2_subjectlist} ${group2_groupname} ${volume_index} ${template_mask} ${channel_nr}
    python3 add_covariates.py -input output.txt -covar ${covariate_file} -output output.csv
    mv output.csv TMP_ch${channel_nr}.csv
    ./relabel_groups.sh TMP_ch${channel_nr}.csv
    rm output.txt
done

python3 merge_channel_dumps.py
mv merged.csv ${output_filename}

rm TMP_ch*.csv
