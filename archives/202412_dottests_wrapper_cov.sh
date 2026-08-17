#!/bin/bash

# main analysis wrapper for 2024-12.
# only hbo
# covariate list only consists of sans saps and dose
# b45-23 model

# hbo
# lrv-pcon
./gen_ttest_cov.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   lrv.hbo.sublist.txt \
			   pcon.hbo.sublist.txt \
			   ../../data_aux/select.cov.1D \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbo_lrv-pcon

# lrv-con
./gen_ttest_cov.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   lrv.hbo.sublist.txt \
			   con.hbo.sublist.txt \
			   ../../data_aux/select.cov.1D \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbo_lrv-con

# con-pcon
./gen_ttest_cov.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   con.hbo.sublist.txt \
			   pcon.hbo.sublist.txt \
			   ../../data_aux/select.cov.1D \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbo_con-pcon
