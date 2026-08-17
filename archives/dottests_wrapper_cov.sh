#!/bin/bash

# hbo
# lrv-pcon
./gen_ttest_cov.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   lrv.hbo.sublist.txt \
			   pcon.hbo.sublist.txt \
			   ../../data_aux/covars.1D \
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
			   ../../data_aux/covars.1D \
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
			   ../../data_aux/covars.1D \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbo_con-pcon

# hbr
# lrv-pcon
./gen_ttest_cov.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   lrv.hbr.sublist.txt \
			   pcon.hbr.sublist.txt \
			   ../../data_aux/covars.1D \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbr_lrv-pcon

# lrv-con
./gen_ttest_cov.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   lrv.hbr.sublist.txt \
			   con.hbr.sublist.txt \
			   ../../data_aux/covars.1D \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbr_lrv-con

# con-pcon
./gen_ttest_cov.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   con.hbr.sublist.txt \
			   pcon.hbr.sublist.txt \
			   ../../data_aux/covars.1D \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbr_con-pcon