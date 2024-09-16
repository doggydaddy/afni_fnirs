#!/bin/bash

# hbr
# lrv-pcon
./gen_ttest.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   lrv.hbr.sublist.txt \
			   pcon.hbr.sublist.txt \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbr_lrv-pcon
# pcon-lrv
./gen_ttest.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   lrv.hbr.sublist.txt \
			   con.hbr.sublist.txt \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbr_lrv-con
# con-pcon
./gen_ttest.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   pcon.hbr.sublist.txt \
			   con.hbr.sublist.txt \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbr_pcon-con

# hbo
# lrv-pcon
./gen_ttest.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   lrv.hbo.sublist.txt \
			   pcon.hbo.sublist.txt \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbo_lrv-pcon
# pcon-lrv
./gen_ttest.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   lrv.hbo.sublist.txt \
			   con.hbo.sublist.txt \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbo_lrv-con
# con-pcon
./gen_ttest.sh ../../data \
			   ../templates/biopac16ch_template_mask.nii \
			   pcon.hbo.sublist.txt \
			   con.hbo.sublist.txt \
			   no \
			   5 \
			   5 \
			   ttest \
	 		   corsi_hbo_pcon-con