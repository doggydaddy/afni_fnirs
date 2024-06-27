#!/bin/bash

./gen_ttest.sh testdata \
		biopac16ch_template_mask.nii \
		con.hbr.sublist.txt \
		lrv.hbr.sublist.txt \
		no \
		4 \
		4 \
		ttest \
	    mem_hbr_con-lrv

