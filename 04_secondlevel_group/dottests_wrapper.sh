#!/bin/bash

./gen_ttest.sh testdata \
		biopac16ch_template_mask.nii \
		test_kon.subjectlist.txt \
		test_pat.subjectlist.txt \
		no \
		4 \
		4 \
		ttest \
		test_mem_friska-pat

