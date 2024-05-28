#!/bin/bash

3dttest++ -prefix test_mem_friska-pat \
          -mask biopac16ch_template_mask.nii \
          -AminusB \
		-setA test_kon \
      stats.2011K.hbr.mem.nii.gz ../03_firstlevel_glm/testdata/friska_kontroller/2011K/stats.2011K.hbr.mem.nii.gz'[4]' \
      stats.2011K.hbo.mem.nii.gz ../03_firstlevel_glm/testdata/friska_kontroller/2011K/stats.2011K.hbo.mem.nii.gz'[4]' \
      stats.2047K.hbo.mem.nii.gz ../03_firstlevel_glm/testdata/friska_kontroller/2047K/stats.2047K.hbo.mem.nii.gz'[4]' \
      stats.2047K.hbr.mem.nii.gz ../03_firstlevel_glm/testdata/friska_kontroller/2047K/stats.2047K.hbr.mem.nii.gz'[4]' \
      stats.2034K.hbo.mem.nii.gz ../03_firstlevel_glm/testdata/friska_kontroller/2034K/stats.2034K.hbo.mem.nii.gz'[4]' \
      stats.2034K.hbr.mem.nii.gz ../03_firstlevel_glm/testdata/friska_kontroller/2034K/stats.2034K.hbr.mem.nii.gz'[4]' \
		-setB test_pat \
      stats.0049P.hbr.mem.nii.gz ../03_firstlevel_glm/testdata/patienter_lrv/0049P/stats.0049P.hbr.mem.nii.gz'[4]' \
      stats.0049P.hbo.mem.nii.gz ../03_firstlevel_glm/testdata/patienter_lrv/0049P/stats.0049P.hbo.mem.nii.gz'[4]' \
      stats.0034P.hbo.mem.nii.gz ../03_firstlevel_glm/testdata/patienter_lrv/0034P/stats.0034P.hbo.mem.nii.gz'[4]' \
      stats.0034P.hbr.mem.nii.gz ../03_firstlevel_glm/testdata/patienter_lrv/0034P/stats.0034P.hbr.mem.nii.gz'[4]' \
      stats.0003P.hbo.mem.nii.gz ../03_firstlevel_glm/testdata/patienter_lrv/0003P/stats.0003P.hbo.mem.nii.gz'[4]' \
      stats.0003P.hbr.mem.nii.gz ../03_firstlevel_glm/testdata/patienter_lrv/0003P/stats.0003P.hbr.mem.nii.gz'[4]' \
