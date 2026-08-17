# Performing 1st-level analysis using AFNI 3dDeconvolve

## model preparations

First, we need to construct our model:

### convert FSL to AFNI format

use 

                ./onsetDur2afniDecon.sh <data directory>

to convert timing.txt files in FSL onset-duration format to .timing.1D format
that 3dDeconvolve can parse.

### split into individual blocks

use:

                ./splitTimingIntoBlocks.sh <data directory> 

to convert the single .timing.1D model into individual blocks. Note that this
script only does splitting and will keep the block number amplitude moderation.

### split into individual blocks with corsi score amplitude modulation

if you want to use the corsi score for to modulate the block amplitudes, use:

                ./remodulateBlocks.sh <data directory> <corsi score csv>

**instead** of *splitTimingIntoBlocks.sh*. The script will read the timing.1D
file and split them into separate blocks. The block amplitudes will be taken
from the csv file. 

Note: If the script doesn't manage to find the patient in the csv file, then it
will use the block number as amplitude modulation as a fall-back automatically.

### split into 4 individual blocks, no amplitude modulation

if you want plain, non-amplitude-modulated block regressors (static
amplitude 1), use:

                ./splitTimingIntoStaticBlocks.sh <data directory>

**instead** of *splitTimingIntoBlocks.sh* / *remodulateBlocks.sh*. The first
entry in each *.timing.1D file is a non-task block and is dropped; the
remaining 4 blocks are written out as na.block1.1D - na.block4.1D (blocks
1-2 = easy, blocks 3-4 = hard).

## running the GLM

For the corsi-score / block-number amplitude-modulated model, use
*first_level_glm.sh* (or *first_level_glm_wrapper.sh* to run it over a
whole data directory).

For the non-amplitude-modulated 4-block model (na.block1-4.1D, produced by
*splitTimingIntoStaticBlocks.sh*), use *first_level_glm_na.sh* (or
*first_level_glm_na_wrapper.sh*). It fits one regressor per block (b1-b4)
and computes two contrasts:

  - na-hard-easy: `SYM: -0.5*b1 -0.5*b2 +0.5*b3 +0.5*b4`
  - na-corsi:     `SYM: +0.25*b1 +0.25*b2 +0.25*b3 +0.25*b4`

Output is written to stats.<subj>.<type>.<task>.na.nii.gz, distinct from the
amplitude-modulated model's stats.<subj>.<type>.<task>.nii.gz.