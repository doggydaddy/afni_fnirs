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