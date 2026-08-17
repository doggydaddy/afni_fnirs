# conversion from fnirsoft export to .nii format

# dependencies

* bash
* python3, with the following packages:
    - numpy
    - pandas
    - argparse
    - glob, os, sys
    - re
    - math
* AFNI programs

Usage:

In terminal, call wrapper script:

        ./wrapper_export2nii.sh <data_directory>

Assuming <data_directory> is organized in the following manner:

- data
    - group\_1
        - subj\_1 
            - files ...
        - subj\_2
    - group\_2
        - subj\_1 
        - subj\_2
    - ...

# Notes

The wrapper script will NOT work if the subject directories lies directly under
main data directory.

To export different data types [hbo, hbr, hbt, oxy], just comment/uncomment
relevant lines in the wrapper script since it is small and one can find it
easily.

This script takes a LONG time to run, mostly because of manipulation of many
temporary files during the conversion process. Please run this script on a fast
hard-drive as loading the temporary data to memory isn't feasible in one of the
crucial steps.

Pass `-skip_nii` to regenerate only timing.txt/onsetdur.txt (skips the slow
per-timepoint .nii reconstruction). Useful for reprocessing timing metadata
alone, e.g. after a fix to the marker/timing logic - the .nii export doesn't
depend on it and doesn't need to be redone.

## timing offset (fixed 2026-08-17)

marker2onsetdur() writes onset times relative to marker -5 (task start),
matching frame 0 of the exported .nii (data2nii() starts the export there).
Downstream 3dDeconvolve calls use -local_times, which expects onsets
relative to the start of the analyzed run - so this alignment matters. See
../../notes.txt (2026-08-17) for the bug this fixed and its impact.