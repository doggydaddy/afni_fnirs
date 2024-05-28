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