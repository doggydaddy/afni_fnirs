## Preprocessing converted .nii data

Call preprocessing script using:

        ./preproc.sh <data directory>

For each .nii file, we perform de-meaning/scaling and prewhitening/temporal
filtering.

the preprocessed files are labelled: 

        prp_<file name>.nii