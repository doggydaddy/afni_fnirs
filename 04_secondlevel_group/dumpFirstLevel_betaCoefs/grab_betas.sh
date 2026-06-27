#!/bin/bash

# from 2 files list (and their group labels)
# grabs values (beta coefficients) from a given volume
# and a channel (mask and channel number)

files_list_1=$1
files_list_1_group=$2
files_list_2=$3
files_list_2_group=$4

volume_nr=$5

channel_mask=$6
channel_nr=$7

sed "s/$/[$volume_nr]/" $files_list_1 > tmp_file_list_1.txt
sed "s/$/[$volume_nr]/" $files_list_2 > tmp_file_list_2.txt

3dcalc -a $channel_mask -expr "within(a,$channel_nr,$channel_nr)" -prefix tmp_channel_mask.nii

3dROIstats -mask tmp_channel_mask.nii `cat tmp_file_list_1.txt` > tmp_g1.txt
sed -i "s/$/\t$files_list_1_group/" tmp_g1.txt

3dROIstats -mask tmp_channel_mask.nii `cat tmp_file_list_2.txt` > tmp_g2.txt
sed -i "s/$/\t$files_list_2_group/" tmp_g2.txt

{ cat tmp_g1.txt; tail -n +2 tmp_g2.txt; } > output.txt

rm tmp_channel_mask.nii
rm tmp_file_list_1.txt
rm tmp_file_list_2.txt
rm tmp_g1.txt tmp_g2.txt