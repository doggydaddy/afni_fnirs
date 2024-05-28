#!/bin/bash

datadir=$1
mask=$2
gA_filelist=$3
gB_filelist=$4
is_paired=$5
volA=$6
volB=$7
outputdir=$8
prefix=$9

touch cmd.sh

if [[ "${is_paired}" == "paired" ]]; then
    cat ttest_paired_header.txt | sed "s/REPLACE/${prefix}/g"  > cmd.sh
else
    cat ttest_header.txt | sed "s/REPLACE/${prefix}/g" > cmd.sh
fi
sed -i 's/MASK/'"${mask}"'/g' cmd.sh

gA=`cat ${gA_filelist}`
gB=`cat ${gB_filelist}`
gA_name=`basename ${gA_filelist}`
gB_name=`basename ${gB_filelist}`
gA_prefix=`echo ${gA_name} | cut -d'.' -f1`
gB_prefix=`echo ${gB_name} | cut -d'.' -f1`

echo "		-setA $gA_prefix \\" >> cmd.sh
for i in $gA; do
    i_path=$(dirname "${i}")
    i_name=$(basename "${i}")
    subj=`echo ${i_name} | cut -d'_' -f2`
    echo "      ${subj} $i'[$volA]' \\" >> cmd.sh
done

echo "		-setB $gB_prefix \\" >> cmd.sh
for i in $gB; do
    i_path=$(dirname "${i}")
    i_name=$(basename "${i}")
    subj=`echo ${i_name} | cut -d'_' -f2`
    echo "      ${subj} $i'[$volB]' \\" >> cmd.sh
done

# finishing touches
sed '$ s/.$//' cmd.sh
chmod +x cmd.sh

# print command
cat cmd.sh
# run command
./cmd.sh

# make output directory if doesn't already exist, 
# and copy t-test results (including the run command) to it
mkdir -p $outputdir/$prefix/
mv cmd.sh $outputdir/$prefix/
mv ${prefix}* $outputdir/$prefix/
