#!/usr/bin/python

import subprocess
import numpy as np
import pandas as pd
import argparse
import glob, os, sys
import re
import math

# script settings
np.set_printoptions(suppress=True)

# ------------------
# "static" variables
# ------------------
# paths to template files (.nii and .txt), resolved relative to this script's
# location rather than hardcoded to one host's mount name (e.g. /mnt/speyside)
# so the script runs unchanged on any mount aliasing the same karim_fnirs tree
# (highlands, bellevue, speyside, ...).
_karim_fnirs_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
template_nii_path =  os.path.join(_karim_fnirs_root, 'afni_fnirs', 'templates', 'biopac16ch_template.nii')
template_txt_path =  os.path.join(_karim_fnirs_root, 'afni_fnirs', 'templates', 'biopac16ch_template.txt')
mirrored_data_list = os.path.join(_karim_fnirs_root, 'data', 'mirror_list.txt')
# TR is set here since all datasets should in theory have the same TR
tr = 0.51

# ---------------
# argument parser
# ---------------
parser = argparse.ArgumentParser(prog='export2nii.py', 
                                 description='parse fnirsoft exports to nifti format')
parser.add_argument('-foldername', 
                    help='input directory name')
parser.add_argument('-tag', 
                    help='valid tags are hbo, hbt, hbr, and oxy')
parser.add_argument('-prefix',
                    help='prefix for output files')
parser.add_argument('-skip_nii', action='store_true',
                    help='skip writing .nii files (only (re)generate timing.txt/onsetdur.txt). '
                         'Useful for reprocessing timing metadata alone without repeating the '
                         'slow per-timepoint .nii reconstruction, e.g. after a marker2onsetdur fix.')
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

# support functions
# -----------------
# find oxy Block/Time/Marker 
def find_pairs(tag):
    tags = ['hbo', 'hbr', 'hbt', 'oxy']
    if tag == 'hbo':
        search_string = ".*hbo.Block*"
    elif tag == 'hbr':
        search_string = ".*hbr.Block*"
    elif tag == 'hbt':
        search_string = ".*hbt.Block*"
    elif tag == 'oxy':
        search_string = ".*oxy.Block*"
    else:
        print("find_pairs: unknown tag specified!")
        return( [] )

    return_object = []

    # search for Block file(s)    
    search_block = re.compile(search_string)
    search_results = list(filter(search_block.match, txt_files))

    # find Time equivalent files 
    # marker file is shared by all tags, but should have the same index
    for x in search_results:
        an_entry = []
        an_entry.append(x)
        time_equivalent = x.replace("Block", "Time")
        time_equivalent_results = [time_equivalent in txt_files]
        if sum(time_equivalent_results) == 1:
            an_entry.append(time_equivalent)
        
        marker = x.replace(tag+".Block", "Marker")
        marker_results = [marker in txt_files]
        if sum(marker_results) == 1:
            an_entry.append(marker)
        
        return_object.append(an_entry)

    return(return_object)

def parse_pair(block_file, time_file):
    blk = np.loadtxt(block_file, delimiter='\t', skiprows=14)
    blk = np.round(blk, decimals=6)
    tim = np.loadtxt(time_file, delimiter='\t', skiprows=14)

    if tim.shape[0] != blk.shape[0]:
        print("sanity check failed, number of samples and number of time-stamps does not match!")
        return 0
    if blk.shape[1] != 16:
        print("sanity check failed, number of channels is not 16")
        return 0
    
    return_object = np.c_[tim, blk]
    return( return_object )

def find_nearest(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def parse_marker(marker_file, time_array):
    mrk = np.loadtxt(marker_file, delimiter='\t', skiprows=14)
    mrk_df = pd.DataFrame({'MarkerTime': mrk[:,0], 'MarkerType': mrk[:,1]})
    filtered_mrk_df = mrk_df.drop_duplicates(subset=['MarkerType'], keep='last')

    nearest_time = []
    for index, row in filtered_mrk_df.iterrows():
        nearest_time_entry = find_nearest(time_array, row['MarkerTime'])
        nearest_time.append(nearest_time_entry)

    filtered_mrk_df.insert(2, "NearestTime", nearest_time)
    return(filtered_mrk_df) 

def trim_data(data, marker):
    start_time = marker.loc[marker['MarkerType']==-5].values[0][2]
    end_time = marker.loc[marker['MarkerType']==108].values[0][2]
    data_startidx = np.where( data[:,0] == start_time )[0][0]
    data_endidx = np.where( data[:,0] == end_time )[0][0]
    trimmed_data = data[data_startidx:data_endidx, :]
    return( trimmed_data )

def trim_marker(marker, start_marker, end_marker):
    start_idx = marker.index[marker['MarkerType'].eq(start_marker)].min()
    tm = marker[start_idx:]
    tr = tm.loc[tm['MarkerType']==end_marker].index.values[0]
    tm = tm.truncate(after=tr)
    return( tm )

def split_data(data, marker, start_marker, end_marker):
    start_time = marker.loc[marker['MarkerType']==start_marker].values[0][2]
    end_time = marker.loc[marker['MarkerType']==end_marker].values[0][2]
    data_startidx = np.where( data[:,0] == start_time )[0][0]
    data_endidx = np.where( data[:,0] == end_time )[0][0]
    output_data = data[data_startidx:data_endidx, :]
    return( output_data )

def marker2onsetdur(marker):
    # data2nii() builds the exported .nii starting at marker -5 (task start),
    # so frame 0 of the .nii corresponds to marker -5's absolute recording
    # time, not absolute time zero. First-level GLM scripts use 3dDeconvolve
    # with -local_times, which interprets onsets as relative to the start of
    # the analyzed run (frame 0) - so onsets here must be re-zeroed to marker
    # -5's time, or every downstream regressor is offset by that (per-subject,
    # 12-232s in this dataset) amount relative to where it actually occurs in
    # the exported data.
    t0 = marker.loc[marker['MarkerType'] == -5, 'NearestTime'].values[0]
    onsets = marker['NearestTime'] - t0
    durations= onsets.diff()
    ons = onsets[:-1].values
    dur = durations.drop([0]).values
    onsetdur = pd.DataFrame({'onset': ons, 'duration': dur, 'stimulus_type': range(6)})
    onsetdur = onsetdur.drop([0])
    return( onsetdur.round(3) )
    
def mirror_data(data):
    data[:, [1, 15]] = data[:, [15, 1]]
    data[:, [2, 16]] = data[:, [16, 2]]
    data[:, [3, 13]] = data[:, [13, 3]]
    data[:, [4, 14]] = data[:, [14, 4]]
    data[:, [5, 11]] = data[:, [11, 5]]
    data[:, [6, 12]] = data[:, [12, 6]]
    data[:, [7, 9]] = data[:, [9, 7]]
    data[:, [8, 10]] = data[:, [10, 8]]
    return( data )

def data2nii(data, output_filename):
    coords = np.loadtxt(template_txt_path)
    idx_out_coords = coords[:, 0:3]
    for i in range(data.shape[0]-1):
        idx_out_data = np.transpose(data)[1:17, i]
        idx_out = np.c_[idx_out_coords, idx_out_data]
        np.savetxt('t.txt', idx_out)
        os.system("3dUndump -master "+template_nii_path+" -ijk -datum float -prefix tmp"+str(i).zfill(6)+".nii t.txt > /dev/null 2>&1") 
        os.system("rm t.txt")

    os.system("find . -name 'tmp*.nii' | sort > tmp.filelist")
    os.system("3dTcat -tr "+str(tr)+" -prefix "+output_filename+" `cat tmp.filelist` > /dev/null 2>&1")
    os.system("3drefit -denote "+output_filename+" > /dev/null 2>&1")

    # cleanup
    os.system("rm tmp.filelist")
    os.system("rm tmp*.nii")

# \support functions
# ------------------

# -------------
# parsing input
# -------------

# -foldername
dirname = args.foldername
isExist = os.path.exists(dirname)
if isExist:
    os.chdir(dirname)
else: 
    print("input error, folder doesn't exist!")
    exit(-1)

# grabbing all *.txt files in specified foldername
# any available data should be in:
# <...>.Block#.txt and <...>.Time#.txt pairs
txt_files = []
for file in glob.glob("*.txt"):
    txt_files.append(file)

# -tag
valid_tags = ['hbo', 'hbt', 'hbr', 'oxy']
thetag = args.tag
if thetag in valid_tags:
    tag_files = find_pairs(thetag)
else:
    print("input error, invalid tag, choose between hbo, hbt, hbr, and oxy.")
    exit(-1)

prefix = args.prefix

if len(tag_files) > 1:
    print("number of acquisitions found: "+str(len(tag_files)))

for i in range(len(tag_files)):

    # attempt to locate tagged txt files in folder
    tag_data = parse_pair(tag_files[i][0], tag_files[i][1])
    tag_marker = parse_marker(tag_files[i][2], tag_data[:,0])

    # trim data to start and end of acquisition tags
    # this is not necessary if mem and gonogo tasks will be exported separately
    #ttag_data = trim_data(tag_data, tag_marker)
    #ttag_marker = trim_marker(tag_marker)

    # checking if the dataset is mirrored, and fix it
    # assuming -prefix option is IDENTICAL to "patient number tag"
    mirrored_list = open(mirrored_data_list, "r")
    mirrored_datasets = mirrored_list.read()

    if prefix in mirrored_datasets:
        print("STATUS: dataset in mirrored dataset list, reorganizing channels...")
        tag_data = mirror_data(tag_data)

    # -------------------------------------------
    # splitting data into go-nogo part and memory
    # -------------------------------------------

    # extract go-nogo part
    export_gonogo_flag = False
    try:
        data_gonogo = split_data(tag_data, tag_marker, 107, 108)
    except IndexError: 
        print("getting gonogo part of the data failed! skipping export")
        export_gonogo_flag = False
        
    # extract memory task data (include baseline, or resting-state)
    export_mem_flag = True
    try: 
        data_mem = split_data(tag_data, tag_marker, -5, 106)
    except IndexError:
        print("getting vs-mem part of the data failed! skipping export")
        export_mem_flag = False
    try:
        marker_mem = trim_marker(tag_marker, -5, 106)
    except IndexError:
        print("getting vs-mem part of the marker failed! skipping export")
        export_mem_flag = False

    # marker2onsetdur() assumes exactly one each of markers -5,101,102,103,
    # 104,105,106 (start, 5 block onsets, end). trim_marker()/drop_duplicates()
    # don't guarantee that - e.g. a spuriously repeated block-onset marker
    # alongside missing ones (both seen in this dataset) passes trim_marker
    # without raising, but leaves marker_mem the wrong shape and would corrupt
    # (or crash) onsetdur.txt. Catch that here instead of downstream.
    if export_mem_flag:
        expected_types = {-5, 101, 102, 103, 104, 105, 106}
        found_types = set(marker_mem['MarkerType'].values)
        if found_types != expected_types:
            missing = sorted(expected_types - found_types)
            unexpected = sorted(found_types - expected_types)
            print(f"vs-mem markers incomplete/malformed (missing={missing}, "
                  f"unexpected={unexpected})! skipping export")
            export_mem_flag = False


    # ------------
    # output files
    # ------------

    # constructing output filenames according to input argument: prefix and tag
    prefix = args.prefix
    gonogo_outfilename = prefix+'.'+thetag+'.gonogo.fnirs.nii'
    mem_outfilename = prefix+'.'+thetag+'.mem.fnirs.nii'
    mem_marker_outfilename = prefix+'.'+thetag+'.mem.timing.txt'
    mem_onsetdur_outfilename = prefix+'.'+thetag+'.mem.onsetdur.txt'

    # saving go-nogo
    if export_gonogo_flag:
        if args.skip_nii:
            print("skipping go-no-go .nii (skip_nii)")
        else:
            print("saving go-no-go")
            data2nii(data_gonogo, gonogo_outfilename)

    # saving mem data/marker pair
    # note: onsetdur output may be used as entry to AFNI's timing_tool.py with -fsl_timing_files flag:
    # i.e. timing_tool.py -fsl_timing_files <onsetdur.txt output from this script> -write_timing <outout>
    #np.savetxt(mem_marker_outfilename, marker_mem, fmt='%1.3f')
    if export_mem_flag:
        print("saving corsi")
        marker_mem.to_csv(mem_marker_outfilename, sep=' ', index=False)
        onsetdur_mem = marker2onsetdur(marker_mem)
        onsetdur_mem.to_csv(mem_onsetdur_outfilename, sep=' ', index=False, header=False )
        if args.skip_nii:
            print("skipping .nii for " + mem_outfilename + " (skip_nii)")
        else:
            print("saving data to: " + mem_outfilename)
            data2nii(data_mem, mem_outfilename)
