#!/usr/bin/python

import numpy as np
import glob
import os
import sys

# Navigate to the problematic folder
folder = "testdata/1008KP"
os.chdir(folder)

# Find marker file
marker_files = glob.glob("*.Marker*.txt")
print(f"Found marker files: {marker_files}")

if marker_files:
    marker_file = marker_files[0]
    print(f"\nReading: {marker_file}")
    
    # Load and display markers
    mrk = np.loadtxt(marker_file, delimiter='\t', skiprows=14)
    print(f"\nMarker shape: {mrk.shape}")
    print("\nUnique marker types found:")
    unique_markers = np.unique(mrk[:,1])
    print(unique_markers)
    
    print("\nAll markers with times:")
    for i in range(mrk.shape[0]):
        print(f"Time: {mrk[i,0]:.3f}, Type: {int(mrk[i,1])}")
