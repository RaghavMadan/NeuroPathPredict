import nilearn
from nilearn import image
import numpy as np
import pandas as pd
import math
import gc
import warnings
import os
from os.path import join, exists, split
import sys
import subprocess

warnings.filterwarnings('ignore')

# input_dir and output_dir
input_dir = "/Users/raghav_m/Desktop/NeuroPathPredict/NeuroPathPredict/EucDistTransform/input"
print(f'input dir: {input_dir}')
out_dir = "/Users/raghav_m/Desktop/NeuroPathPredict/NeuroPathPredict/EucDistTransform/output/1_thres"
print(f'output dir: {out_dir}')
mv_pr_dir = "/Users/raghav_m/Desktop/NeuroPathPredict/NeuroPathPredict/EucDistTransform/processed"
# in_dir_bin = out_dir
# out_dir_bin = "/Users/raghav_m/Desktop/NeuroPathPredict/NeuroPathPredict/EucDistTransform/output/2_bin"
# in_dir_edt = out_dir_bin
# out_dir_edt = "/Users/raghav_m/Desktop/NeuroPathPredict/NeuroPathPredict/EucDistTransform/output/3_edt"
# in_dir_edt_m = out_dir_edt
# out_dir_edt_m = "/Users/raghav_m/Desktop/NeuroPathPredict/NeuroPathPredict/EucDistTransform/output/4_edt_m"

# Assign threshold
thr = '0'

#Batch input files
data = []
files = os.listdir(input_dir)
files.sort()
print(files)
print("total files to be processed", len(files))
counter = 0

for file in files:
    name = file.split('_proc_0.5.nii')[0]
    input_file = join(input_dir,file)
    print(input_file)

    out_fl_nm = f'{name}_thr_{thr}.nii'
    out_fl = join(out_dir,out_fl_nm)

    mv_pr_fl = join(mv_pr_dir,file)

    subprocess.run(["fslmaths", input_file, "-thr", thr, out_fl])
    subprocess.run(["mv", input_file, mv_pr_fl])
    counter = counter +1
    print(counter)


gc.collect()
print("cache cleared")