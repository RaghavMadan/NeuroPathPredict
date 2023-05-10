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
input_dir = "/Users/raghav_m/Desktop/NeuroPathPredict/NeuroPathPredict/EucDistTransform/output/2_bin"
print(f'input dir: {input_dir}')
out_dir_edt = "/Users/raghav_m/Desktop/NeuroPathPredict/NeuroPathPredict/EucDistTransform/output/3_edt"
out_dir_edt_m = "/Users/raghav_m/Desktop/NeuroPathPredict/NeuroPathPredict/EucDistTransform/output/4_edt_m"
mask_fl = "MNI152b_brain_mask_0.1.nii"

#Batch input files
data = []
files = os.listdir(input_dir)
files.sort()
print(files)
print("total files to be processed", len(files))
counter = 0

for file in files:
    name = file.split('.nii')[0]
    input_file = join(input_dir,file)
    print(input_file)

    out_fl_nm_edt = f'dxFrom_{name}.nii'
    out_fl_nm_edt_m = f'dxFrom_{name}_masked.nii'
    out_fl_edt = join(out_dir_edt,out_fl_nm_edt)
    out_fl_edt_m = join(out_dir_edt_m,out_fl_nm_edt_m)

    subprocess.run(["niimath", input_file, "-binv", "-edt", out_fl_edt])
    subprocess.run(["fslmaths", out_fl_edt, "-mas", mask_fl, out_fl_edt_m])
    counter = counter +1
    print(counter)

gc.collect()
print("cache cleared")