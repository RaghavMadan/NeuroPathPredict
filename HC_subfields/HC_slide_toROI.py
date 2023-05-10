import nilearn
from nilearn import image
import numpy as np
import numpy.linalg as npl
import pandas as pd
import nibabel as nib
from nibabel.affines import apply_affine
import math
import gc
import warnings
import os
from os.path import join, exists, split
import sys
warnings.filterwarnings('ignore')

#Inputs

# input_dir
homedir = os.getcwd()
input_dir = join(homedir,'coord_ext/HC_toROI')
out_dir = join(homedir,'coord_ext')
files = os.listdir(input_dir)
files.sort()
print(files)

df_hc = pd.read_csv(join(homedir,'coord_ext/ROI_MNI_xyz_coord_HC_mask.csv'))
print(df_hc.shape)
print(df_hc.head)
df_hc_roi = df_hc
df_hc_roi['roi_label'] = ""


#Batch input files
print("total files to be extracted", len(files))
counter = 0
for file in files:
    name = file.split('.csv')[0]
    name = name.split('ROI_MNI_xyz_coord_')[1]
    input_file = join(input_dir,file)
    print(name)
    df_roi = pd.read_csv(join(input_file))
    df_roi = df_roi.drop(df_roi.columns[0],axis=1)
    #print(df_roi.shape)
    df_roi['roi_label'] = name

    df_hc_roi = pd.merge(df_hc_roi, df_roi, on = ['0','1','2'], how = 'outer')

    counter = counter + 1
    #print(counter)

df_hc_roi = df_hc_roi.drop(df_hc_roi.columns[0],axis=1)



print(df_hc_roi.head)
df_hc_roi.to_csv(join(out_dir,'ROI_HC_xyz_coord_labels.csv'))
