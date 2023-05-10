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
input_dir = join(homedir,'coord_ext/input')
out_dir = join(homedir,'coord_ext/output')
files = os.listdir(input_dir)
files.sort()
print(files)

#Batch input files
print("total files to be extracted", len(files))
counter = 0
for file in files:
    data = []
    name = file.split('_HC_masked_proc_0.5.nii.gz')[0]
    input_file = join(input_dir,file)
    #print(input_file)
    img=nilearn.image.load_img(input_file, wildcards=True, dtype=None)
    data_obj = img.get_data()
    print(name)
    print(data_obj.shape)
    filename_ijk = f'ROI_MNI_vox_coord_{name}.csv'
    filename_xyz = f'ROI_MNI_xyz_coord_{name}.csv'
    print(filename_ijk)
    print(filename_xyz)
    
    vox_center = (np.array(data_obj.shape) - 1) // 2
    print(vox_center)

    n_i, n_j, n_k = data_obj.shape
    vox_extract = np.empty([0,3])
    xyz_extract = np.empty([0,3])
    index = 0;
    for i in range(n_i):
        for j in range(n_j):
            for k in range(n_k):
                if data_obj[i,j,k] == 1:
                    vox_extract = np.append(vox_extract, np.array([[i,j,k]]), axis = 0)
                    xyz_extract = np.append(xyz_extract, (apply_affine(img.affine, np.array([[i,j,k]]))), axis =0)
                index = index+1;
    df_ijk = pd.DataFrame(vox_extract)
    df_xyz = pd.DataFrame(xyz_extract)
    print(df_ijk)
    print(df_xyz)
    
    df_ijk.to_csv(join(out_dir,filename_ijk))
    df_xyz.to_csv(join(out_dir,filename_xyz))
    
    del df_ijk
    del df_xyz
    del vox_extract
    del xyz_extract
    del data_obj
    gc.collect()
    print("cache cleared")
    counter = counter + 1
    print(counter)