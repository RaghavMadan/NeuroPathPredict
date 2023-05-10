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

warnings.filterwarnings('ignore')

#Inputs

# input_dir
input_dir = "/Users/raghav_m/Desktop/NeuroPathPredict/input_files/nifti"
print(input_dir)

#Batch input files
data = []
files = os.listdir(input_dir)
files.sort()
print(files)
print("total files to be extracted", len(files))
counter = 0
for file in files:
    name = file.split('.nii.gz')[0]
    input_file = join(input_dir,file)
    print(input_file)
    img=nilearn.image.load_img(input_file, wildcards=True, dtype=None)
    data_obj = img.get_data()
    print(name)
    print(data_obj.shape)
    data.append(img)
    filename_txt = f'CFN_MNI_{name}.txt'
    key1 = f"CFN_{name}_int"

    n_i, n_j, n_k = data_obj.shape
    center_i = (n_i - 1) // 2  # // for integer division
    center_j = (n_j - 1) // 2
    center_k = (n_k - 1) // 2
    center_i, center_j, center_k
    center_vox_value = data_obj[center_i, center_j, center_k]
    print(center_vox_value, center_i, center_j, center_k)

    #Extracting coordinates and their values
    rows = n_i*n_j*n_k
    vox_extract = np.empty([rows, 4])
    index = 0;
    for i in range(n_i):
        for j in range(n_j):
            for k in range(n_k):
                vox_extract[index][0] = math.floor(i)
                vox_extract[index][1] = math.floor(j)
                vox_extract[index][2] = math.floor(k)
                vox_extract[index][3] = data_obj[i,j,k]
                index = index+1;
    df = pd.DataFrame(vox_extract)

    #Check df before saving to file
    df = df.astype({0:int,1:int,2:int})

    df = df.rename(columns={0:"X",1:"Y",2:"Z",3:key1})
    print(name)
    print(df)
    print("----min----")
    print(df[["X","Y","Z",key1]].min())
    print("----max----")
    print(df[["X","Y","Z",key1]].max())

    count = (df[key1] == 0).sum()
    print ('Vox with zero value:', count)

    #Saving data to a text file
    df.to_csv(filename_txt)

    del df
    del vox_extract
    del data_obj
    gc.collect()
    print("cache cleared")
    counter = counter + 1
    print(counter)