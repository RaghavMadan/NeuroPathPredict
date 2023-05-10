
import os
from os.path import join
import sys
import nilearn
from nilearn import image
from nilearn.masking import apply_mask
import numpy as np
import pandas as pd
import shutil
import gc

homedir = os.getcwd()
print(homedir)
indir = join(homedir,'processed/4_edt_m')
print(indir)
outdir = join(homedir,'EDT_output_data')
mask_dir = join(homedir,'MNI_ROI_masks')
ROI = ("AM","HC","IPL","MFG","SMTG")

files = os.listdir(indir)
print(f'number of files to be processed: {len(files)}')
files.sort()
print(files)
counter = 0

for file in files:
    
    infile = join(indir,file)
    name = file.split('_masked.nii.gz')[0]
    # a = name.split('mni_icbm152_CerebrA_tal_nlin_sym_09c_')[0]
    # b = name.split('mni_icbm152_CerebrA_tal_nlin_sym_09c_')[1]
    # name = f'{a}CerebrA_{b}'
    img=nilearn.image.load_img(infile, wildcards=True, dtype=None)
    data_obj = img.get_fdata()
    print(name)
    print(data_obj.shape)
    for roi in ROI:
        fl_nm = f'edt_{name}_{roi}.csv'
        #print(fl_nm)
        out_fldr = join(outdir,roi)
        out_fl = join(out_fldr,fl_nm)
        mask = join(mask_dir,f'{roi}_mask.nii.gz')

        masked_data = apply_mask(img,mask,dtype='f',smoothing_fwhm=None)
        print(masked_data.shape)

        df = pd.DataFrame(masked_data)

        df.to_csv(out_fl)
        #shutil.move(fl_path,out_path)
        del fl_nm, out_fldr, out_fl, mask, masked_data,df
    
    del name, img, data_obj
    gc.collect()