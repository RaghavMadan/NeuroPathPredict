#Script to upsample nifti files to 2009b

import os
from os.path import join, exists, split
import numpy as np
import sys
import subprocess
import warnings
import scipy
warnings.filterwarnings('ignore')
import nibabel as nib
import nibabel
from nibabel import processing
import matplotlib.pyplot as plt
import logging
import shutil
import gc
import nilearn
from nilearn import image

#Input directory
homedir = split(os.getcwd())[0]
print(homedir)
# create directory for saving data
datadir = join(homedir,'upsample/input_data')
outputdir = join(homedir,'upsample/output_data')

out_fld_1 = join(outputdir, '1_RegToMNI2009b_upsampled')
out_fld_2 = join(outputdir, '2_reg_upsampled_masked')

#Destination to move input files once processed
dst_path = join(homedir, 'upsample/proc_data')

#Loading ref image (MNI 2009b 0.5mm)
img_ref = 'mni_icbm152_t1_tal_nlin_sym_09b_hires.nii'

#Loading mask image (MNI 2009b 0.5mm)
nb_img_ref = nib.load(img_ref)
img_mask = 'MNI152b_brain_mask_0.1.nii'
nb_img_mask = nib.load(img_mask)
nb_img_mask_data = nb_img_mask.get_data()
nb_img_mask_data = nb_img_mask_data.astype(bool)

#Function to process the atlas
def proc_atlas(file):
    img_toreg = join(datadir,file)
    src_path = img_toreg
    logging.info(img_toreg)
    img_tr = nib.load(img_toreg)
    nb_img_toreg = img_tr.slicer[...,0]
    print(nb_img_toreg.shape)
    name = file.split('.nii')[0]
    print(f'starting process for {name}')
    out_fl_1 = join(out_fld_1, f'{name}_reg_up.nii')
    out_fl_2 = join(out_fld_2, f'{name}_proc_0.5.nii.gz')
    print(f'MNI2009b_shape:{(nb_img_ref).shape}')
    print(f'Atlas_shape:{(nb_img_toreg).shape}')
    
    
    out_img = nibabel.processing.resample_from_to(nb_img_toreg, nb_img_ref, order=1, mode='nearest')
    nib.save(out_img, out_fl_1)
    out_img_data = out_img.get_data()
    out_img_masked = out_img_data[nb_img_mask_data]
    unmasked_data = np.zeros(out_img_data.shape, dtype=out_img_masked.dtype)
    unmasked_data[nb_img_mask_data] = out_img_masked
    plt.matshow(unmasked_data[:, :, 100], cmap='bone', vmin='0', vmax='111')
    out_img_header = nb_img_ref.header
    out_img_header.set_data_dtype(out_img_masked.dtype)
    out_img_affine = nb_img_ref.affine
    out_img2 = nib.Nifti1Image(unmasked_data, out_img_affine, out_img_header)
    nib.save(out_img2, out_fl_2)
    logging.info(nb_img_ref.affine == out_img2.affine)
    logging.info(f'shape of input atlas: {nb_img_toreg.shape}')
    logging.info(f'Affine of input atlas: {nb_img_toreg.affine}')
    logging.info(f'shape of processed atlas: {out_img2.shape}')
    shutil.move(src_path, dst_path)
    print(f'{name} is processed')
    del img_toreg, src_path, out_img, out_img2
    gc.collect()

#Generating log file
log_fl = join(out_fld_2, "log.txt")
logging.basicConfig(filename=log_fl, level=logging.INFO,
                    format="%(asctime)s %(message)s")
                    
#Looping through all atlases
files = os.listdir(datadir)
files.sort()
print(files)
for file in files:
    proc_atlas(file)
    

print("done")
