import os
from os.path import join, exists, split
import shutil

homedir = os.getcwd()
print(homedir)
indir = join(homedir,"bz_output_data/consolidated")
outdir = join(homedir,"bz_output_data/bz_output_roi")

ROI = ("AM","HC","IPL","MFG","SMTG")

files = os.listdir(indir)
print(f'Number of files to be moved: {len(files)}')
counter = 0

for file in files:
    
    fl_path = join(indir,file)
    
    for roi in ROI:
        if roi in file:
            out_fldr = f'{roi}'
            out_path = join(outdir,out_fldr)
            shutil.move(fl_path,out_path)
    counter = counter+1
    print(f'files moved: {counter}')
        
print ('Move complete')

