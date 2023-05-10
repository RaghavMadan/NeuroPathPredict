

import nilearn
from nilearn import image
import numpy as np
import pandas as pd
import math
import gc
import os
from os.path import join, exists, split

var = "MNI"
ROI = (["AM",25856],["HC",25752],["IPL",24000],["MFG",25232],["SMTG",25272])
Radii = [1,2,5,7,10,12.5,15,20,25,30,40]

homedir = os.getcwd()
print(homedir)
indir = join(homedir,"bz_output_data")
print(indir)
wdir = join(indir,var)
print(wdir)
fldrs = os.listdir(wdir)

def variable_extractor(fldrs,ROI,Radii):
    for fldr in fldrs:
        print(fldr)
        path = join(wdir,fldr)
        files = os.listdir(path)
        files.sort()
        print (f'total files to processed for {fldr} : {len(files)}')
        df_AM, df_HC, df_IPL, df_MFG, df_SMTG = covariate_extractor(files,ROI,Radii)

        out_fl_AM = os.join(wdir,f'bz_{fldr}_AM')
        df_AM.to_csv(out_fl_AM)

        out_fl_HC = os.join(wdir,f'bz_{fldr}_HC')
        df_HC.to_csv(out_fl_HC)

        out_fl_IPL = os.join(wdir,f'bz_{fldr}_IPL')
        df_IPL.to_csv(out_fl_IPL)

        out_fl_MFG = os.join(wdir,f'bz_{fldr}_MFG')
        df_MFG.to_csv(out_fl_MFG)

        out_fl_SMTG = os.join(wdir,f'bz_{fldr}_SMTG')
        df_SMTG.to_csv(out_fl_SMTG)

def covariate_extractor(files,ROI,Radii):
    index = 0
    for file in files:
        roi = file.split('_')[2]
        r = file.split('_')[3]
        df_AM, df_HC, df_IPL, df_MFG, df_SMTG = append(roi,r,file,ROI,Radii)

    return df_AM, df_HC, df_IPL, df_MFG, df_SMTG

def df_setup(ROI,Radii):
    col_names = []
    for rad in Radii:
            col_name = f'bz_rad_{rad}'
            col_names.append(col_name)
    print (col_names)

    for roi,rows in ROI:
        A = np.empty(rows,len(Radii))
        
        if roi == "AM":
            df_AM = pd.DataFrame(A,columns =col_names)

        if roi == "HC":
            df_HC = pd.DataFrame(A,columns =col_names)

        if roi == "IPL":
            df_IPL = pd.DataFrame(A,columns =col_names)

        if roi == "MFG":
            df_MFG = pd.DataFrame(A,columns =col_names)

        if roi == "SMTG":
            df_SMTG = pd.DataFrame(A,columns =col_names)

    return df_AM, df_HC, df_IPL, df_MFG, df_SMTG

def append(roi,r,file,ROI,Radii):
    dataset = pd.read_csv(file)
    df = pd.DataFrame(dataset)
    df = df[df.columns[1]]
    col_head = f'bz_rad_{r}'

    df_AM, df_HC, df_IPL, df_MFG, df_SMTG = df_setup(ROI,Radii)

    if roi == "AM":
        df_AM[col_head] = df.values

    if roi == "HC":
        df_HC[col_head] = df.values

    if roi == "IPL":
        df_IPL[col_head] = df.values

    if roi == "MFG":
        df_MFG[col_head] = df.values

    if roi == "SMTG":
        df_SMTG[col_head] = df.values
    
    return df_AM, df_HC, df_IPL, df_MFG, df_SMTG
        