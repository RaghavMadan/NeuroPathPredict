
import nilearn
from nilearn import image
import numpy as np
import pandas as pd
import math
import gc
import os
from os.path import join, exists, split

var = "WMT"
ROI = (["AM",25856],["HC",25752],["IPL",24000],["MFG",25232],["SMTG",25272])
Radii = [1,2,5,7,10,12.5,15,20,25,30,40]

homedir = os.getcwd()
print(homedir)
indir = join(homedir,"bz_output_data")
print(indir)
wdir = join(indir,var)
print(wdir)
fldrs = os.listdir(wdir)

def append_val(roi,r, file, ROI, Radii, path, df_AM, df_HC, df_IPL, df_MFG, df_SMTG):
    fl_nm = join(path,file)
    dataset = pd.read_csv(fl_nm)
    df = pd.DataFrame(dataset)
    df = df[df.columns[1]]
    col_head = f'bz_rad_{r}'

    if roi == "AM":
        df_AM[col_head] = df.values
    elif roi == "HC":
        df_HC[col_head] = df.values
    elif roi == "IPL":
        df_IPL[col_head] = df.values
    elif roi == "MFG":
        df_MFG[col_head] = df.values
    elif roi == "SMTG":
        df_SMTG[col_head] = df.values

    return df_AM, df_HC, df_IPL, df_MFG, df_SMTG

def df_setup(ROI,Radii):
    col_names = []
    for rad in Radii:
            col_name = f'bz_rad_{rad}'
            col_names.append(col_name)

    for roi,rows in ROI:
        A = np.empty([rows,len(Radii)])
        
        if roi == "AM":
            df_AM = pd.DataFrame(A,columns =col_names)

        elif roi == "HC":
            df_HC = pd.DataFrame(A,columns =col_names)

        elif roi == "IPL":
            df_IPL = pd.DataFrame(A,columns =col_names)

        elif roi == "MFG":
            df_MFG = pd.DataFrame(A,columns =col_names)

        elif roi == "SMTG":
            df_SMTG = pd.DataFrame(A,columns =col_names)

    return df_AM, df_HC, df_IPL, df_MFG, df_SMTG

def covariate_extractor(files,ROI,Radii,path,df_AM, df_HC, df_IPL, df_MFG, df_SMTG):
    for file in files:
        roi = file.split('_')[-2]
        r = file.split('_')[-1]
        r = r.split('.csv')[0]
        df_AM, df_HC, df_IPL, df_MFG, df_SMTG = append_val(roi,r,file,ROI,Radii,path, df_AM, df_HC, df_IPL, df_MFG, df_SMTG)

    return df_AM, df_HC, df_IPL, df_MFG, df_SMTG

def variable_extractor(fldrs,ROI,Radii,wdir):
    for fldr in fldrs:
        print(fldr)
        path = join(wdir,fldr)
        files = os.listdir(path)
        #files.sort()
        print (f'total files to processed for {fldr} : {len(files)}')
        df_AM, df_HC, df_IPL, df_MFG, df_SMTG = df_setup(ROI,Radii)
        df_AM, df_HC, df_IPL, df_MFG, df_SMTG = covariate_extractor(files,ROI,Radii,path, df_AM, df_HC, df_IPL, df_MFG, df_SMTG)
        out_fl_AM = join(wdir,f'bz_{fldr}_AM.csv')
        df_AM.to_csv(out_fl_AM)

        out_fl_HC = join(wdir,f'bz_{fldr}_HC.csv')
        df_HC.to_csv(out_fl_HC)

        out_fl_IPL = join(wdir,f'bz_{fldr}_IPL.csv')
        df_IPL.to_csv(out_fl_IPL)

        out_fl_MFG = join(wdir,f'bz_{fldr}_MFG.csv')
        df_MFG.to_csv(out_fl_MFG)

        out_fl_SMTG = join(wdir,f'bz_{fldr}_SMTG.csv')
        df_SMTG.to_csv(out_fl_SMTG)

        print(f'{fldr} is processed')
        del df_AM, df_HC, df_IPL, df_MFG, df_SMTG
        gc.collect()

variable_extractor(fldrs, ROI, Radii, wdir)