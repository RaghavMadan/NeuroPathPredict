
import numpy as np
import numpy.linalg as npl
import pandas as pd
import gc
import warnings
import os
from os.path import join, exists, split
import sys
import subprocess
warnings.filterwarnings('ignore')

# input/output dir
homedir = os.getcwd()
input_dir = join(homedir,'Vox')
out_dir = join(homedir,'roi_avg/BZ')
ROI = ("AM","HC","MFG","SMTG","IPL")

for roi in ROI:
    print(roi)

    in_fldr= join(input_dir,f'bz_output_data/bz_output_roi/{roi}')
    out_fl = join(out_dir,f'{roi}_roi_mean_val.csv')

    files =os.listdir(in_fldr)
    files.sort()
    #print(files)
    print(f'total files to be processed for {roi} are {len(files)}.')
    col = np.array([])
    val= np.array([])

    

    for file in files:
        name = file.split(f'_{roi}.csv')[0]

        if name[-1] == '_':
            name = name[:len(name)-1]

        df = pd.read_csv(join(in_fldr,file))
        bz_R = ('1','2','5','7','10','12.5','15','20','25','30','40')

        for bz_r in bz_R:
            col = col = np.hstack((col,f'{name}_bz_r{bz_r}'))
            mean_val = df[f"bz_rad_{bz_r}"].mean()
            val = np.hstack((val,mean_val))
        

    df = pd.concat([pd.DataFrame(col),pd.DataFrame(val)], axis =1)
    df.columns = ['var',f'{roi}']
    print(df)

    df.to_csv(out_fl)
    del df,col, val, mean_val, name, file, files
