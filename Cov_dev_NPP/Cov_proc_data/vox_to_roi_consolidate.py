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
input_dir = join(homedir,'roi_avg')
out_dir = join(homedir,'roi_avg')
fldr = 'Var_ext'
in_fldr= join(input_dir,fldr)

files =os.listdir(in_fldr)
files.sort()
print(files)
print(len(files))

data = pd.DataFrame(columns=['var'])

for file in files:
    name = file.split('_roi_mean_val.csv')[0]

    df = pd.read_csv(join(in_fldr,file))
    df = df.drop(columns=['Unnamed: 0'])
    cols = [x for x in df.columns if x in df.columns or x == 'var']
    print(cols)
    data = pd.merge(data, df[cols], on='var', how='outer', suffixes=['',''])

print(data)
out_fl = (join(out_dir,f'{fldr}_roi_mean_val.csv'))
data.to_csv(out_fl)





