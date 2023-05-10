import os
from os.path import join, exists, split

homedir = join(split(os.getcwd())[0],split(os.getcwd())[1])
print(homedir)

def WMT(file):
    src = join(datadir,file)
    name = file.split('_proc_0.5.nii.gz')[0]
    dst = join(datadir,f'{name}.nii.gz')
    os.rename(src,dst)

def CerebrA(file):
    src = join(datadir,file)
    name = file.split('mni_icbm152_CerebrA_tal_nlin_sym_09c_')[-1]
    dst = join(datadir,f'CerebrA_{name}')
    os.rename(src,dst)

datadir = join(homedir,'var')
files = os.listdir(datadir)
print(len(files))
files.sort()

for file in files:
    #CerebrA(file)
    WMT(file)

del files