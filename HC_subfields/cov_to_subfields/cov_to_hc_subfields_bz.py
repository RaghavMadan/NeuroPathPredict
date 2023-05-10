
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
input_dir = join(homedir,'Input')
out_dir = join(homedir,'Output')
proc_dir = join(homedir,'Proc')
files = os.listdir(input_dir)
files.sort()
print(files)

df_hc = pd.read_csv(join(homedir,'ROI_HC_subfield_labels.csv'))
df_hc = np.array(df_hc)
print(df_hc.shape)
n_hc_i,n_hc_j = df_hc.shape

#Batch input files
print("total files to be processed", len(files))
counter = 0

for file in files:
    
    input_file = join(input_dir,file)
    name = file.split('HC.csv')[0]
    print (name)
    df_cov = pd.read_csv(join(input_file))
    df_cov = np.array(df_cov)
    print(df_cov.shape)
    print(df_cov[0][1])
    n_cov_i, n_cov_j = df_cov.shape
    assert n_hc_i == n_cov_i , "Number of rows not equal in the two files."

    df_CA1_lh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_CA1_rh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_CA3_lh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_CA3_rh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_CA4_lh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_CA4_rh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_DG_lh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_DG_rh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_SB_lh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_SB_rh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_Pre_SB_lh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_Pre_SB_rh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_Para_SB_lh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])
    df_Para_SB_rh = pd.DataFrame(columns=['index', 'r1','r2','r5','r7','r10','r12.5','r15','r20','r25','r30','r40','x','y','z'])

    for i in range(n_hc_i):
        ind = df_cov[i][0]
        r1 = df_cov[i][1]
        r2 = df_cov[i][2]
        r5 = df_cov[i][3]
        r7 = df_cov[i][4]
        r10 = df_cov[i][5]
        r12_5 = df_cov[i][6]
        r15 = df_cov[i][7]
        r20 = df_cov[i][8]
        r25 = df_cov[i][9]
        r30 = df_cov[i][10]
        r40 = df_cov[i][11]
        x = df_hc[i][1]
        y = df_hc[i][2]
        z = df_hc[i][3]

        if df_hc[i][4] == "lh_HC_CA1":
            df_CA1_lh = df_CA1_lh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)

        elif df_hc[i][4] == "rh_HC_CA1":
            df_CA1_rh = df_CA1_rh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)

        elif df_hc[i][4] == "lh_HC_CA3":
            df_CA3_lh = df_CA3_lh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)

        elif df_hc[i][4] == "rh_HC_CA3":
            df_CA3_rh = df_CA3_rh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)

        elif df_hc[i][4] == "lh_HC_CA4":
            df_CA4_lh = df_CA4_lh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)

        elif df_hc[i][4] == "rh_HC_CA4":
            df_CA4_rh = df_CA4_rh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)

        elif df_hc[i][4] == "lh_HC_DG":
            df_DG_lh = df_DG_lh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)

        elif df_hc[i][4] == "rh_HC_DG":
            df_DG_rh = df_DG_rh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)    

        elif df_hc[i][4] == "lh_HC_SB":
            df_SB_lh = df_SB_lh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)  

        elif df_hc[i][4] == "rh_HC_SB":
            df_SB_rh = df_SB_rh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)       

        elif df_hc[i][4] == "lh_HC_pre_SB":
            df_Pre_SB_lh = df_Pre_SB_lh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)   

        elif df_hc[i][4] == "rh_HC_pre_SB":
            df_Pre_SB_rh = df_Pre_SB_rh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)    

        elif df_hc[i][4] == "lh_HC_para_SB":
            df_Para_SB_lh = df_Para_SB_lh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)  

        elif df_hc[i][4] == "rh_HC_para_SB":
            df_Para_SB_rh = df_Para_SB_rh.append({'index':ind,'r1':r1,'r2':r2,'r5':r5,'r7':r7,'r10':r10,'r12.5':r12_5,'r15':r15,
            'r20':r20,'r25':r25,'r30':r30,'r40':r40,'x':x,'y':y,'z':z},ignore_index = True)     


    df_CA1_lh.to_csv(join(out_dir,f"CA1_lh/{name}CA1_lh.csv"))
    df_CA1_rh.to_csv(join(out_dir,f"CA1_rh/{name}CA1_rh.csv"))
    df_CA3_lh.to_csv(join(out_dir,f"CA3_lh/{name}CA3_lh.csv"))
    df_CA3_rh.to_csv(join(out_dir,f"CA3_rh/{name}CA3_rh.csv"))
    df_CA4_lh.to_csv(join(out_dir,f"CA4_lh/{name}CA4_lh.csv"))
    df_CA4_rh.to_csv(join(out_dir,f"CA4_rh/{name}CA4_rh.csv"))
    df_DG_lh.to_csv(join(out_dir,f"DG_lh/{name}DG_lh.csv"))
    df_DG_rh.to_csv(join(out_dir,f"DG_rh/{name}DG_rh.csv"))
    df_SB_lh.to_csv(join(out_dir,f"SB_lh/{name}SB_lh.csv"))
    df_SB_rh.to_csv(join(out_dir,f"SB_rh/{name}SB_rh.csv"))
    df_Pre_SB_lh.to_csv(join(out_dir,f"Pre_SB_lh/{name}Pre_SB_lh.csv"))
    df_Pre_SB_rh.to_csv(join(out_dir,f"Pre_SB_rh/{name}Pre_SB_rh.csv"))
    df_Para_SB_lh.to_csv(join(out_dir,f"Para_SB_lh/{name}Para_SB_lh.csv"))
    df_Para_SB_rh.to_csv(join(out_dir,f"Para_SB_rh/{name}Para_SB_rh.csv"))

    proc_file = join(proc_dir,file)
    subprocess.run(["mv", input_file, proc_file])

    del df_CA1_lh, df_CA1_rh,df_CA3_lh,df_CA3_rh,df_CA4_lh,df_CA4_rh,i, df_DG_lh,df_DG_rh,df_SB_lh,df_SB_rh,df_Para_SB_lh,df_Para_SB_rh,df_Pre_SB_lh,df_Pre_SB_rh
    del ind,x,y,z,r1,r2,r5,r7,r10,r12_5,r15,r20,r25,r30,r40
    del name, df_cov

    gc.collect()
    print("cache cleared")
    counter = counter + 1
    print(counter)
