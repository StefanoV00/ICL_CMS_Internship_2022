# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 14:57:55 2022

@author: Stefano
"""
import pandas as pd
import numpy as np

def mypd_shuffle(data):
    idx = np.random.permutation(data.index)
    data = data.reindex(idx)
    return data

def mypd_removenans(data):
    #Remove NaNs
    i = 0
    while i < len(data):
        row_i = data.loc[i]
        if row_i.isnull().values.any():
            data = data.drop(i)
        i += 1
    return data

path = "/vols/cms/magnan/HGCAL/Prod_Stefano/SinglePhotons2/"
name1 = "step3ticl_En40to400_eta21_run"
name2 = "_FlatTracksters"
ext = ".root"
allcomand = 'scp -rp sv519@lx02.hep.ph.ic.ac.uk:"'
for i in range(425, 501):
    allcomand += path+name1+f"{i}"+name2+ext+" "
allcomand += '" C:/root_v6.26.02/root-on-vscode/myTTrees/RemoteFiles'
print(allcomand)