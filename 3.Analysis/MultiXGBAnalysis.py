# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 12:19:03 2022

@author: Stefano
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import random as rand

from tqdm import tqdm
from timeit import default_timer as timer
import cProfile, pstats 
import multiprocessing

from scipy.stats import chisquare, kstwo
from scipy.optimize import curve_fit, NonlinearConstraint
from iminuit import Minuit

from uncertainties import ufloat

ps = {"text.usetex": True,
      "font.size" : 16,
      "font.family" : "Times New Roman",
      "axes.labelsize": 16,
      "legend.fontsize": 14,
      "xtick.labelsize": 14,
      "ytick.labelsize": 14,
      "figure.figsize": [7.5, 6],
      "mathtext.default": "default"
       }
pl.rcParams.update(ps)
del ps

import xgboost as xgb
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error as MSE
from sklearn.multioutput import MultiOutputRegressor as MOR

from XGBMultiCV import *


#%% DOWNLOAD & ARRANGE DATA
#______________________________________________________________________________
Elevel = "200min"
Nout = 4
balance = 1
infolayers = 1

data = pd.read_csv("../myData/TSTree_SimFromCP" + Elevel + ".txt", sep = ";")
data = data.apply(pd.to_numeric, errors='coerce')

if balance == True:
    data = data.sort_values("DeltaCPVtxXY", ascending = False)
    #data.reindex()
    count0 = 0
    for i in data["DeltaCPVtxXY"]:
        if i == 0: count0 += 1
    diff = 2*count0 - len(data)
    if diff > 0:
        data = data.drop(data.index[-diff:])
    else:
        data = data.drop(data.index[:-diff])

if not infolayers:
    i = 1
    while i > 0:
        try:
            data = data.drop([f"E{i}", f"LN{i}"], axis = 1)
            i += 1
        except KeyError:
            break
data = data.iloc[:, 2:]

#%%
# SET UP MULTI XGBOOST WITH SKLEARN
#______________________________________________________________________________
xgb_ps = {
          "learning_rate"    : 0.3, #learning rate, small->high num_round
          "max_depth"        : 4,
          "n_estimators"     :100, #boost
          "objective"        : "reg:squarederror",
          "tree_method"      : "hist", #or "gpu-hist", fastest ones
          "subsample"        : 1,
          "sampling_method"  : "uniform", #gradient_based might be better
          "colsample_bytree" : 1,
          "reg_lambda"       : 5, #regularization term (could go up to 10)
          "reg_alpha"        : 5, #regularization term (could go up to 10)
          "n_jobs"           : 2
         }


myMOR, rmses, corrects = XGBMultiRegCV(data, 
                                     Nout = 4, 
                                     nfold = 5, 
                                     correct_thr = 10,
                                     printcv = False,
                                     **xgb_ps)

myMOR, rmses, corrects = XGBMultiRegCV(data, 
                                     Nout = 4, 
                                     nfold = 5, 
                                     correct_thr = 0.2,
                                     printcv = False,
                                     **xgb_ps)

myMOR, rmses, corrects = XGBMultiRegCV(data, 
                                     Nout = 4, 
                                     nfold = 5, 
                                     correct_thr = 0.1,
                                     printcv = False,
                                     **xgb_ps)

myMOR, rmses, corrects = XGBMultiRegCV(data, 
                                     Nout = 4, 
                                     nfold = 5, 
                                     correct_thr = 0.05,
                                     printcv = False,
                                     **xgb_ps)



# xgb_r = xgb.XGBRegressor(learning_rate    = 0.3, #small->high num_round
#                          max_depth        = 4,
#                          n_estimators     = 100, #boost
#                          subsample        = 1,
#                          sampling_method  = "uniform", #gradient_based might be better
#                          colsample_bytree = 1,
#                          reg_lambda       = 5, #could go to 10, regularises
#                          reg_alpha        = 5, #could go to 10, regularises
#                          seed = 123)
# MOR, rmses2, corrects2 = MultiCV(data, 4, xgb_r, 5)



