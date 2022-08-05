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
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error as MSE
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.metrics import plot_confusion_matrix
from sklearn.preprocessing import StandardScaler


def is_correct(pred_i, y_i, thr = 0.2):
    if y_i == 1 or pred_i == 1:
        correct_E = pred_i - y_i == 0
    else:
        correct_E = abs(pred_i - y_i) / ( 1 - y_i) < thr
    return correct_E

def count_correct(pred, y, thr = 0.2):
    count = 0
    all_e = 0
    for (p_i, y_i) in zip(pred, y):
        count += int(is_correct(p_i, y_i, thr))
        all_e += 1
    return count / all_e

def count_correct_norm (pred, y, thr = 0.2):
    count = 0
    all_e = 0
    for (p_i, y_i) in zip(pred, y):
        if y_i == 1:
            count += int(is_correct(p_i, y_i, thr))
            all_e += 1
        else:
            count += int(is_correct(p_i, y_i, thr))
            count += int(is_correct(1-p_i, 1-y_i, thr))
            all_e += 2
    return count / all_e

#%%
Nout = 4
importances_c_all = [[], [], []]
importances_r_all = [[], [], []]

xgb_rps = {
          "learning_rate"    : 0.3, 
          "max_depth"        : 4,
          "n_estimators"     :25, 
          "objective"        : "reg:squarederror",
          "tree_method"      : "hist",
          "subsample"        : 1,
          "sampling_method"  : "gradient_based", 
          "colsample_bytree" : 1,
          "reg_lambda"       : 5, 
          "reg_alpha"        : 5,
         }
xgb_cps = {
          "eta"              : 0.3, 
          "max_depth"        : 4,
          "n_estimators"     : 25,
          "objective"        : "binary:logistic",
          "tree_method"      : "hist", 
          "subsample"        : 1,
          "sampling_method"  : "gradient_based", 
          "lambda"           : 5, 
          "alpha"            : 5, 
         }


for Ei, Elevel in enumerate(["100max", "200max", "200min"]):
    
    # DOWNLOAD DATA ___________________________________________________________
    data = pd.read_csv("../myData/TSTree_SimFromCP"+Elevel+".txt", sep = ";",
                       low_memory = False)
    data = data.apply(pd.to_numeric, errors='coerce')
    columns = data.columns
    
    #Remove NaNs
    i = 0
    while i < len(data):
        #if i%100000 == 0: print(i)
        row_i = data.loc[i]
        nans = np.isnan(row_i.values)
        if nans.any():
            data = data.drop(i)
            #print(i, f" : nan in {columns[nans==1]}")
        else:
            infs = np.isinf(row_i.values)
            if infs.any():
                data = data.drop(i)
                #print(i, f" : inf in {columns[infs==1]}")
            else:
                superbig = row_i.values > 10**38
                supersmall = row_i.values < - 10**38
                if superbig.any():
                    data.drop(i)
                    #print(i, f" : superbig in {columns[superbig==1]}")
                elif supersmall.any():
                    data.drop(i)
                    #print(i, f" : supersmall in {columns[supersmall==1]}")
        i += 1
    del i

    data = data.reset_index().iloc[:, 1:]
    
    # ARRANGE DATA ____________________________________________________________
    x = data.iloc[:, 2:-Nout]
    y = data.iloc[:, -2]
    x = x.reset_index().iloc[:, 1:]
    y = y.reset_index().iloc[:, 1]
    yclass = pd.Series(name="N Photons", index=np.arange(len(data)), dtype=int)
    for i in range(len(data)):
        try:
            yclass.loc[i] = 0 if y[i] == 1 else 1
        except KeyError:
            print("Continued")
            continue

    data = data.sort_values("DeltaCPVtxXY", ascending = False)
    count1 = len(y[y==1])
    diff = 2*count1- len(y)
    if diff > 0:
        x = x.drop(x.index[-diff:])
        y = y.drop(y.index[-diff:])
        yclass = yclass.drop(yclass.index[-diff:])
    else:
        x = x.drop(x.index[:-diff])
        y = y.drop(y.index[:-diff])
        yclass = yclass.drop(yclass.index[:-diff])
    x = x.reset_index().iloc[:, 1:]
    y = y.reset_index().iloc[:, 1]
    yclass = yclass.reset_index().iloc[:, 1]
    
    print("NOW REMOVE LCs ***************************************************")
    i = 1
    while i < 201:
        try:
            data = data.drop([f"E{i}", f"LN{i}"], axis = 1)
            x = x.drop([f"E{i}", f"LN{i}"], axis = 1)
            i += 1
        except KeyError:
            break
    
    count1 = len(y[y==1])
    count2 = len(y) - count1
    
    for j in tqdm([0]*20, "Repetitions with "+Elevel+f" & {[count1, count2]}"):
        seed = int(rand.random()*1000)
        trainc_x, testc_x, trainc_y, testc_y= train_test_split(x, yclass,
                                                              test_size = 0.25, 
                                                           random_state = seed)
        trainr_x, testr_x, trainr_y, testr_y= train_test_split(x, y,
                                                               test_size = 0.25, 
                                                           random_state = seed)
       
        xgb_c = xgb.XGBClassifier(**xgb_cps, use_label_encoder = False)
        xgb_r = xgb.XGBRegressor(**xgb_rps)
        xgb_c.fit(trainc_x, trainc_y)
        xgb_r.fit(trainr_x, trainr_y)

        predsclass = xgb_c.predict(testc_x)
        pred       = xgb_r.predict(testr_x)
        
        importances_c_all[Ei].append(xgb_c.feature_importances_)
        importances_r_all[Ei].append(xgb_r.feature_importances_)
        
del Ei, i, j
del x, y, yclass
del count1, count2
del seed, trainc_x, testc_x, trainc_y, testc_y
del trainr_x, trainr_y, testr_x, testr_y
del pred, predsclass

columns = data.columns[2:-4]

#%% INDIVIDUAL CASES AVERAGED FSCORES
nf = 25
importances_c = np.array(importances_c_all)
importances_r = np.array(importances_r_all)
for i in range(len(importances_c)):
    importances_c[i] = np.mean(importances_c_all[i], axis = 0)
    importances_r[i] = np.mean(importances_r_all[i], axis = 0)
    
    sort_c, columns_sort_c = zip(*sorted(zip(importances_c[i], columns))) 
    sort_r, columns_sort_r = zip(*sorted(zip(importances_r[i], columns))) 
    most_c = sort_c[-nf:]
    most_r = sort_r[-nf:]
    columns_most_c = columns_sort_c[-nf:]
    columns_most_r = columns_sort_r[-nf:]
    
    pl.figure()
    pl.title("Classifier Highest Weights (no lcs)")
    pl.xlabel("Average Importance (on 25 runs)")
    y_pos = np.arange(nf)
    pl.barh(y_pos, most_c, align = 'center')
    pl.yticks(y_pos, labels=columns_most_c)
    pl.grid()
    pl.tight_layout()
    pl.show()
    
    pl.figure()
    pl.title("Regressor Highest Weights (no lcs)")
    pl.xlabel("Average Importance (on 25 runs)")
    y_pos = np.arange(nf)
    pl.barh(y_pos, most_r, align = 'center')
    pl.yticks(y_pos, labels=columns_most_r)
    pl.grid()
    pl.tight_layout()
    pl.show()
    
  
#%% ALL TOGETHER
nf = 30
importances_sum_c = np.mean(importances_c, axis = 0)
importances_sum_r = np.mean(importances_r, axis = 0)
xgb_c.feature_importances = importances_sum_c
xgb_r.feature_importances = importances_sum_r

sort_c, columns_sort_c = zip(*sorted(zip(importances_sum_c, columns))) 
sort_r, columns_sort_r = zip(*sorted(zip(importances_sum_r, columns)))  

most_c = sort_c[-nf:]
most_r = sort_r[-nf:]
columns_most_c = columns_sort_c[-nf:]
columns_most_r = columns_sort_r[-nf:]
least_c = sort_c[:nf]
least_r = sort_r[:nf]
columns_least_c = columns_sort_c[:nf]
columns_least_r = columns_sort_r[:nf]

pl.figure()
pl.title("Classifier Highest Weights")
pl.xlabel("Average Importance (on 25 runs)")
y_pos = np.arange(nf)
pl.barh(y_pos, most_c, align = 'center')
pl.yticks(y_pos, labels=columns_most_c)
pl.grid()
pl.tight_layout()
pl.show()

pl.figure()
pl.title("Regressor Highest Weights")
pl.xlabel("Average Importance (on 25 runs)")
y_pos = np.arange(nf)
pl.barh(y_pos, most_r, align = 'center')
pl.yticks(y_pos, labels=columns_most_r)
pl.grid()
pl.tight_layout()
pl.show()

pl.figure()
pl.title("Classifier Lowest Weights")
pl.xlabel("Average Importance (on 25 runs)")
y_pos = np.arange(nf)
pl.barh(y_pos, least_c, align = 'center')
pl.yticks(y_pos, labels=columns_least_c)
pl.gca().invert_yaxis()
pl.grid()
pl.tight_layout()
pl.show()

pl.figure()
pl.title("Regressor Lowest Weights")
pl.xlabel("Average Importance (on 25 runs)")
y_pos = np.arange(nf)
pl.barh(y_pos, least_r, align = 'center')
pl.yticks(y_pos, labels=columns_least_r)
pl.gca().invert_yaxis()
pl.grid()
pl.tight_layout()
pl.show()


    # ConfusionMatrixDisplay(
    #     confusion_matrix = confusion_matrix(test_y, predsclass, normalize = "true"),
    #     display_labels = ["Single", "Double"]
    #     ).plot( cmap = "Greys", colorbar = False)
    # xgb.plot_importance(xgb_c, importance_type = "weight", max_num_features= 25,
    #                     title = "Weight Importance", show_values = False)
    # pl.tight_layout()
    # xgb.plot_importance(xgb_c, importance_type = "gain", max_num_features= 25,
    #                     title = "Gain Importance", show_values = False)
    # pl.tight_layout()
    