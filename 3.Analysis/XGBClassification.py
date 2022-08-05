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
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error as MSE
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.metrics import plot_confusion_matrix


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
# USE XGBOOST - CLASSIFICATION
#______________________________________________________________________________
Nout = 4
xgb_ps = {
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
for Elevel in ["100max", "200max"]:#,  "200min"]: 
    
    data = pd.read_csv("../myData/TSTree_SimFromCP"+Elevel+".txt", sep = ";")
    data = data.apply(pd.to_numeric, errors='coerce')
    x = data.iloc[:, 2:-Nout]
    y = data.iloc[:, -2]
    
    yclass = pd.Series(name = "N Photons", index=np.arange(len(data)), 
                       dtype = int)
    for i in range(len(data)):
        yclass.loc[i] = 0 if y[i] == 1 else 1
    count0 = 0
    for yi in y:
        if yi < 1: count0 += 1
    count2 = len(y) - count0
    
    train_x, test_x, train_y, test_y = train_test_split(x, yclass,
                                                        test_size = 0.2, 
                                                        random_state = 123)
    # start = timer()
    # xgb_c = xgb.XGBClassifier(**xgb_ps, use_label_encoder = False)
    # xgb_c.fit(train_x, train_y)
    # predsclass = xgb_c.predict(test_x)
    # end = timer()
    # print(f"\n{count0} events with ONE photon\n{count2} with TWO photons")
    # print("Time taken: {:.9f}s".format(end-start)); del start, end
    # print("Accuracy (classification): %.4f" %accuracy_score(test_y, predsclass))
    # xgb.plot_importance(xgb_c, importance_type = "weight", max_num_features= 25,
    #                     title = "Weight Importance", show_values = False)
    # pl.tight_layout()
    # xgb.plot_importance(xgb_c, importance_type = "gain", max_num_features= 25,
    #                     title = "Gain Importance", show_values = False)
    # pl.tight_layout()
    
    
    print("\n****************************************************************")
    print("NOW BALANCE ******************************************************")
    data = data.sort_values("DeltaCPVtxXY", ascending = False)
    count0 = 0
    for yi in y:
        if yi == 1: count0 += 1
    count2 = len(y) - count0
    diff = 2*count0 - len(y)
    if diff > 0:
        x = x.drop(x.index[-diff:])
        y = y.drop(y.index[-diff:])
        yclass = yclass.drop(yclass.index[-diff:])
    else:
        x = x.drop(x.index[:-diff])
        y = y.drop(y.index[:-diff])
        yclass = yclass.drop(yclass.index[:-diff])
    y = y.reset_index().iloc[:, 1]
    
    count0 = 0
    for yi in y:
        if yi == 1: count0 += 1
    count2 = len(y) - count0
    train_x, test_x, train_y, test_y = train_test_split(x, yclass,
                                                        test_size = 0.2, 
                                                        random_state = 123)
    start = timer()
    xgb_c = xgb.XGBClassifier(**xgb_ps, use_label_encoder = False)
    xgb_c.fit(train_x, train_y)
    predsclass = xgb_c.predict(test_x)
    end = timer()
    print(f"\n{count0} events with ONE photon\n{count2} with TWO photons")
    print("Time taken: {:.9f}s".format(end-start)); del start, end
    print("Accuracy (classification): %.4f" %accuracy_score(test_y, predsclass))
    ConfusionMatrixDisplay(
        confusion_matrix = confusion_matrix(test_y, predsclass, normalize = "true"),
        display_labels = ["Single", "Double"]
        ).plot( cmap = "Greys", colorbar = False)
    xgb.plot_importance(xgb_c, importance_type = "weight", max_num_features= 25,
                        title = "Weight Importance", show_values = False)
    pl.tight_layout()
    xgb.plot_importance(xgb_c, importance_type = "gain", max_num_features= 25,
                        title = "Gain Importance", show_values = False)
    pl.tight_layout()
    
    
    print("\n****************************************************************")
    print("NOW REMOVE DUBBIOUS **********************************************")
    data = data.drop(["ProjACnvx", "ProjMaxDist", "ProjMaxDistCentre"], axis=1)
    x    = x   .drop(["ProjACnvx", "ProjMaxDist", "ProjMaxDistCentre"], axis=1)
    i = 1
    while i > 0:
        try:
            data = data.drop([f"AConvex{i}", f"MaxDist{i}"], axis = 1)
            x    = x   .drop([f"AConvex{i}", f"MaxDist{i}"], axis = 1)
            i += 1
        except KeyError:
            break
    i = 0
    while i > 0:
        try:
            data = data.drop([f"NLaysWNLCs{i}",f"NLaysWNLCssmaller{i}"],axis=1)
            x    = x   .drop([f"NLaysWNLCs{i}",f"NLaysWNLCssmaller{i}"],axis=1)
            i += 1
        except KeyError:
            break
    count0 = 0
    for yi in y:
        if yi < 1: count0 += 1
    count2 = len(y) - count0
    train_x, test_x, train_y, test_y = train_test_split(x, yclass,
                                                        test_size = 0.2, 
                                                        random_state = 123)
    start = timer()
    xgb_c = xgb.XGBClassifier(**xgb_ps, use_label_encoder = False)
    xgb_c.fit(train_x, train_y)
    predsclass = xgb_c.predict(test_x)
    end = timer()
    print(f"\n{count0} events with ONE photon\n{count2} with TWO photons")
    print("Time taken: {:.9f}s".format(end-start)); del start, end
    print("Accuracy (classification): %.4f" %accuracy_score(test_y, predsclass))
    ConfusionMatrixDisplay(
        confusion_matrix = confusion_matrix(test_y, predsclass, normalize = "true"),
        display_labels = ["Single", "Double"]
        ).plot( cmap = "Greys", colorbar = False)
    xgb.plot_importance(xgb_c, importance_type = "weight", max_num_features= 25,
                        title = "Weight Importance (no dubbious)", 
                        show_values = False)
    pl.tight_layout()
    xgb.plot_importance(xgb_c, importance_type = "gain", max_num_features= 25,
                        title = "Gain Importance (no dubbious)",
                        show_values = False)
    pl.tight_layout()
    
    
    print("\n****************************************************************")
    print("NOW REMOVE LCs ***************************************************")
    i = 1
    while i > 0:
        try:
            data = data.drop([f"E{i}", f"LN{i}"], axis = 1)
            x = x.drop([f"E{i}", f"LN{i}"], axis = 1)
            i += 1
        except KeyError:
            break
    count0 = 0
    for yi in y:
        if yi < 1: count0 += 1
    count2 = len(y) - count0
    train_x, test_x, train_y, test_y = train_test_split(x, yclass,
                                                        test_size = 0.2, 
                                                        random_state = 123)
    start = timer()
    xgb_c = xgb.XGBClassifier(**xgb_ps, use_label_encoder = False)
    xgb_c.fit(train_x, train_y)
    predsclass = xgb_c.predict(test_x)
    end = timer()
    print(f"\n{count0} events with ONE photon\n{count2} with TWO photons")
    print("Time taken: {:.9f}s".format(end-start)); del start, end
    print("Accuracy (classification): %.4f" %accuracy_score(test_y, predsclass))
    ConfusionMatrixDisplay(
        confusion_matrix = confusion_matrix(test_y, predsclass, normalize = "true"),
        display_labels = ["Single", "Double"]
        ).plot( cmap = "Greys", colorbar = False)
    xgb.plot_importance(xgb_c, importance_type = "weight", max_num_features= 25,
                        title = "Weight Importance (no lcs)",
                        show_values = False)
    pl.tight_layout()
    xgb.plot_importance(xgb_c, importance_type = "gain", max_num_features= 25,
                        title = "Gain Importance (no lcs)",
                        show_values = False)
    pl.tight_layout()


















