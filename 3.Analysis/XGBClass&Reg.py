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
importances_c = []
importances_r = []

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


for Elevel in ["100max", "200max", "200min"]:
    
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
    # yclass.loc[y==1] = 0
    # yclass.loc[y==0] = 1
    for i in range(len(data)):
        try:
            yclass.loc[i] = 0 if y[i] == 1 else 1
        except KeyError:
            print("Continued")
            continue
    
    
    print("\n****************************************************************")
    print("NOW BALANCE ******************************************************")
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
    
    count1 = len(y[y==1])
    count2 = len(y) - count1
    seed = int(rand.random()*1000)
    trainc_x, testc_x, trainc_y, testc_y= train_test_split(x, yclass,
                                                           test_size = 0.25, 
                                                           random_state = seed)
    trainr_x, testr_x, trainr_y, testr_y= train_test_split(x, y,
                                                           test_size = 0.25, 
                                                           random_state = seed)
    
    # print(f"\n{count1} events with ONE photon\n{count2} with TWO photons")
    # start = timer()
    # xgb_c = xgb.XGBClassifier(**xgb_cps, use_label_encoder = False)
    # xgb_r = xgb.XGBRegressor(**xgb_rps)
    # xgb_c.fit(trainc_x, trainc_y)
    # xgb_r.fit(trainr_x, trainr_y)
    # end = timer()
    # print("Training (75%) time: {:.9f}s".format(end-start))
   
    # start = timer()
    # predsclass = xgb_c.predict(testc_x)
    # end = timer()
    # print("Classifier (25%) time: {:.9f}s".format(end-start))
   
    # start = timer()
    # pred       = xgb_r.predict(testr_x)
    # end = timer()
    # print("Regressor (25%) time: {:.9f}s".format(end-start))
    # del start, end
    
    # conf_matrix= confusion_matrix(testc_y, predsclass, normalize = "true")
    # print("\nAccuracy (single class): %.4f" %conf_matrix[0,0])
    # print("Accuracy (double class): %.4f" %conf_matrix[1,1])
    # pred[predsclass == 0] = 1
    # #correctsph = count_correct(pred, testr_y, 10.00)
    # corrects20 = count_correct(pred, testr_y)
    # corrects10 = count_correct(pred, testr_y, 0.10)
    # corrects05 = count_correct(pred, testr_y, 0.05)
    # corrects20norm = count_correct_norm(pred, testr_y)
    # corrects10norm = count_correct_norm(pred, testr_y, 0.10)
    # corrects05norm = count_correct_norm(pred, testr_y, 0.05)
    # corrects20mono = count_correct(pred[testr_y == 1], testr_y[testr_y == 1])
    # corrects10mono = count_correct(pred[testr_y == 1], testr_y[testr_y == 1], 0.10)
    # corrects05mono = count_correct(pred[testr_y == 1], testr_y[testr_y == 1], 0.05)
    # corrects20bi = count_correct(pred[testr_y < 1], testr_y[testr_y < 1])
    # corrects10bi = count_correct(pred[testr_y < 1], testr_y[testr_y < 1], 0.10)
    # corrects05bi = count_correct(pred[testr_y < 1], testr_y[testr_y < 1], 0.05)
    # compare = abs(pred - np.array(testr_y)) / np.array(testr_y)
    # compare_multi = compare[testr_y != 1]
    # compare_mono  = compare[testr_y == 1]
    # #print("\n\nAccuracy (Nph) : %.4f" %(correctsph))
    # print("\nEvent  Accuracy (0.20): %.4f" %(corrects20))
    # print("Photon Accuracy (0.20): %.4f" %(corrects20norm))
    # print("One-Ph Accuracy (0.20): %.4f" %(corrects20mono))
    # print("Two-Ph Accuracy (0.20): %.4f" %(corrects20bi))
    # print("\nEvent  Accuracy (0.10): %.4f" %(corrects10))
    # print("Photon Accuracy (0.10): %.4f" %(corrects20norm))
    # print("One-Ph Accuracy (0.10): %.4f" %(corrects10mono))
    # print("Two-Ph Accuracy (0.10): %.4f" %(corrects10bi))
    # print("\nEvent  Accuracy (0.05): %.4f" %(corrects05))
    # print("Photon Accuracy (0.05): %.4f" %(corrects20norm))
    # print("One-Ph Accuracy (0.05): %.4f" %(corrects05mono))
    # print("Two-Ph Accuracy (0.05): %.4f" %(corrects05bi))
    # xgb.plot_importance(xgb_c, importance_type = "weight", max_num_features= 25,
    #                     title = "Classifier Weight",
    #                     show_values = False)
    # pl.tight_layout()
    # xgb.plot_importance(xgb_r, importance_type = "weight", max_num_features= 25,
    #                     title = "Regressor Weight", 
    #                     show_values = False)
    # pl.tight_layout()
    
    
    print("\n****************************************************************")
    print("NOW REMOVE LCs ***************************************************")
    i = 1
    while i < 201:
        try:
            data = data.drop([f"E{i}", f"LN{i}"], axis = 1)
            trainc_x = trainc_x.drop([f"E{i}", f"LN{i}"], axis = 1)
            trainr_x = trainr_x.drop([f"E{i}", f"LN{i}"], axis = 1)
            testc_x  = testc_x.drop([f"E{i}", f"LN{i}"], axis = 1)
            testr_x  = testr_x.drop([f"E{i}", f"LN{i}"], axis = 1)
            i += 1
        except KeyError:
            break
    print(f"\n{count1} events with ONE photon\n{count2} with TWO photons")
    start = timer()
    xgb_c = xgb.XGBClassifier(**xgb_cps, use_label_encoder = False)
    xgb_r = xgb.XGBRegressor(**xgb_rps)
    xgb_c.fit(trainc_x, trainc_y)
    xgb_r.fit(trainr_x, trainr_y)
    end = timer()
    print("Training (75%) time: {:.9f}s".format(end-start))
    
    start = timer()
    predsclass = xgb_c.predict(testc_x)
    end = timer()
    print("Classifier (25%) time: {:.9f}s".format(end-start))
    
    start = timer()
    pred       = xgb_r.predict(testr_x)
    end = timer()
    print("Regressor (25%) time: {:.9f}s".format(end-start))
    del start, end
    
    conf_matrix= confusion_matrix(testc_y, predsclass, normalize = "true")
    print("\nAccuracy (single class): %.4f" %conf_matrix[0,0])
    print("Accuracy (double class): %.4f" %conf_matrix[1,1])
    pred[predsclass == 0] = 1
    #correctsph = count_correct(pred, testr_y, 10.00)
    corrects20 = count_correct(pred, testr_y)
    corrects10 = count_correct(pred, testr_y, 0.10)
    corrects05 = count_correct(pred, testr_y, 0.05)
    corrects20norm = count_correct_norm(pred, testr_y)
    corrects10norm = count_correct_norm(pred, testr_y, 0.10)
    corrects05norm = count_correct_norm(pred, testr_y, 0.05)
    corrects20mono = count_correct(pred[testr_y == 1], testr_y[testr_y == 1])
    corrects10mono = count_correct(pred[testr_y == 1], testr_y[testr_y == 1], 0.10)
    corrects05mono = count_correct(pred[testr_y == 1], testr_y[testr_y == 1], 0.05)
    corrects20bi = count_correct(pred[testr_y < 1], testr_y[testr_y < 1])
    corrects10bi = count_correct(pred[testr_y < 1], testr_y[testr_y < 1], 0.10)
    corrects05bi = count_correct(pred[testr_y < 1], testr_y[testr_y < 1], 0.05)
    compare = abs(pred - np.array(testr_y)) / np.array(testr_y)
    compare_multi = compare[testr_y != 1]
    compare_mono  = compare[testr_y == 1]
    #print("\n\nAccuracy (Nph) : %.4f" %(correctsph))
    print("\nEvent  Accuracy (0.20): %.4f" %(corrects20))
    print("Photon Accuracy (0.20): %.4f" %(corrects20norm))
    print("One-Ph Accuracy (0.20): %.4f" %(corrects20mono))
    print("Two-Ph Accuracy (0.20): %.4f" %(corrects20bi))
    print("\nEvent  Accuracy (0.10): %.4f" %(corrects10))
    print("Photon Accuracy (0.10): %.4f" %(corrects20norm))
    print("One-Ph Accuracy (0.10): %.4f" %(corrects10mono))
    print("Two-Ph Accuracy (0.10): %.4f" %(corrects10bi))
    print("\nEvent  Accuracy (0.05): %.4f" %(corrects05))
    print("Photon Accuracy (0.05): %.4f" %(corrects20norm))
    print("One-Ph Accuracy (0.05): %.4f" %(corrects05mono))
    print("Two-Ph Accuracy (0.05): %.4f" %(corrects05bi))
    xgb.plot_importance(xgb_c, importance_type = "weight", max_num_features= 25,
                        title = "Classifier Weight (no lcs)",
                        show_values = False)
    pl.tight_layout()
    xgb.plot_importance(xgb_r, importance_type = "weight", max_num_features= 25,
                        title = "Regressor Weight (no lcs)", 
                        show_values = False)
    pl.tight_layout()
    importances_c.append(xgb_c.feature_importances_)
    importances_r.append(xgb_r.feature_importances_)
    continue
    
    print("\n****************************************************************")
    print("NOW REMOVE CORRELATED*********************************************")
    data = data.drop(["ProjRatio6890", "ProjRatio6895", "ProjRatio9095"],
                     axis=1)
    trainc_x= trainc_x.drop(["ProjRatio6890", "ProjRatio6895", "ProjRatio9095"],
                     axis=1)
    trainr_x= trainr_x.drop(["ProjRatio6890", "ProjRatio6895", "ProjRatio9095"],
                     axis=1)
    testc_x = testc_x .drop(["ProjRatio6890", "ProjRatio6895", "ProjRatio9095"],
                     axis=1)
    testr_x = testr_x .drop(["ProjRatio6890", "ProjRatio6895", "ProjRatio9095"],
                     axis=1)
    print(f"\n{count1} events with ONE photon\n{count2} with TWO photons")
    start = timer()
    xgb_c = xgb.XGBClassifier(**xgb_cps, use_label_encoder = False)
    xgb_r = xgb.XGBRegressor(**xgb_rps)
    xgb_c.fit(trainc_x, trainc_y)
    xgb_r.fit(trainr_x, trainr_y)
    end = timer()
    print("Training (75%) time: {:.9f}s".format(end-start))
    
    start = timer()
    predsclass = xgb_c.predict(testc_x)
    end = timer()
    print("Classifier (25%) time: {:.9f}s".format(end-start))
    
    start = timer()
    pred       = xgb_r.predict(testr_x)
    end = timer()
    print("Regressor (25%) time: {:.9f}s".format(end-start))
    del start, end
    
    conf_matrix= confusion_matrix(testc_y, predsclass, normalize = "true")
    print("\nAccuracy (single class): %.4f" %conf_matrix[0,0])
    print("Accuracy (double class): %.4f" %conf_matrix[1,1])
    pred[predsclass == 0] = 1
    #correctsph = count_correct(pred, testr_y, 10.00)
    corrects20 = count_correct(pred, testr_y)
    corrects10 = count_correct(pred, testr_y, 0.10)
    corrects05 = count_correct(pred, testr_y, 0.05)
    corrects20norm = count_correct_norm(pred, testr_y)
    corrects10norm = count_correct_norm(pred, testr_y, 0.10)
    corrects05norm = count_correct_norm(pred, testr_y, 0.05)
    corrects20mono = count_correct(pred[testr_y == 1], testr_y[testr_y == 1])
    corrects10mono = count_correct(pred[testr_y == 1], testr_y[testr_y == 1], 0.10)
    corrects05mono = count_correct(pred[testr_y == 1], testr_y[testr_y == 1], 0.05)
    corrects20bi = count_correct(pred[testr_y < 1], testr_y[testr_y < 1])
    corrects10bi = count_correct(pred[testr_y < 1], testr_y[testr_y < 1], 0.10)
    corrects05bi = count_correct(pred[testr_y < 1], testr_y[testr_y < 1], 0.05)
    compare = abs(pred - np.array(testr_y)) / np.array(testr_y)
    compare_multi = compare[testr_y != 1]
    compare_mono  = compare[testr_y == 1]
    #print("\n\nAccuracy (Nph) : %.4f" %(correctsph))
    print("\nEvent  Accuracy (0.20): %.4f" %(corrects20))
    print("Photon Accuracy (0.20): %.4f" %(corrects20norm))
    print("One-Ph Accuracy (0.20): %.4f" %(corrects20mono))
    print("Two-Ph Accuracy (0.20): %.4f" %(corrects20bi))
    print("\nEvent  Accuracy (0.10): %.4f" %(corrects10))
    print("Photon Accuracy (0.10): %.4f" %(corrects20norm))
    print("One-Ph Accuracy (0.10): %.4f" %(corrects10mono))
    print("Two-Ph Accuracy (0.10): %.4f" %(corrects10bi))
    print("\nEvent  Accuracy (0.05): %.4f" %(corrects05))
    print("Photon Accuracy (0.05): %.4f" %(corrects20norm))
    print("One-Ph Accuracy (0.05): %.4f" %(corrects05mono))
    print("Two-Ph Accuracy (0.05): %.4f" %(corrects05bi))
    
    
    print("\n****************************************************************")
    print("NOW REMOVE DUBBIOUS **********************************************")
    data = data.drop(["ProjACnvx", "ProjMaxDist", "ProjMaxDistCentre"], axis=1)
    trainc_x= trainc_x.drop(["ProjACnvx", "ProjMaxDist", "ProjMaxDistCentre"], 
                            axis=1)
    trainr_x= trainr_x.drop(["ProjACnvx", "ProjMaxDist", "ProjMaxDistCentre"], 
                            axis=1)
    testc_x = testc_x .drop(["ProjACnvx", "ProjMaxDist", "ProjMaxDistCentre"], 
                            axis=1)
    testr_x = testr_x .drop(["ProjACnvx", "ProjMaxDist", "ProjMaxDistCentre"], 
                            axis=1)
    i = 1
    while i > 0:
        try:
            data = data.drop([f"AConvex{i}", f"MaxDist{i}"], axis = 1)
            trainc_x= trainc_x.drop([f"AConvex{i}", f"MaxDist{i}"], axis = 1)
            trainr_x= trainr_x.drop([f"AConvex{i}", f"MaxDist{i}"], axis = 1)
            testc_x = testc_x .drop([f"AConvex{i}", f"MaxDist{i}"], axis = 1)
            testr_x = testr_x .drop([f"AConvex{i}", f"MaxDist{i}"], axis = 1)
            i += 1
        except KeyError:
            break
    i = 0
    while i > 0:
        try:
            data = data.drop([f"NLaysWNLCs{i}",f"NLaysWNLCssmaller{i}"],axis=1)
            trainc_x= trainc_x.drop([f"NLaysWNLCs{i}",f"NLaysWNLCssmaller{i}"],
                                    axis=1)
            trainr_x= trainr_x.drop([f"NLaysWNLCs{i}",f"NLaysWNLCssmaller{i}"],
                                    axis=1)
            testc_x = testc_x.drop([f"NLaysWNLCs{i}",f"NLaysWNLCssmaller{i}"],
                                    axis=1)
            testr_x = testr_x.drop([f"NLaysWNLCs{i}",f"NLaysWNLCssmaller{i}"],
                                    axis=1)
            i += 1
        except KeyError:
            break
    print(f"\n{count1} events with ONE photon\n{count2} with TWO photons")
    start = timer()
    xgb_c = xgb.XGBClassifier(**xgb_cps, use_label_encoder = False)
    xgb_r = xgb.XGBRegressor(**xgb_rps)
    xgb_c.fit(trainc_x, trainc_y)
    xgb_r.fit(trainr_x, trainr_y)
    end = timer()
    print("Training (75%) time: {:.9f}s".format(end-start))
    
    start = timer()
    predsclass = xgb_c.predict(testc_x)
    end = timer()
    print("Classifier (25%) time: {:.9f}s".format(end-start))
    
    start = timer()
    pred       = xgb_r.predict(testr_x)
    end = timer()
    print("Regressor (25%) time: {:.9f}s".format(end-start))
    del start, end
    
    conf_matrix= confusion_matrix(testc_y, predsclass, normalize = "true")
    print("\nAccuracy (single class): %.4f" %conf_matrix[0,0])
    print("Accuracy (double class): %.4f" %conf_matrix[1,1])
    pred[predsclass == 0] = 1
    #correctsph = count_correct(pred, testr_y, 10.00)
    corrects20 = count_correct(pred, testr_y)
    corrects10 = count_correct(pred, testr_y, 0.10)
    corrects05 = count_correct(pred, testr_y, 0.05)
    corrects20norm = count_correct_norm(pred, testr_y)
    corrects10norm = count_correct_norm(pred, testr_y, 0.10)
    corrects05norm = count_correct_norm(pred, testr_y, 0.05)
    corrects20mono = count_correct(pred[testr_y == 1], testr_y[testr_y == 1])
    corrects10mono = count_correct(pred[testr_y == 1], testr_y[testr_y == 1], 0.10)
    corrects05mono = count_correct(pred[testr_y == 1], testr_y[testr_y == 1], 0.05)
    corrects20bi = count_correct(pred[testr_y < 1], testr_y[testr_y < 1])
    corrects10bi = count_correct(pred[testr_y < 1], testr_y[testr_y < 1], 0.10)
    corrects05bi = count_correct(pred[testr_y < 1], testr_y[testr_y < 1], 0.05)
    compare = abs(pred - np.array(testr_y)) / np.array(testr_y)
    compare_multi = compare[testr_y != 1]
    compare_mono  = compare[testr_y == 1]
    #print("\n\nAccuracy (Nph) : %.4f" %(correctsph))
    print("\nEvent  Accuracy (0.20): %.4f" %(corrects20))
    print("Photon Accuracy (0.20): %.4f" %(corrects20norm))
    print("One-Ph Accuracy (0.20): %.4f" %(corrects20mono))
    print("Two-Ph Accuracy (0.20): %.4f" %(corrects20bi))
    print("\nEvent  Accuracy (0.10): %.4f" %(corrects10))
    print("Photon Accuracy (0.10): %.4f" %(corrects20norm))
    print("One-Ph Accuracy (0.10): %.4f" %(corrects10mono))
    print("Two-Ph Accuracy (0.10): %.4f" %(corrects10bi))
    print("\nEvent  Accuracy (0.05): %.4f" %(corrects05))
    print("Photon Accuracy (0.05): %.4f" %(corrects20norm))
    print("One-Ph Accuracy (0.05): %.4f" %(corrects05mono))
    print("Two-Ph Accuracy (0.05): %.4f" %(corrects05bi))
    
    
    del corrects20, corrects20norm, corrects20mono, corrects20bi
    del corrects10, corrects10norm, corrects10mono, corrects10bi
    del corrects05, corrects05norm, corrects05mono, corrects05bi
    del conf_matrix



#%%
importances_sum_c = np.sum(importances_c, axis = 0)
importances_sum_r = np.sum(importances_r, axis = 0)
xgb_c.feature_importances = importances_sum_c
xgb_r.feature_importances = importances_sum_r

xgb.plot_importance(xgb_c, importance_type = "weight", max_num_features= 25,
                    title = "Classifier Weight (no lcs)",
                    show_values = False)
pl.tight_layout()
xgb.plot_importance(xgb_r, importance_type = "weight", max_num_features= 25,
                    title = "Regressor Weight (no lcs)", 
                    show_values = False)
pl.tight_layout()
#%%
nf = 25
most_c, columns_most_c = zip(*sorted(zip(importances_sum_c, columns))) 
most_r, columns_most_r = zip(*sorted(zip(importances_sum_r, columns))) 
least_c, columns_least_c = zip(*sorted(zip(importances_sum_c, columns))) 
least_r, columns_least_r = zip(*sorted(zip(importances_sum_r, columns))) 

most_c = most_c[-nf:]
most_r = most_r[-nf:]
columns_most_c = columns_most_c[-nf:]
columns_most_r = columns_most_r[-nf:]
least_c = least_c[:nf]
least_r = least_r[:nf]
columns_least_c = columns_least_c[:nf]
columns_least_r = columns_least_r[:nf]

pl.figure()
pl.title("Classifier Highest Weights")
pl.xlabel("F Score")
y_pos = np.arange(nf)
pl.barh(y_pos, most_c, align = 'center')
pl.yticks(y_pos, labels=columns_most_c)
#pl.gca().invert_yaxis()
pl.grid()
pl.tight_layout()
pl.show()

pl.figure()
pl.title("Regressor Highest Weights")
pl.xlabel("F Score")
y_pos = np.arange(nf)
pl.barh(y_pos, most_r, align = 'center')
pl.yticks(y_pos, labels=columns_most_r)
#pl.gca().invert_yaxis()
pl.grid()
pl.tight_layout()
pl.show()

pl.figure()
pl.title("Classifier Lowest Weights")
pl.xlabel("F Score")
y_pos = np.arange(nf)
pl.barh(y_pos, least_c, align = 'center')
pl.yticks(y_pos, labels=columns_least_c)
pl.gca().invert_yaxis()
pl.grid()
pl.tight_layout()
pl.show()

pl.figure()
pl.title("Regressor Lowest Weights")
pl.xlabel("F Score")
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
    