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

#%% DOWNLOAD & ARRANGE DATA
#______________________________________________________________________________
Elevel = "200max"
Nout = 4
balance = 0
infolayers = 1

data = pd.read_csv("../myData/TSTree_SimFromCP" + Elevel + ".txt", sep = ";")
data = data.apply(pd.to_numeric, errors='coerce')

deltas = np.array(sorted(np.array(data['DeltaCPVtxXY']), reverse = True))
deltas = deltas[deltas==0]
Nsingle = len(deltas)
Ndouble = len(data) - len(deltas)
data_double = data.drop(data.index[Ndouble:])
data_single = data.drop(data.index[:-Nsingle+1])
idx = np.random.permutation(data_double.index)
data_double = data_double.reindex(idx)

if balance == True: #same number of events per class
    diff = Nsingle - Ndouble
    if diff > 0:
        data_single = data_single.drop(data.index[diff:])
    else:
        data_double = data_double.drop(data.index[diff:])

data = pd.concat((data_single, data_double))
idx = np.random.permutation(data.index)
data = data.reindex(idx)


if not infolayers: #do not use Ei, LNi columns
    i = 1
    while i > 0:
        try:
            data = data.drop([f"E{i}", f"LN{i}"], axis = 1)
            i += 1
        except KeyError:
            break

#Remove NaNs
i = 0
while i < len(data):
    row_i = data.loc[i]
    if row_i.isnull().values.any():
        data = data.drop(i)
    i += 1


deltas = np.array(data['DeltaCPVtxXY'])
deltas = deltas[deltas==0]
count0 = len(deltas)
count2 = len(data) - len(deltas)

x = data.iloc[:, 2:-Nout]
y = data.iloc[:, -2]
deltacp = data.iloc[:, 0]
deltalay = data.iloc[:, 1]
train_x, test_x, train_y, test_y = train_test_split(x, y,
                                                    test_size = 0.2, 
                                                    random_state = 123)
deltacp_test = deltacp.loc[test_y.index]
data_xgb = xgb.DMatrix(data = x, label = y)

#%%
# USE XGBOOST - CLASSIFICATION + REGRESSION
#______________________________________________________________________________
pthrs = [0.99, 0.95, 0.90]

xgb_rps = {
          "learning_rate"    : 0.3, #learning rate, small->high num_round
          "max_depth"        : 4,
          "n_estimators"     :25, #boost
          "objective"        : "reg:squarederror",
          "tree_method"      : "hist", #or "gpu-hist", fastest ones
          "subsample"        : 1,
          "sampling_method"  : "uniform", #gradient_based might be better
          "colsample_bytree" : 1,
          "reg_lambda"       : 5, #regularization term (could go up to 10)
          "reg_alpha"        : 5, #regularization term (could go up to 10)
          "n_jobs"           : 1
         }
xgb_cps = {
          "eta"              : 0.3, #learning rate, small->high num_round
          "max_depth"        : 4,
          "n_estimators"     : 25, #boost
          "objective"        : "binary:logistic",
          "tree_method"      : "hist", #or "gpu-hist", fastest ones
          "subsample"        : 1,
          "sampling_method"  : "gradient_based", #or uniform
          "lambda"           : 5, #regularization term (could go up to 10)
          "alpha"            : 5, #regularization term (could go up to 10)
         }

xgb_c = xgb.XGBClassifer(**xgb_cps)
xgb_r = xgb.XGBRegressor(**xgb_rps)


start = timer()
xgb_r.fit(train_x, train_y)
pred_orig = xgb_r.predict(test_x)
end = timer()
print(f"\n{count0} events with ONE photon\n{count2} with TWO photons")
print("Time taken: {:.9f}s".format(end-start)); del start, end
rmse = np.sqrt(MSE(test_y, pred_orig))
print("RMSE          : % f" %(rmse))

for pthr in pthrs:
    pred = pred_orig * 1
    for i in range(len(pred)): 
        if pred_orig[i]>pthr: pred[i]=1
    correctsph = count_correct(pred, test_y, 10.00)
    corrects20 = count_correct(pred, test_y)
    corrects10 = count_correct(pred, test_y, 0.10)
    corrects05 = count_correct(pred, test_y, 0.05)
    corrects20norm = count_correct_norm(pred, test_y)
    corrects10norm = count_correct_norm(pred, test_y, 0.10)
    corrects05norm = count_correct_norm(pred, test_y, 0.05)
    corrects20mono = count_correct(pred[test_y == 1], test_y[test_y == 1])
    corrects10mono = count_correct(pred[test_y == 1], test_y[test_y == 1], 0.10)
    corrects05mono = count_correct(pred[test_y == 1], test_y[test_y == 1], 0.05)
    corrects20bi = count_correct(pred[test_y < 1], test_y[test_y < 1])
    corrects10bi = count_correct(pred[test_y < 1], test_y[test_y < 1], 0.10)
    corrects05bi = count_correct(pred[test_y < 1], test_y[test_y < 1], 0.05)
    compare = abs(pred - np.array(test_y)) / np.array(test_y)
    compare_multi = compare[test_y != 1]
    compare_mono  = compare[test_y == 1]
    print("Accuracy (Nph) : %.4f" %(correctsph))
    print("Event  Accuracy (0.20): %.4f" %(corrects20))
    print("Photon Accuracy (0.20): %.4f" %(corrects20norm))
    print("One-Ph Accuracy (0.20): %.4f" %(corrects20mono))
    print("Two-Ph Accuracy (0.20): %.4f" %(corrects20bi))
    print("Event  Accuracy (0.10): %.4f" %(corrects10))
    print("Photon Accuracy (0.10): %.4f" %(corrects20norm))
    print("One-Ph Accuracy (0.10): %.4f" %(corrects10mono))
    print("Two-Ph Accuracy (0.10): %.4f" %(corrects10bi))
    print("Event  Accuracy (0.05): %.4f" %(corrects05))
    print("Photon Accuracy (0.05): %.4f" %(corrects20norm))
    print("One-Ph Accuracy (0.05): %.4f" %(corrects05mono))
    print("Two-Ph Accuracy (0.05): %.4f" %(corrects05bi))
    print(f"Note: say ONE photon if\npredicted fraction above {pthr}")

    xgb.plot_importance(xgb_r, importance_type = "weight", max_num_features= 25,
                        title = "Weight Importance", show_values = False)
    pl.tight_layout()
    xgb.plot_importance(xgb_r, importance_type = "gain", max_num_features= 25,
                        title = "Gain Importance", show_values = False)
    pl.tight_layout()
    xgb.plot_importance(xgb_r, importance_type = "cover", max_num_features= 25,
                        title = "Cover Importance", show_values = False)
    pl.tight_layout()




#******************************************************************************
#%% DO ALL
#______________________________________________________________________________
Elevel = "200min"
Nout = 4

data = pd.read_csv("../myData/TSTree_SimFromCP" + Elevel + ".txt", sep = ";")
data = data.apply(pd.to_numeric, errors='coerce')
#Remove NaNs
i = 0
while i < len(data):
    row_i = data.loc[i]
    if row_i.isnull().values.any():
        data = data.drop(i)
    i += 1

count0 = 0
for i in data["DeltaCPVtxXY"]:
    if i == 0: count0 += 1
count2 = len(data) - count0
x = data.iloc[:, 2:-Nout]
y = data.iloc[:, -2]
deltacp = data.iloc[:, 0]
deltalay = data.iloc[:, 1]
train_x, test_x, train_y, test_y = train_test_split(x, y,
                                                    test_size = 0.2, 
                                                    random_state = 123)
deltacp_test = deltacp.loc[test_y.index]
data_xgb = xgb.DMatrix(data = x, label = y)

xgb_rps = {
          "learning_rate"    : 0.3, #learning rate, small->high num_round
          "max_depth"        : 4,
          "n_estimators"     :25, #boost
          "objective"        : "reg:squarederror",
          "tree_method"      : "hist", #or "gpu-hist", fastest ones
          "subsample"        : 1,
          "sampling_method"  : "uniform", #gradient_based might be better
          "colsample_bytree" : 1,
          "reg_lambda"       : 5, #regularization term (could go up to 10)
          "reg_alpha"        : 5, #regularization term (could go up to 10)
          "n_jobs"           : 1
         }
xgb_cps = {
          "eta"              : 0.3, #learning rate, small->high num_round
          "max_depth"        : 4,
          "n_estimators"     : 25, #boost
          "objective"        : "binary:logistic",
          "tree_method"      : "hist", #or "gpu-hist", fastest ones
          "subsample"        : 1,
          "sampling_method"  : "gradient_based", #or uniform
          "lambda"           : 5, #regularization term (could go up to 10)
          "alpha"            : 5, #regularization term (could go up to 10)
         }

xgb_c = xgb.XGBClassifer(**xgb_cps)
xgb_r = xgb.XGBRegressor(**xgb_rps)

# start = timer()
# xgb_r.fit(train_x, train_y)
# pred_orig = xgb_r.predict(test_x)
# end = timer()
# print(f"\n{count0} events with ONE photon\n{count2} with TWO photons")
# print("Time taken: {:.9f}s".format(end-start)); del start, end
# rmse = np.sqrt(MSE(test_y, pred_orig))
# print("RMSE          : %f" %(rmse))

# for pthr in pthrs:
#     pred = pred_orig * 1
#     for i in range(len(pred)): 
#         if pred_orig[i]>pthr: pred[i]=1
#     correctsph = count_correct(pred, test_y, 10.00)
#     corrects20 = count_correct(pred, test_y)
#     corrects10 = count_correct(pred, test_y, 0.10)
#     corrects05 = count_correct(pred, test_y, 0.05)
#     corrects20norm = count_correct_norm(pred, test_y)
#     corrects10norm = count_correct_norm(pred, test_y, 0.10)
#     corrects05norm = count_correct_norm(pred, test_y, 0.05)
#     corrects20mono = count_correct(pred[test_y == 1], test_y[test_y == 1])
#     corrects10mono = count_correct(pred[test_y == 1], test_y[test_y == 1], 0.10)
#     corrects05mono = count_correct(pred[test_y == 1], test_y[test_y == 1], 0.05)
#     corrects20bi = count_correct(pred[test_y < 1], test_y[test_y < 1])
#     corrects10bi = count_correct(pred[test_y < 1], test_y[test_y < 1], 0.10)
#     corrects05bi = count_correct(pred[test_y < 1], test_y[test_y < 1], 0.05)
#     compare = abs(pred - np.array(test_y)) / np.array(test_y)
#     compare_multi = compare[test_y != 1]
#     compare_mono  = compare[test_y == 1]
#     print("\n\nAccuracy (Nph) : %.4f" %(correctsph))
#     print("\nEvent  Accuracy (0.20): %.4f" %(corrects20))
#     print("Photon Accuracy (0.20): %.4f" %(corrects20norm))
#     print("One-Ph Accuracy (0.20): %.4f" %(corrects20mono))
#     print("Two-Ph Accuracy (0.20): %.4f" %(corrects20bi))
#     print("\nEvent  Accuracy (0.10): %.4f" %(corrects10))
#     print("Photon Accuracy (0.10): %.4f" %(corrects20norm))
#     print("One-Ph Accuracy (0.10): %.4f" %(corrects10mono))
#     print("Two-Ph Accuracy (0.10): %.4f" %(corrects10bi))
#     print("\nEvent  Accuracy (0.05): %.4f" %(corrects05))
#     print("Photon Accuracy (0.05): %.4f" %(corrects20norm))
#     print("One-Ph Accuracy (0.05): %.4f" %(corrects05mono))
#     print("Two-Ph Accuracy (0.05): %.4f" %(corrects05bi))
#     print(f"Note: say ONE photon if\npredicted fraction above {pthr}")

# xgb.plot_importance(xgb_r, importance_type = "weight", max_num_features= 25,
#                     title = "Weight Importance", show_values = False)
# pl.tight_layout()
# xgb.plot_importance(xgb_r, importance_type = "gain", max_num_features= 25,
#                     title = "Gain Importance", show_values = False)
# pl.tight_layout()


print("\n********************************************************************")
print("NOW BALANCE **********************************************************")
deltas = np.array(data['DeltaCPVtxXY'])
deltas = deltas[deltas==0]
Nsingle = len(deltas)
Ndouble = len(data) - len(deltas)
data_double = data.drop(data.index[Ndouble:])
data_single = data.drop(data.index[:-Nsingle])
idx = np.random.permutation(data_double.index)
data_double = data_double.reindex(idx)
diff = Nsingle - Ndouble

if diff > 0:
    data_single = data_single.reset_index().iloc[:,1:]
    data_single = data_single.drop(data.index[:diff])
else:
    #lendrop = len(data_single)+diff
    data_double = data_double.reset_index().iloc[:,1:]
    data_double = data_double.drop(data.index[:-diff])

data = pd.concat((data_single, data_double))
idx = np.random.permutation(data.index)
data = data.reindex(idx)

count0 = 0
for i in data["DeltaCPVtxXY"]:
    if i == 0: count0 += 1
count2 = len(data) - count0
x = data.iloc[:, 2:-Nout]
y = data.iloc[:, -2]
deltacp = data.iloc[:, 0]
deltalay = data.iloc[:, 1]
train_x, test_x, train_y, test_y = train_test_split(x, y,
                                                    test_size = 0.2, 
                                                    random_state = 123)
deltacp_test = deltacp.loc[test_y.index]
data_xgb = xgb.DMatrix(data = x, label = y)
xgb_r = xgb.XGBRegressor(**xgb_ps)

start = timer()
xgb_r.fit(train_x, train_y)
pred_orig = xgb_r.predict(test_x)
end = timer()
print(f"\n{count0} events with ONE photon\n{count2} with TWO photons")
print("Time taken: {:.9f}s".format(end-start)); del start, end
rmse = np.sqrt(MSE(test_y, pred_orig))
print("RMSE          : % f" %(rmse))

for pthr in pthrs:
    pred = pred_orig * 1
    for i in range(len(pred)): 
        if pred_orig[i]>pthr: pred[i]=1
    correctsph = count_correct(pred, test_y, 10.00)
    corrects20 = count_correct(pred, test_y)
    corrects10 = count_correct(pred, test_y, 0.10)
    corrects05 = count_correct(pred, test_y, 0.05)
    corrects20norm = count_correct_norm(pred, test_y)
    corrects10norm = count_correct_norm(pred, test_y, 0.10)
    corrects05norm = count_correct_norm(pred, test_y, 0.05)
    corrects20mono = count_correct(pred[test_y == 1], test_y[test_y == 1])
    corrects10mono = count_correct(pred[test_y == 1], test_y[test_y == 1], 0.10)
    corrects05mono = count_correct(pred[test_y == 1], test_y[test_y == 1], 0.05)
    corrects20bi = count_correct(pred[test_y < 1], test_y[test_y < 1])
    corrects10bi = count_correct(pred[test_y < 1], test_y[test_y < 1], 0.10)
    corrects05bi = count_correct(pred[test_y < 1], test_y[test_y < 1], 0.05)
    compare = abs(pred - np.array(test_y)) / np.array(test_y)
    compare_multi = compare[test_y != 1]
    compare_mono  = compare[test_y == 1]
    print("\n\nAccuracy (Nph) : %.4f" %(correctsph))
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
    print(f"Note: say ONE photon if\npredicted fraction above {pthr}")

xgb.plot_importance(xgb_r, importance_type = "weight", max_num_features= 25,
                    title = "Weight Importance", show_values = False)
pl.tight_layout()
xgb.plot_importance(xgb_r, importance_type = "gain", max_num_features= 25,
                    title = "Gain Importance", show_values = False)
pl.tight_layout()


print("\n********************************************************************")
print("NOW REMOVE LCs *******************************************************")
i = 1
while i > 0:
    try:
        data = data.drop([f"E{i}", f"LN{i}"], axis = 1)
        i += 1
    except KeyError:
        break
    
count0 = 0
for i in data["DeltaCPVtxXY"]:
    if i == 0: count0 += 1
count2 = len(data) - count0
x = data.iloc[:, 2:-Nout]
y = data.iloc[:, -2]
deltacp = data.iloc[:, 0]
deltalay = data.iloc[:, 1]
train_x, test_x, train_y, test_y = train_test_split(x, y,
                                                    test_size = 0.2, 
                                                    random_state = 123)
deltacp_test = deltacp.loc[test_y.index]
data_xgb = xgb.DMatrix(data = x, label = y)

xgb_r = xgb.XGBRegressor(**xgb_ps)
start = timer()
xgb_r.fit(train_x, train_y)
pred_orig = xgb_r.predict(test_x)
end = timer()
print(f"\n{count0} events with ONE photon\n{count2} with TWO photons")
print("Time taken: {:.9f}s".format(end-start)); del start, end
rmse = np.sqrt(MSE(test_y, pred_orig))
print("RMSE          : % f" %(rmse))

for pthr in pthrs:
    pred = pred_orig * 1
    for i in range(len(pred)): 
        if pred_orig[i]>pthr: pred[i]=1
    correctsph = count_correct(pred, test_y, 10.00)
    corrects20 = count_correct(pred, test_y)
    corrects10 = count_correct(pred, test_y, 0.10)
    corrects05 = count_correct(pred, test_y, 0.05)
    corrects20norm = count_correct_norm(pred, test_y)
    corrects10norm = count_correct_norm(pred, test_y, 0.10)
    corrects05norm = count_correct_norm(pred, test_y, 0.05)
    corrects20mono = count_correct(pred[test_y == 1], test_y[test_y == 1])
    corrects10mono = count_correct(pred[test_y == 1], test_y[test_y == 1], 0.10)
    corrects05mono = count_correct(pred[test_y == 1], test_y[test_y == 1], 0.05)
    corrects20bi = count_correct(pred[test_y < 1], test_y[test_y < 1])
    corrects10bi = count_correct(pred[test_y < 1], test_y[test_y < 1], 0.10)
    corrects05bi = count_correct(pred[test_y < 1], test_y[test_y < 1], 0.05)
    compare = abs(pred - np.array(test_y)) / np.array(test_y)
    compare_multi = compare[test_y != 1]
    compare_mono  = compare[test_y == 1]
    print("\n\nAccuracy (Nph) : %.4f" %(correctsph))
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
    print(f"Note: say ONE photon if\npredicted fraction above {pthr}")

xgb.plot_importance(xgb_r, importance_type = "weight", max_num_features= 25,
                    title = "Weight Importance", show_values = False)
pl.tight_layout()
xgb.plot_importance(xgb_r, importance_type = "gain", max_num_features= 25,
                    title = "Gain Importance", show_values = False)
pl.tight_layout()
# xgb.plot_importance(xgb_r, importance_type = "cover", max_num_features= 25,
#                     title = "Cover Importance", show_values = False)
# pl.tight_layout()


#%%
# USE XGBOOST REGRESSION - CROSS VALIDATION
#______________________________________________________________________________
data_xgb = xgb.DMatrix(data = x, label = y)
xgb_ps = {
          "eta"              : 0.3, #learning rate, small->high num_round
          "max_depth"        : 4,
          "objective"        : "reg:squarederror",
          "tree_method"      : "hist", #or "gpu-hist", fastest ones
          "subsample"        : 1,
          "sampling_method"  : "uniform", #gradient_based might be better
          "lambda"           : 5, #regularization term (could go up to 10)
          "alpha"            : 5, #regularization term (could go up to 10)
          "n_jobs"           : multiprocessing.cpu_count()// 2
          }

cv_results = xgb.cv(dtrain = data_xgb,
                    params = xgb_ps,
                    num_boost_round = 25,
                    nfold = 4,
                    early_stopping_rounds = 10,
                    metrics = "rmse",
                    as_pandas = True)

pl.figure()
rounds = np.array(cv_results.index)
rmses  = np.array(cv_results["test-rmse-mean"])
rmstd  = np.array(cv_results["test-rmse-std"])
pl.errorbar(rounds, rmses, rmstd, fmt = '.', ls = 'dashed', 
           capsize = 4, label = 'Test Imporvement')
pl.legend()
pl.xlabel("Round n.")
pl.ylabel("RMSE")
pl.show()





#%%
# USE XGBOOST - CLASSIFICATION
#______________________________________________________________________________
Nout = 4

for Elevel in ["100max", "200max", "200min"]: 
    
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
                        title = "Weight Importance", show_values = False)
    pl.tight_layout()
    xgb.plot_importance(xgb_c, importance_type = "gain", max_num_features= 25,
                        title = "Gain Importance", show_values = False)
    pl.tight_layout()


















