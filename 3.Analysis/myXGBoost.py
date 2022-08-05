# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 15:52:44 2022

@author: Stefano
"""

import numpy as np
import pandas as pd
import random as rand

import multiprocessing
from uncertainties import ufloat


import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error as MSE
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.metrics import plot_confusion_matrix

xgb_ps = {
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
# could also pass more parameters, as shown in 
# https://xgboost.readthedocs.io/en/stable/python/python_api.html
xgb_r = xgb.XGBRegressor(**xgb_ps)


def myCRloss():
    grad = 0
    hessian = 0
    return grad, hessian
