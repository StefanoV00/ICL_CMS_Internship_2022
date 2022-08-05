# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 12:18:27 2022

@author: Stefano

File for ML classes.
"""

import numpy as np
import pandas as pd


import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error as MSE
from sklearn.multioutput import MultiOutputRegressor as MOR




#******************************************************************************
#* "CORRECTNESS" FUNCTIONS
#______________________________________________________________________________

def is_correct(pred_i, y_i, thr):
    p_score1, p_score2, p_fr1, p_fr2 = pred_i;
    correct_n = ((p_score1-p_score2)*(-1)**(y_i[0]))<0
    correct_E = abs(p_fr1 - y_i[2]) / y_i[2] < 0.2
    if y_i[1]: #there are 2 photons
        correct_E = correct_E and abs(p_fr2 - y_i[3]) / y_i[3] < thr
    return correct_n and correct_E

def count_correct(pred, y, thr):
    count = 0
    all_e = 0
    for (p_i, y_i) in zip(pred, y):
        count += int(is_correct(p_i, y_i, thr))
        all_e += 1
    return count / all_e




#******************************************************************************
#* MULTI-OUTPUT-XGB-REGRESSION
#******************************************************************************
#______________________________________________________________________________
def XGBMultiRegCV(df, Nout, nfold=3, count_correct = count_correct,
                  correct_thr = 0.2, printcv = True, **kwargs):
    """
    Algorithm for CrossValidation of a MultiOutput Regression using xgboost's 
    Regressor wrapped by sklearn.multioutput.MultiOutputRegressor. \n
    Note: the multiple outputs are treated separately, so it is not a perfect
    multioutput regressor. It likely fails to capture relations between outputs.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with inputs and outputs (which are in the LAST Nout columns).
    Nout : int
        Number of outputs. Everything else is input.
    nfold : int, default is 3
        Number of k-folds in the cross validation.
    count_correct : function, default is helperML.count_correct()
        Function taking as argument prediction array, expected arrays, threshold
        and returning fraction of times the outputs result in correct 
        classification.
        If 0, instead, don't count.
    correct_thr : float in [0,1]
        Third argument of count_correct, it is the relative error between 
        results under which the prediction is considered correct.
    printcv : bool
        Print "CV fold i"
    **kwargs : xgboost.XGBRegressor kwargs
        MultiCV defines a XGBRegressor by itself using given kwargs.
        
        
    Returns
    -------
    MOR  : sklearn.multioutput.MultiOutputRegressor
        The trained Regressor
    RMSEs : np.ndarray (nfold, Nout)
        RMSEs of each testing
    Accuracies : np.ndarray (nfold), if count_correct is not 0
        Fraction of correct classification of each testing.
        It requires the count_correct function to be defined.
    """
    k = nfold
    print(f'{nfold} Cross Validation start')
    n = int(np.floor( len(df) / k))
    
    rmses = []
    corrects = []
    for i in range(k):
        if printcv: print(f"CV fold {i+1}")
        if i == k - 1:
            #at last one also include residual entries, thus weirder
            training = df.drop(df.iloc[n*i : n*(i+1) -(n*(i+1)-len(df)), :].index)
            testing  = df.iloc[n*i : n*(i+1) - (n*(i+1)-len(df)), :]
        else:
            training = df.drop(df.iloc[n*i : n*(i+1), :].index)
            testing  = df.iloc[n * i:n * (i + 1), :]

        train_x = training.iloc[:, :-Nout]
        train_y = training.iloc[:, -Nout:]
        test_x  = testing.iloc [:, :-Nout]
        test_y  = testing.iloc [:, -Nout:]
        
        xgb_r = xgb.XGBRegressor(**kwargs)
        multi_reg = MOR(xgb_r,n_jobs = -1).fit(train_x, train_y)   
        
        pred = multi_reg.predict(test_x)
        #rmses.append(MSE(test_y,pred,multioutput='raw_values,squared=False))
        rmses.append(np.sqrt(np.mean(np.power(pred-test_y, 2),axis = 0)))
        if count_correct:
            corrects.append(count_correct(pred, test_y.values, correct_thr))
     
    results = [multi_reg, np.array(rmses)]
    print('Cross Validation end')
    rmse  = np.mean(rmses, axis = None)
    rmstd = np.std(rmses, axis = None)
    print (f'RMSE    = {rmse:.5f} +- {rmstd:.5f}')
    if count_correct:
        correctmean = np.mean(corrects)
        correctstd  = np.std(corrects)
        print (f'Accuracy: {correctmean:.5f} +- {correctstd:.5f}')
        results.append(np.array(corrects))
    return results





#following yields exact same results as above, with xgb_r defined at beginning
def MultiCV(df, Nout, xgb_r , nfold=3, count_correct = count_correct):
    """
    Algorithm for CrossValidation of a MultiOutput Regression using xgboost's 
    Regressor wrapped by sklearn.multioutput.MultiOutputRegressor. \n
    Note: the multiple outputs are treated separately, so it is not a perfect
    multioutput regressor. It likely fails to capture relations between outputs.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with inputs and outputs (which are in the LAST Nout columns).
    Nout : int
        Number of outputs. Everything else is input.
    xgb_r : xgboost.XGBRegressor
        The configured regressor
    nfold : int, default is 3
        Number of k-folds in the cross validation.
    count_correct : function, default is helperML.count_correct()
        Function taking as argument prediction and expected arrays and
        returning fraction of times the Nout results in correct
        classification.
        If 0, instead, do not count. 

    Returns
    -------
    MOR  : sklearn.multioutput.MultiOutputRegressor
        The trained Regressor
    RMSEs : np.ndarray (nfold, Nout)
        RMSEs of each testing
    Accuracies : np.ndarray (nfold)
        Fraction of correct classification of each testing.
        It requires the count_correct function to be defined.
    """
    k = nfold
    print('training start')
    n = int(np.floor( len(df) / k))
    
    rmses = []
    corrects = []
    for i in range(k):
        if i == k - 1:
            #at last one also include residual entries, thus weirder
            training=df.drop(df.iloc[n*i : n*(i+1) -(n*(i+1)-len(df)), :].index)
            testing  = df.iloc[n*i : n*(i+1) - (n*(i+1)-len(df)), :]
        else:
            training = df.drop(df.iloc[n*i : n*(i+1), :].index)
            testing  = df.iloc[n * i:n * (i + 1), :]

        train_x = training.iloc[:, :-Nout]
        train_y = training.iloc[:, -Nout:]
        test_x  = testing .iloc[:, :-Nout]
        test_y  = testing .iloc[:, -Nout:]
        multi_reg = MOR(xgb_r,n_jobs = -1).fit(train_x, train_y)   
        
        pred = multi_reg.predict(test_x)
        rmses.append(np.sqrt(np.mean(np.power(pred-test_y, 2), axis = 0)))
        if count_correct:
            corrects.append(count_correct(pred, test_y))
     
    results = [MOR, np.array(rmses)]
    print('Training end')
    rmse  = np.mean(rmses, axis = None)
    rmstd = np.std(rmses, axis = None)
    print ('RMSE = ', rmse, "+-", rmstd)
    if count_correct:
        correctmean = np.mean(corrects)
        correctstd  = np.std(corrects)
        print ('Accuracy : ', correctmean, "+-", correctstd)
        results.append(np.array(corrects))
    return results























