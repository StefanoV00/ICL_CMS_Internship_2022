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

import torch
from torch.optim import SGD, Adam
import torch.nn as nn
import torch.nn.functional as nnf
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error as MSE

from PyTorchNetwork import MyNetwork, uniform_seq, uniform_norm_seq
from PyTorchNetwork import count_correct, count_correct1, get_tensors

#%% DOWNLOAD & ARRANGE DATA
#______________________________________________________________________________
Elevel = "200min_norm"
Nout = 4
balance = 0
infolayers = 1

data = pd.read_csv("../myData/TSTree_SimFromCP" + Elevel + ".txt", sep = ";")
data = data.apply(pd.to_numeric, errors='coerce')

#Remove NaNs
i = 0
while i < len(data):
    row_i = data.loc[i]
    if row_i.isnull().values.any():
        data = data.drop(i)
    i += 1


# Balance number of events per class, if required
if balance == True: 
    data = data.sort_values("DeltaCPVtxXY", ascending = False)
    count0 = 0
    for i in data["DeltaCPVtxXY"]:
        if i == 0: count0 += 1
    diff = 2*count0 - len(data)
    if diff > 0:
        data = data.drop(data.index[-diff:])
    else:
        data = data.drop(data.index[:-diff])
    data = data.reset_index().iloc[:,1:]

#Remove Ei and LNi columns, if required
if not infolayers: 
    i = 1
    while i > 0:
        try:
            data = data.drop([f"E{i}", f"LN{i}"], axis = 1)
        except KeyError:
            break

count0 = 0
for i in data["DeltaCPVtxXY"]:
    if i == 0: count0 += 1
count2 = len(data) - count0

x = data.iloc[:, 2:-Nout]
y = data.iloc[:, -Nout:]
deltacp = data.iloc[:, 0]
deltalay = data.iloc[:, 1]
train_x, test_x, train_y, test_y = train_test_split(x, y,
                                                    test_size = 0.2, 
                                                    random_state = 123)
deltacp_test = deltacp.loc[test_y.index]

x1 = data.iloc[:, 2:-Nout]
y1 = data.iloc[:, -2]
train_x1, test_x1, train_y1, test_y1 = train_test_split(x1, y1,
                                                    test_size = 0.2, 
                                                    random_state = 123)

#%% NETWORK ANALYSIS
#______________________________________________________________________________
device = "cuda" if torch.cuda.is_available() else "cpu"
Nin = len(x.iloc[0])
H = 2
loss_fn = nn.MSELoss()
opt_name = "Adam"
learning_rate = 1e-3

Nh = int(min( (0.66*Nin+Nout), Nin) // H)
while Nh * Nin + Nh ** H + Nh * Nout > 0.5 * len(train_x) and Nh>4:
    Nh -= 2
corrects_adam = [[], []]   

for Nh in  [16]:
    start = timer()
    print(f"NN has 4 OUTPUTS &\n-{H} hidden layers\n-with {Nh} neurons each")
    print(f"-{opt_name} Optimizer")
    uni_arch = uniform_seq(Nin, Nout, [Nh]*H, prob = 0.0)
    network = MyNetwork(uni_arch)
    optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                         lr = learning_rate)
    test_loss,train_loss, test_correct,train_correct =network.train_test(x, y, 
                                                            loss_fn, 
                                                            optimizer, 
                                                            test_frac = 0.2,
                                                            batch_size = 100,
                                                            epochs = 5, 
                                                            correct_thr = 0.1)
    corrects_adam[0].append(test_correct)
    end = timer()
    print("Time taken: {:.9f}s".format(end-start)); del start, end

    start = timer()
    print(f"NN has 1 OUTPUT &\n-{H} hidden layers\n-with {Nh} neurons each")
    print(f"-{opt_name} Optimizer")
    uni_arch = uniform_seq(Nin, 1, [Nh]*H, prob = 0.0)
    network = MyNetwork(uni_arch)
    optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                         lr = learning_rate)
    test_loss,train_loss,test_correct,train_correct =network.train_test(x1, y1, 
                                                loss_fn, 
                                                optimizer, 
                                                test_frac = 0.2,
                                                batch_size = 100,
                                                epochs = 5,
                                                count_correct=count_correct1,
                                                correct_thr = 0.1)
    corrects_adam[1].append(test_correct)
    end = timer()
    print("Time taken: {:.9f}s".format(end-start)); del start, end

if balance:
    print("\n-----------------------------------------------------------------")
    print("----------------- NOW NO LCs INFO -------------------------------")
    i = 1
    while i > 0:
        try:
            data = data.drop([f"E{i}", f"LN{i}"], axis = 1)
        except KeyError:
            break
    for Nh in  [4, 8, 16, 32]:
        start = timer()
        print(f"NN has 4 OUTPUTS &\n-{H} hidden layers\n-with {Nh} neurons each")
        print(f"-{opt_name} Optimizer")
        uni_arch = uniform_seq(Nin, Nout, [Nh]*H, prob = 0.0)
        network = MyNetwork(uni_arch)
        optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                             lr = learning_rate)
        test_loss,train_loss, test_correct,train_correct =network.train_test(x, y, 
                                                            loss_fn, 
                                                            optimizer, 
                                                            test_frac = 0.2,
                                                            batch_size = 100,
                                                            epochs = 5, 
                                                            correct_thr = 0.1)
        corrects_adam[0].append(test_correct)
        end = timer()
        print("Time taken: {:.9f}s".format(end-start)); del start, end
    
        start = timer()
        print(f"NN has 1 OUTPUT &\n-{H} hidden layers\n-with {Nh} neurons each")
        print(f"-{opt_name} Optimizer")
        uni_arch = uniform_seq(Nin, 1, [Nh]*H, prob = 0.0)
        network = MyNetwork(uni_arch)
        optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                             lr = learning_rate)
        test_loss,train_loss,test_correct,train_correct =network.train_test(x1, y1, 
                                                loss_fn, 
                                                optimizer, 
                                                test_frac = 0.2,
                                                batch_size = 100,
                                                epochs = 5,
                                                count_correct=count_correct1,
                                                correct_thr = 0.1)
        corrects_adam[1].append(test_correct)
        end = timer()
        print("Time taken: {:.9f}s".format(end-start)); del start, end

#%% FOR DEBUGGING
# device = "cuda" if torch.cuda.is_available() else "cpu"
# Nin = len(x.iloc[0])
# H = 2
# loss_fn = nn.MSELoss()
# opt_name = "Adam"
# learning_rate = 0.1

# Nh = int(min( (0.66*Nin+Nout), Nin) // H)
# while Nh * Nin + Nh ** H + Nh * Nout > 0.5 * len(train_x) and Nh>4:
#     Nh -= 2
# corrects_adam = [[], []]   

# # DEBUGGING
# Nh = 8
# uni_arch = uniform_seq(Nin, Nout, [Nh]*H, 
#                        actfunc = nn.ReLU, prob = 0.0)
# network = MyNetwork(uni_arch, resultfunc = torch.sigmoid)
# optimizer = locals()[f"{opt_name}"]( network.parameters(), 
#                                      lr = learning_rate)
# losses, corrects =network.train(train_x, train_y, 
#                                             loss_fn, 
#                                             optimizer, 
#                                             batch_size = 100)










# train_x, train_y = get_tensors(train_x, train_y1, dtype = torch.float32)

# train_data   = torch.utils.data.TensorDataset(train_x, train_y)   
# train_loader = torch.utils.data.DataLoader(train_data, 
#                                            batch_size=100)
    
# losses = []; correct = []
# for batch_n, (X, Y) in enumerate(train_loader):
 
#     # Compute prediction and loss
#     pred = network(X)
#     predsize = pred.size()[1]
#     Ysize = len(Y.size())
#     pred = pred[:,0]
#     loss = loss_fn(pred, Y)
#     lossn = loss.item()
#     break
#%% DO ALL
#______________________________________________________________________________
Elevels = ["100max", "100max_norm", "200max", "200max_norm", "200min", 
           "200min_norm"]
def doallNN(Elevel = "200min_norm", Nout = 4):
    print("******************************************************************")
    print("* ", Elevel, "                                                  *")
    print("******************************************************************")
    from torch.optim import Adam
    from torch.optim import SGD
    opt_name = "Adam"
    H = 3
    if opt_name[:4] == "Adam":
        learning_rate = 1e-3
    elif opt_name == "SGD":
        learning_rate = 1

    
    data = pd.read_csv("../myData/TSTree_SimFromCP"+Elevel+".txt", sep = ";")
    data = data.apply(pd.to_numeric, errors='coerce')
    #Remove NaNs
    i = 0
    while i < len(data):
        row_i = data.loc[i]
        if row_i.isnull().values.any():
            data = data.drop(i)
        i += 1
    data = data.reset_index().iloc[:,1:]
    
    count0 = 0
    for i in data["DeltaCPVtxXY"]:
        if i == 0: count0 += 1
    count2 = len(data) - count0
    
    x = data.iloc[:, 2:-Nout]
    y = data.iloc[:, -Nout:]
    deltacp = data.iloc[:, 0]
    deltalay = data.iloc[:, 1]
    train_x, test_x, train_y, test_y = train_test_split(x, y,
                                                        test_size = 0.2)
    deltacp_test = deltacp.loc[test_y.index]
    
    x1 = data.iloc[:, 2:-Nout]
    y1 = data.iloc[:, -2]
    train_x1, test_x1, train_y1, test_y1 = train_test_split(x1, y1,
                                                        test_size = 0.2)
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    Nin = len(x.iloc[0])
    loss_fn = nn.MSELoss()
    
    for Nh in  [8, 16, 32]:
        start = timer()
        print(f"\nNN has 4 OUTPUTS &\n-{H} hidden layers\n-with {Nh} neurons each")
        print(f"-{opt_name} Optimizer")
        uni_arch = uniform_seq(Nin, Nout, [Nh]*H, prob = 0.0)
        network = MyNetwork(uni_arch)
        optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                             lr = learning_rate)
        network.full_train(train_x, train_y,loss_fn, optimizer, 
                           batch_size = 100, epochs = 20, correct_thr = 0.1)
        network.test(test_x, test_y, loss_fn,
                     correct_thr = 10.0, printcorrect = True)
        network.test(test_x, test_y, loss_fn,
                     correct_thr = 0.20, printcorrect = True)
        network.test(test_x, test_y, loss_fn,
                     correct_thr = 0.10, printcorrect = True)
        network.test(test_x, test_y, loss_fn,
                     correct_thr = 0.05, printcorrect = True)
        end = timer()
        print("Time taken: {:.9f}s".format(end-start)); del start, end
    
    
    print("\n----------------------------------------------------------------")
    print("------------NOW IT IS BALANCED ---------------------------------")
    print("---------------------------------------------------------------\n")
    data = data.sort_values("DeltaCPVtxXY", ascending = False)
    count0 = 0
    for i in data["DeltaCPVtxXY"]:
        if i == 0: count0 += 1
    diff = 2*count0 - len(data)
    if diff > 0:
        data = data.drop(data.index[-diff:])
    else:
        data = data.drop(data.index[:-diff])
    data = data.reset_index().iloc[:,1:]
    x = data.iloc[:, 2:-Nout]
    y = data.iloc[:, -Nout:]
    deltacp = data.iloc[:, 0]
    deltalay = data.iloc[:, 1]
    train_x, test_x, train_y, test_y = train_test_split(x, y,
                                                        test_size = 0.2)
    
    
    for Nh in  [8, 16, 32]:
        start = timer()
        print(f"\nNN has 4 OUTPUTS &\n-{H} hidden layers\n-with {Nh} neurons each")
        print(f"-{opt_name} Optimizer")
        uni_arch = uniform_seq(Nin, Nout, [Nh]*H, prob = 0.0)
        network = MyNetwork(uni_arch)
        optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                             lr = learning_rate)
        network.full_train(train_x, train_y, loss_fn, optimizer, 
                           batch_size = 100, epochs = 20, correct_thr = 0.1)
        network.test(test_x, test_y, loss_fn,
                     correct_thr = 10.0, printcorrect = True)
        network.test(test_x, test_y, loss_fn, 
                     correct_thr = 0.20, printcorrect = True)
        network.test(test_x, test_y, loss_fn, 
                     correct_thr = 0.10, printcorrect = True)
        network.test(test_x, test_y, loss_fn,
                     correct_thr = 0.05, printcorrect = True)
        end = timer()
        print("Time taken: {:.9f}s".format(end-start)); del start, end

    print("\n--------------------------------------------------------------")
    print("----------------- NOW NO LCs INFO ----------------------------")
    i = 1
    while i > 0:
        try:
            data = data.drop([f"E{i}", f"LN{i}"], axis = 1)
        except KeyError:
            break
    x = data.iloc[:, :-Nout]
    #y = data.iloc[:, -Nout:]
    train_x, test_x, train_y, test_y = train_test_split(x, y,
                                                        test_size = 0.2)
    for Nh in  [8, 16, 32]:
        start = timer()
        print(f"\nNN has 4 OUTPUTS &\n-{H} hidden layers\n-with {Nh} neurons each")
        print(f"-{opt_name} Optimizer")
        uni_arch = uniform_seq(Nin, Nout, [Nh]*H, prob = 0.0)
        network = MyNetwork(uni_arch)
        optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                             lr = learning_rate)
        network.full_train(train_x, train_y, loss_fn, optimizer, 
                           batch_size = 100, epochs = 20, correct_thr = 0.1)
        network.test(test_x, test_y, loss_fn, 
                     correct_thr = 10.0, printcorrect = True)
        network.test(test_x, test_y, loss_fn,  
                     correct_thr = 0.20, printcorrect = True)
        network.test(test_x, test_y, loss_fn, 
                     correct_thr = 0.10, printcorrect = True)
        network.test(test_x, test_y, loss_fn, 
                     correct_thr = 0.05, printcorrect = True)
        end = timer()
        print("Time taken: {:.9f}s".format(end-start)); del start, end

for elev in Elevels:
    #if elev[-1] == "m":
    doallNN(elev)






#%% DO ALL
#______________________________________________________________________________
Elevels = ["100max", "100max_norm", "200max", "200max_norm", "200min", 
           "200min_norm"]
def doallNN(Elevel = "200min_norm", Nout = 4):
    print("******************************************************************")
    print("* ", Elevel, "                                                   *")
    print("******************************************************************")
    from torch.optim import Adam
    balance = 0
    infolayers = 1
    
    data = pd.read_csv("../myData/TSTree_SimFromCP"+Elevel+".txt", sep = ";")
    data = data.apply(pd.to_numeric, errors='coerce')
    #Remove NaNs
    i = 0
    while i < len(data):
        row_i = data.loc[i]
        if row_i.isnull().values.any():
            data = data.drop(i)
        i += 1
    data = data.reset_index().iloc[:,1:]
    
    # Balance number of events per class, if required
    if balance == True: 
        data = data.sort_values("DeltaCPVtxXY", ascending = False)
        count0 = 0
        for i in data["DeltaCPVtxXY"]:
            if i == 0: count0 += 1
        diff = 2*count0 - len(data)
        if diff > 0:
            data = data.drop(data.index[-diff:])
        else:
            data = data.drop(data.index[:-diff])
        data = data.reset_index().iloc[:,1:]
    
    #Remove Ei and LNi columns, if required
    if not infolayers: 
        i = 1
        while i > 0:
            try:
                data = data.drop([f"E{i}", f"LN{i}"], axis = 1)
            except KeyError:
                break
    
    count0 = 0
    for i in data["DeltaCPVtxXY"]:
        if i == 0: count0 += 1
    count2 = len(data) - count0
    
    x = data.iloc[:, 2:-Nout]
    y = data.iloc[:, -Nout:]
    deltacp = data.iloc[:, 0]
    deltalay = data.iloc[:, 1]
    train_x, test_x, train_y, test_y = train_test_split(x, y,
                                                        test_size = 0.2, 
                                                        random_state = 123)
    deltacp_test = deltacp.loc[test_y.index]
    
    x1 = data.iloc[:, 2:-Nout]
    y1 = data.iloc[:, -2]
    train_x1, test_x1, train_y1, test_y1 = train_test_split(x1, y1,
                                                        test_size = 0.2, 
                                                        random_state = 123)
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    Nin = len(x.iloc[0])
    H = 2
    loss_fn = nn.MSELoss()
    opt_name = "Adam"
    learning_rate = 0.1 
    
    for Nh in  [4, 8, 16, 32]:
        start = timer()
        print(f"NN has 4 OUTPUTS &\n-{H} hidden layers\n-with {Nh} neurons each")
        print(f"-{opt_name} Optimizer")
        uni_arch = uniform_seq(Nin, Nout, [Nh]*H, prob = 0.0)
        network = MyNetwork(uni_arch)
        optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                             lr = learning_rate)
        network.full_train(train_x, train_y,loss_fn, optimizer, 
                           batch_size = 100, epochs = 3, correct_thr = 0.1)
        network.test(test_x, test_y, loss_fn,
                     correct_thr = 10.0, printcorrect = True)
        network.test(test_x, test_y, loss_fn,
                     correct_thr = 0.20, printcorrect = True)
        network.test(test_x, test_y, loss_fn,
                     correct_thr = 0.10, printcorrect = True)
        network.test(test_x, test_y, loss_fn,
                     correct_thr = 0.05, printcorrect = True)
        end = timer()
        print("Time taken: {:.9f}s".format(end-start)); del start, end
    
        start = timer()
        print(f"\nNN has 1 OUTPUT &\n-{H} hidden layers\n-with {Nh} neurons each")
        print(f"-{opt_name} Optimizer")
        uni_arch = uniform_seq(Nin, 1, [Nh]*H, prob = 0.0)
        network = MyNetwork(uni_arch)
        optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                             lr = learning_rate)
        network.full_train(train_x, train_y1,loss_fn, optimizer, 
                           batch_size = 100, epochs = 3, 
                           count_correct = count_correct1, correct_thr = 0.1)
        network.test(test_x, test_y1, loss_fn, count_correct1,
                     correct_thr = 10.0, printcorrect = True)
        network.test(test_x, test_y1, loss_fn, count_correct1, 
                     correct_thr = 0.20, printcorrect = True)
        network.test(test_x, test_y1, loss_fn, count_correct1,
                     correct_thr = 0.10, printcorrect = True)
        network.test(test_x, test_y1, loss_fn, count_correct1,
                     correct_thr = 0.05, printcorrect = True)
        end = timer()
        print("Time taken: {:.9f}s".format(end-start)); del start, end
    
    
    print("\n----------------------------------------------------------------")
    print("------------NOW IT IS BALANCED ---------------------------------")
    print("---------------------------------------------------------------\n")
    data = data.sort_values("DeltaCPVtxXY", ascending = False)
    count0 = 0
    for i in data["DeltaCPVtxXY"]:
        if i == 0: count0 += 1
    diff = 2*count0 - len(data)
    if diff > 0:
        data = data.drop(data.index[-diff:])
    else:
        data = data.drop(data.index[:-diff])
    data = data.reset_index().iloc[:,1:]
    
    for Nh in  [4, 8, 16, 32]:
        start = timer()
        print(f"NN has 4 OUTPUTS &\n-{H} hidden layers\n-with {Nh} neurons each")
        print(f"-{opt_name} Optimizer")
        uni_arch = uniform_seq(Nin, Nout, [Nh]*H, prob = 0.0)
        network = MyNetwork(uni_arch)
        optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                             lr = learning_rate)
        network.full_train(train_x, train_y,loss_fn, optimizer, 
                           batch_size = 100, epochs = 3, correct_thr = 0.1)
        network.test(test_x, test_y, loss_fn,
                     correct_thr = 10.0, printcorrect = True)
        network.test(test_x, test_y, loss_fn, 
                     correct_thr = 0.20, printcorrect = True)
        network.test(test_x, test_y, loss_fn, 
                     correct_thr = 0.10, printcorrect = True)
        network.test(test_x, test_y, loss_fn,
                     correct_thr = 0.05, printcorrect = True)
        end = timer()
        print("Time taken: {:.9f}s".format(end-start)); del start, end
    
        start = timer()
        print(f"\nNN has 1 OUTPUT &\n-{H} hidden layers\n-with {Nh} neurons each")
        print(f"-{opt_name} Optimizer")
        uni_arch = uniform_seq(Nin, 1, [Nh]*H, prob = 0.0)
        network = MyNetwork(uni_arch)
        optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                             lr = learning_rate)
        network.full_train(train_x, train_y1,loss_fn, optimizer, 
                           batch_size = 100, epochs = 3, 
                           count_correct = count_correct1, correct_thr = 0.1)
        network.test(test_x, test_y1, loss_fn, count_correct1,
                     correct_thr = 10.0, printcorrect = True)
        network.test(test_x, test_y1, loss_fn, count_correct1,
                     correct_thr = 0.20, printcorrect = True)
        network.test(test_x, test_y1, loss_fn, count_correct1,
                     correct_thr = 0.10, printcorrect = True)
        network.test(test_x, test_y1, loss_fn, count_correct1, 
                     correct_thr = 0.05, printcorrect = True)
        end = timer()
        print("Time taken: {:.9f}s".format(end-start)); del start, end


    print("\n--------------------------------------------------------------")
    print("----------------- NOW NO LCs INFO ----------------------------")
    i = 1
    while i > 0:
        try:
            data = data.drop([f"E{i}", f"LN{i}"], axis = 1)
        except KeyError:
            break
    for Nh in  [4, 8, 16, 32]:
        start = timer()
        print(f"NN has 4 OUTPUTS &\n-{H} hidden layers\n-with {Nh} neurons each")
        print(f"-{opt_name} Optimizer")
        uni_arch = uniform_seq(Nin, Nout, [Nh]*H, prob = 0.0)
        network = MyNetwork(uni_arch)
        optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                             lr = learning_rate)
        network.full_train(train_x, train_y,loss_fn, optimizer, 
                           batch_size = 100, epochs = 3, correct_thr = 0.1)
        network.test(test_x, test_y, loss_fn, 
                     correct_thr = 10.0, printcorrect = True)
        network.test(test_x, test_y, loss_fn,  
                     correct_thr = 0.20, printcorrect = True)
        network.test(test_x, test_y, loss_fn, 
                     correct_thr = 0.10, printcorrect = True)
        network.test(test_x, test_y, loss_fn, 
                     correct_thr = 0.05, printcorrect = True)
        end = timer()
        print("Time taken: {:.9f}s".format(end-start)); del start, end
    
        start = timer()
        print(f"\nNN has 1 OUTPUT &\n-{H} hidden layers\n-with {Nh} neurons each")
        print(f"-{opt_name} Optimizer")
        uni_arch = uniform_seq(Nin, 1, [Nh]*H, prob = 0.0)
        network = MyNetwork(uni_arch)
        optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                             lr = learning_rate)
        network.full_train(train_x, train_y1,loss_fn, optimizer, 
                           batch_size = 100, epochs = 3, 
                           count_correct = count_correct1, correct_thr = 0.1)
        network.test(test_x, test_y1, loss_fn, count_correct1, 
                     correct_thr = 10.0, printcorrect = True)
        network.test(test_x, test_y1, loss_fn, count_correct1,
                     correct_thr = 0.20, printcorrect = True)
        network.test(test_x, test_y1, loss_fn, count_correct1,
                     correct_thr = 0.10, printcorrect = True)
        network.test(test_x, test_y1, loss_fn, count_correct1,
                     correct_thr = 0.05, printcorrect = True)
        end = timer()
        print("Time taken: {:.9f}s".format(end-start)); del start, end

for elev in Elevels:
    if elev[-1] == "m":
        doallNN(elev)
