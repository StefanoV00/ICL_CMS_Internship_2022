# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 12:38:10 2022

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
from keras.datasets import mnist
from torch.optim import SGD, Adam, Adamax
import torch.nn as nn
import torch.nn.functional as nnf
from sklearn.metrics import mean_squared_error
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error as MSE

from PyTorchNetwork import MyNetwork, uniform_seq, uniform_norm_seq
from PyTorchNetwork import count_correct, count_correct1, get_tensors

#%%
(train_X, train_y), (test_X, test_y) = mnist.load_data()
x = np.concatenate((train_X, test_X))
x = x.reshape(len(x), 28*28)
x = x.astype(np.float64)#/1000
y = np.concatenate((train_y, test_y))

Nin = len(x[0])
Nout = 10
yout = []
for i in range(len(y)):
    yout.append( np.zeros(10) )
    yout[-1][y[i]] += 1
yout = np.array(yout)
H = 2
loss_fn = nn.MSELoss()
opt_name = "Adam"
learning_rate = 1e-3
Nh = 32


uni_arch = uniform_seq(Nin, Nout, [Nh, Nh], prob = 0.0)
network = MyNetwork(uni_arch, nn.Softmax(dim=1))
optimizer = locals()[f"{opt_name}"]( network.parameters(), 
                                     lr = learning_rate)
def digits_is_correct(pred_i, y_i):
    pred_ij = np.argmax(pred_i)
    return y_i[pred_ij] == 1
def digits_count_correct(pred, y, thr = 0):
    count = 0
    all_e = 0
    for (p_i, y_i) in zip(pred, y):
        c = int(digits_is_correct(p_i, y_i))
        if not np.isnan(c):
            count += c
        all_e += 1
    return count / all_e

train_x, test_x, train_y, test_y = train_test_split(x, yout,
                                         test_size = 0.2)


for param in network.parameters():
  print(param.data)

results = network.full_train(train_x, train_y, 
                                    loss_fn, 
                                    optimizer, 
                                    batch_size = 1000,
                                    epochs = 20, 
                                    count_correct = digits_count_correct,
                                    correct_thr = 0.)
for param in network.parameters():
  print(param.data)

test_x, test_y = get_tensors(test_x, test_y, dtype = torch.float32)

print(network.test(test_x, test_y, loss_fn, digits_count_correct))






