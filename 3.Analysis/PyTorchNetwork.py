# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 20:21:46 2022

@author: Stefano
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from timeit import default_timer as timer

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
import torch.nn as nn
import torch.nn.functional as nnf

from sklearn.model_selection import train_test_split




#******************************************************************************
#* "HELPING" FUNCTIONS
#______________________________________________________________________________

def is_correct(pred_i, y_i, thr):
    p_score1, p_score2, p_fr1, p_fr2 = pred_i;
    correct_n = ((p_score1-p_score2)*(-1)**(y_i[0]))<0
    correct_E = abs(p_fr1 - y_i[2]) / y_i[2] < thr
    if y_i[1]: #there are 2 photons
        correct_E = correct_E and abs(p_fr2 - y_i[3]) / y_i[3] < thr
    return correct_n and correct_E
def count_correct(pred, y, thr):
    count = 0
    all_e = 0
    for (p_i, y_i) in zip(pred, y):
        c = int(is_correct(p_i, y_i, thr))
        if not np.isnan(c):
            count += c
        all_e += 1
    return count / all_e


def is_correct1(pred_i, y_i, thr):
    if y_i == 1 or pred_i == 1:
        correct_E = pred_i - y_i == 0
    else:
        correct_E = abs(pred_i - y_i) / (1 - y_i) < thr
    return correct_E
def count_correct1(pred, y, thr):
    count = 0
    all_e = 0
    for (p_i, y_i) in zip(pred, y):
        count += int(is_correct1(p_i, y_i, thr))
        all_e += 1
    return count / all_e


def get_tensors(x, y, dtype = torch.float32):
    if isinstance(x, pd.core.frame.DataFrame) \
    or isinstance(x, pd.core.series.Series):
        x = torch.tensor(x.values, dtype=torch.float32)
    elif not isinstance(x, torch.Tensor):
        x = torch.tensor(x, dtype=torch.float32)
    if isinstance(y, pd.core.frame.DataFrame) \
    or isinstance(y, pd.core.series.Series):
        y = torch.tensor(y.values, dtype=torch.float32)
    elif not isinstance(y, torch.Tensor):
        y = torch.tensor(y, dtype=torch.float32)
    return x, y




#******************************************************************************
#* PYTORCH NEURAL NETWORK
#******************************************************************************
#______________________________________________________________________________
class MyNetwork(nn.Module):
    """
    Generate a model of a Neural Network inheriting from torch.nn.Module. 

    Parameters
    ----------
    sec_architecture : torch.nn.modules.container.Sequential
        The full architecture of the network, expressed as a sequential 
        container of layers and activations. \n
        A simple example could be:\n
            
            nn.Sequential(nn.Linear(n_in, n_hidden_1),\n
                          nn.sigmoid(),\n
                          nn.Linear(n_hidden_1, n_hidden_2),\n
                          nn.sigmoid(),\n
                          nn.Linear(n_hidden_2, n_out))

    returnfunc : callable or 0, default is torch.sigmoid
        Funcion operated by forward method before on outputs before
        returning them. Another alternative could be nn.Softmax(dim=1).\n
        If 0, don't apply anything.
    
    clip : float, default is 0, i.e. no clipping
        clipping value of gradients

    """
    def __init__(self, seq_architecture, resultfunc = torch.sigmoid, clip = 0):

        super().__init__()
        self.flatten = nn.Flatten()
        self.sequence = seq_architecture
        self.resultfunc = resultfunc
        self.clip = clip
    
    
    def forward(self, X):
        #X = self.flatten(X)
        out = self.sequence(X)
        if self.resultfunc:
            out = self.resultfunc(out) #X = nnf.softmax(X, dim=1)
        return out
    
    
    def backward(self, loss, optimizer): 
        optimizer.zero_grad()
        loss.backward()
        if self.clip:
            torch.nn.utils.clip_grad_norm_(self.parameters(), self.clip)
        optimizer.step()
 
    
    def train (self, train_x, train_y, loss_fn, optimizer, 
               batch_size = 100, count_correct = count_correct, 
               correct_thr = 0.2):
        """
        One epoch of training.

        Parameters
        ----------
        train_x : list, ndarray, DataFrame or tensor
        train_y : list, ndarray, DataFrame or tensor
        loss_fn : PyTorch loss function object
            It is for instance
            - loss_fn = nn.MSELoss()          #NOTE PARENTHESES
            - loss_fn = nn.CrossEntropyLoss() #NOTE PARENTHESES
        optimizer : PyTorch optimizer object
        batch_size : int, default is 100.
        count_correct : function, default is helperML.count_correct()
            Function taking as argument prediction array, expected arrays, threshold
            and returning fraction of times the outputs result in correct 
            classification.
            If 0, instead, don't count.
        correct_thr : float in [0,1]
            Third argument of count_correct, it is the relative error between 
            results under which the prediction is considered correct.

        Returns
        -------
            training losses : np.ndarray (number of batches)
                loss value per batch
            training accuracy : np.ndarray (number of batches)
                correct fraction per batch, returned if count_correct != 0 
        """
        
        train_x, train_y = get_tensors(train_x, train_y, dtype = torch.float32)
        
        train_data   = torch.utils.data.TensorDataset(train_x, train_y)   
        train_loader = torch.utils.data.DataLoader(train_data, 
                                                   batch_size=batch_size)
               
        losses = []; correct = []
        for batch_n, (X, y) in enumerate(train_loader):
    	    
            # Compute prediction and loss
            pred = self(X)
            if pred.size()[1] and len(y.size()) == 1:
                pred = pred[:,0]
            loss = loss_fn(pred, y)
            
            # Record
            losses.append(loss.item())
            if count_correct:
                correct.append(count_correct(pred.detach().numpy(), 
                                             y.detach().numpy(),
                                             correct_thr) )
            #Update
            self.backward(loss, optimizer)
            
        results = np.array(losses)
        if count_correct:
            results = [results]
            results.append( np.array(correct))
        return results
    
    
    def full_train (self, train_x, train_y, loss_fn, optimizer, 
                    batch_size = 100, epochs = 100,
                    count_correct = count_correct, correct_thr = 0.2):
        """
        Training over <epochs> epochs.

        Parameters
        ----------
        train_x : list, ndarray, DataFrame or tensor
        train_y : list, ndarray, DataFrame or tensor
        loss_fn : PyTorch loss function object
            It is for instance
            - loss_fn = nn.MSELoss()          #NOTE PARENTHESES
            - loss_fn = nn.CrossEntropyLoss() #NOTE PARENTHESES
        optimizer : PyTorch optimizer object
        batch_size : int, default 100
        epochs : int, default 100
        count_correct : function, default is helperML.count_correct()
            Function taking as argument prediction array, expected arrays, threshold
            and returning fraction of times the outputs result in correct 
            classification.
            If 0, instead, don't count.
        correct_thr : float in [0,1]
            Third argument of count_correct, it is the relative error between 
            results under which the prediction is considered correct.
        
        Returns
        -------
            training losses : np.ndarray (epochs)
                loss value per epoch
            training accuracy : np.ndarray (epochs)
                correct fraction per epoch, returned if count_correct != 0 
        """
        losses = []
        correct = []
        train_x, train_y = get_tensors(train_x, train_y, dtype = torch.float32)
        for i in range(epochs):
            #shuffle training datasets for new different batches
            idx = torch.tensor(np.random.permutation(len(train_x)))
            train_x = train_x.index_select(0, idx)
            train_y = train_y.index_select(0, idx)
            # idx = np.random.permutation(train_x.index)
            # train_x.reindex(idx)
            # train_y.reindex(idx)
            # go on
            e_results = self.train(train_x, train_y, loss_fn,
                                    optimizer, batch_size,
                                    count_correct = count_correct, 
                                    correct_thr = correct_thr)
            if count_correct:
                e_losses  = e_results[0]
                e_correct = e_results[1]
                losses.append(np.mean(e_losses))
                correct.append(np.mean(e_correct))
                if (epochs<10) or i%5 == 4:
                    print(f"epoch {i+1} ",
                          "- accuracy: {:.4f}".format(np.mean(e_correct)))
            else:
                losses.append(np.mean(e_results))
                if (epochs<10) or i%5 == 4:
                    print(f"epoch {i+1} ")
                
            
        results = np.array(losses)
        if count_correct:
            results = [results]
            results.append( np.array(correct))
        return results
            

    def test(self, test_x, test_y, loss_fn, count_correct = count_correct,
             correct_thr = 0.2, printcorrect = False):
        
        test_x, test_y = get_tensors(test_x, test_y, dtype = torch.float32)
        
        pred = self(test_x)
        if pred.size()[1] and len(test_y.size()) == 1:
            pred = pred[:,0]
        loss = loss_fn(pred, test_y)
        results = loss.item()
        if count_correct:
            correct_frac = count_correct(pred  .detach().numpy(), 
                                         test_y.detach().numpy(),
                                         correct_thr)
            results = [results]
            results.append(correct_frac)
            if printcorrect:
                print(f"Accuracy (to {correct_thr:.2f}): {correct_frac:.4f} ")
        return results
    
    
    def train_test (self, x, y, loss_fn, optimizer, 
                    test_frac = 0.2, batch_size = 100, epochs = 100,
                    count_correct = count_correct, correct_thr = 0.2):
        """
        Training over <epochs> epochs + testing. 

        Parameters
        ----------
        x : list, ndarray, DataFrame or tensor
            Complete x dataset
        y : list, ndarray, DataFrame or tensor
            Complete y dataset
        loss_fn : PyTorch loss function object
            Loss function to be minimised with backward() method. 
            It is for instance
                - loss_fn = nn.MSELoss()          #NOTE PARENTHESES
                - loss_fn = nn.CrossEntropyLoss() #NOTE PARENTHESES
        optimizer : PyTorch optimizer object
            Optimization method, with zero_grad() and step() methods
        test_frac : float in [0,1], default is 0.2
            Fraction of data to be used for testing
        batch_size : int, default 100
            Size of each batch
        epochs : int, default 100
            Number of epochs to loop over
        count_correct : function, default is helperML.count_correct()
            Function taking as argument prediction array, expected arrays, threshold
            and returning fraction of times the outputs result in correct 
            classification.
            If 0, instead, don't count.
        correct_thr : float in [0,1]
            Third argument of count_correct, it is the relative error between 
            results under which the prediction is considered correct.
            
        Returns
        -------
            test loss : float
                loss value of testing
            training losses : np.ndarray (epochs)
                loss value per epoch
            test correct : float
                correct fraction of testing, returned if count_correct != 0
            training accuracy : np.ndarray (epochs)
                correct fraction per epoch, returned if count_correct != 0  
        """
        train_x, test_x, train_y, test_y = train_test_split(x, y,
                                                 test_size = test_frac)
        train_x, train_y = get_tensors(train_x, train_y, dtype = torch.float32)
        
        start = timer()
        train_results = self.full_train(train_x, train_y, 
                                        loss_fn, optimizer, 
                                        batch_size, epochs,
                                        count_correct,
                                        correct_thr)
        end = timer()
        print("Training Time: {:.9f}s".format(end-start)); del start, end

        start = timer()
        test_results = self.test(test_x, test_y, loss_fn,
                                 count_correct = count_correct)
        end = timer()
        print("Testing Time: {:.9f}s".format(end-start)); del start, end
        
        if count_correct:
            train_losses  = train_results[0]
            train_correct = train_results[1]
            test_loss     = test_results [0]
            test_correct  = test_results [1]
            results = [test_loss   , train_losses, 
                       test_correct, train_correct]
        else:
            train_losses = train_results
            test_loss    = test_results
            results = [test_loss   , train_losses]
        
        print (f'LOSS    = {test_loss:.5f}')
        if count_correct:
            print (f'Accuracy ({correct_thr:.3}): {test_correct:.5f}')
        return results
    
    
    def plot_importance(self, nfeatures = 25, least = False, **kwargs):
        """
        Parameters
        ----------
        nfeatures : int, default is 25.
            Number of features to put on plot
        least : bool, default is False.
            Plot most important if False, plot least important if True.
        **kwargs : TYPE
            DESCRIPTION.
    
        Returns
        -------
        0
        """
        layer = self.fc0
        weights = []
        for i in layer:
            weights.append(sum(layer.weight[i]))
        inputs = np.arange(len(weights))
        weights, inputs = zip(*sorted(zip(weights, inputs))) 
        if least:
            weights = weights[:nfeatures]
            inputs = inputs  [:nfeatures]
        if not least:
            weights = weights[-nfeatures:]
            inputs = inputs  [-nfeatures:]
        pl.figure()
        pl.title("Inputs' Importance")
        pl.xlabel("Inputs' Neuron Number")
        pl.ylabel("Weight")
        return 0
    
    
    
    
def uniform_seq (n_in, n_out, n_hidden = [8,8], 
                 propfunc = nn.Linear, actfunc = nn.Sigmoid, prob = 0):
    """
    Generate nn.Sequential object with same propagation and activation 
    everywhere.
    The default values generate a Linear 3-layers Perceptron.

    Parameters
    ----------
    n_in : int
        Number of inputs
    n_out : int
        Number of outputs.
    n_hidden : list or ndarray, deafult is [4,4]
        len(n_hidden) = number of hidden layers
        n_hidden[i]   = number of neurons in ith hidden_layer
    propfunc : callable, default is torch.nn.Linear
        Forward propagation function
    actfunc : callable, default is torch.nn.Sigmoid
        Non-linear activation function
    prob : float in [0,1], default is 0
        Probability for torch.nn.Dropout after a layer info is propagated
    
    Return
    ---------
    seq_architecture: torch.nn.modules.container.Sequential
    """
    layer_list = []
    for nh_i in n_hidden:
        #layer_list.append(nn.BatchNorm1d(n_in))
        layer_list.append(propfunc(int(n_in), int(nh_i)))
        layer_list.append(nn.Dropout(p=prob))
        layer_list.append(actfunc())
        n_in = nh_i *1  
    layer_list.append(propfunc(int(n_hidden[-1]), int(n_out)))
    return nn.Sequential(*layer_list)



def uniform_norm_seq (n_in, n_out, n_hidden = [4,4], 
                      normplace = 0, norm_layer = nn.BatchNorm1d, 
                      propfunc = nn.Linear, actfunc = nn.Sigmoid, prob = 0):
    """
    Generate nn.Sequential object with same propagation and activation 
    everywhere.
    The default values generate a Linear 3-layers Perceptron.

    Parameters
    ----------
    n_in : int
        Number of inputs
    n_out : int
        Number of outputs.
    n_hidden : list or ndarray, deafult is [4,4]
        len(n_hidden) = number of hidden layers
        n_hidden[i]   = number of neurons in ith hidden_layer
    normplace : float or list
        after which layer put the normalisation: 0 means after input, i after
        ith hidden layer
    norm_layer : callable, default is torch.nn.BatchNorm1d
        Normalization layer type.
    propfunc : callable, default is torch.nn.Linear
        Forward propagation function
    actfunc : callable, default is torch.nn.Sigmoid
        Non-linear activation function
    prob : float in [0,1], default is 0
        Probability for torch.nn.Dropout after a layer info is propagated
    
    Return
    ---------
    seq_architecture: torch.nn.modules.container.Sequential
    """
    if not hasattr(normplace, "__len__"):
        normplace = [normplace]
    layer_list = []
    for i, nh_i in enumerate(n_hidden):
        if i in normplace:
            layer_list.append(norm_layer(int(n_in)))
        #layer_list.append(nn.BatchNorm1d(n_in))
        layer_list.append(propfunc(int(n_in), int(nh_i)))
        layer_list.append(nn.Dropout(p=prob))
        layer_list.append(actfunc())
        n_in = nh_i *1  
    layer_list.append(propfunc(int(n_hidden[-1]), int(n_out)))
    return nn.Sequential(*layer_list)
    
    
    
# class UniformNetwork(nn.Module, MyNetwork):
#     """
#     Generate a model of a Neural Network inheriting from torch.nn.Module, where
#     all layers are identical but in the number of neurons.

#     Parameters
#     ----------
#     n_in : int
#         Number of inputs
#     n_out : int
#         Number of outputs.
#     n_hidden : list or ndarray, deafult is [4,4]
#         len(n_hidden) = number of hidden layers
#         n_hidden[i]   = number of neurons in ith hidden_layer
#     propfunc : callable, default is torch.nn.Linear
#         Forward propagation function
#     actfunc : callable, default is torch.nn.Sigmoid
#         Non-linear activation function
#     actfunc : callable or 0, default is torch.sigmoid
#         Funcion operated by forward method before on outputs before
#         returning them. If 0, don't apply anything.
#     prob : float in [0,1], default is 0
#         Probability for torch.nn.Dropout after a layer info is propagated

#     """
#     def __init__(self, n_in, n_out, n_hidden = [4,4], 
#                  propfunc = nn.Linear, actfunc = nn.Sigmoid, 
#                  resultfunc = torch.sigmoid, prob = 0):

#         super().__init__()
#         self.flatten = nn.Flatten()
        
#         layer_list = []
#         for nh_i in n_hidden:
#             #layer_list.append(nn.BatchNorm1d(n_in))
#             layer_list.append(propfunc(int(n_in), int(nh_i)))
#             layer_list.append(nn.Dropout(p=prob))
#             layer_list.append(actfunc())
#             n_in = nh_i *1  
#         layer_list.append(propfunc(int(n_hidden[-1]), int(n_out)))
#         self.sequence = nn.Sequential(*layer_list)
        
#         self.resultfunc = resultfunc
    
    
    
                
    


