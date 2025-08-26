import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
import pandas as pd
import configparser
import seaborn as sb
import pickle
from scipy.stats import wasserstein_distance, ks_2samp
from utils import *
import sys

def main () :
    file = "test_MLP_11_09_trends.pickle"
    with open(file, 'rb') as f :
        trends = pickle.load (f)
        
    config = configparser.ConfigParser ()
    config.read ('test_MLP_11.cfg')

    epoch_number = config.getint('training', 'epochs')
    print(epoch_number)
    
    #train_loss_acc e valid_loss_acc sono liste [(loss, accuracy, precision, recall, F1), (...)]x   
    train_loss = [t[0] for t in trends[0]]
    valid_loss = [t[0] for t in trends[1]]
    
    train_acc = [t[1] for t in trends[0]]
    valid_acc = [t[1] for t in trends[1]]
    
    train_precision = [t[2] for t in trends[0]]
    valid_precision = [t[2] for t in trends[1]]
    
    train_F1 = [t[3] for t in trends[0]]
    valid_F1 = [t[3] for t in trends[1]]
    
    x_sample = np.linspace(0, epoch_number, len(train_loss))
    
    sb.set(style = "whitegrid")
    fig, ax = plt.subplots(2, 2, figsize = (16, 5))
    
    sb.lineplot(x = x_sample, y = train_loss, color='royalblue', label = 'Training loss value', ax = ax[0, 0], linewidth=2)
    sb.lineplot(x = x_sample, y = valid_loss, color='firebrick', label = 'Validation loss value', ax = ax[0, 0], linewidth=2)
    ax[0, 0].set_xlabel('Epoch number')
    ax[0, 0].set_ylabel('Loss value')
    ax[0, 0].legend()
    
    sb.lineplot(x = x_sample, y = train_acc, color='royalblue', label = 'Training accuracy value', ax = ax[0, 1], linewidth=2)
    sb.lineplot(x = x_sample, y = valid_acc, color='firebrick', label = 'Validation accuracy value', ax = ax[0, 1], linewidth=2)
    ax[0, 1].set_xlabel('Epoch number')
    ax[0, 1].set_ylabel('Accuracy value')
    ax[0, 1].legend()
    
    sb.lineplot(x = x_sample, y = train_precision, color='royalblue', label = 'Training precision value', ax = ax[1, 0], linewidth=2)
    sb.lineplot(x = x_sample, y = valid_precision, color='firebrick', label = 'Validation precision value', ax = ax[1, 0], linewidth=2)
    ax[1, 0].set_xlabel('Epoch number')
    ax[1, 0].set_ylabel('Precision value')
    ax[1, 0].legend()
    
    sb.lineplot(x = x_sample, y = train_F1, color='royalblue', label = 'Training F1 value', ax = ax[1, 1], linewidth=2)
    sb.lineplot(x = x_sample, y = valid_F1, color='firebrick', label = 'Validation F1 value', ax = ax[1, 1], linewidth=2)
    ax[1, 1].set_xlabel('Epoch number')
    ax[1, 1].set_ylabel('F1 value')
    ax[1, 1].legend()
       
    
    plt.savefig("loss_func_graph.png", dpi=300)
    plt.show()

if __name__ == '__main__':
    main()  
           
