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
from sklearn.metrics import confusion_matrix, roc_curve, roc_auc_score
from utils import *
import sys

def main () :
    file = "BDT_scores.pkl"
    with open(file, 'rb') as f :
        scores = pickle.load (f) #scores è il dizionario BDT_prob_dict
        
    y_prob_bdt = scores['longitudinale']
    y_valid = scores['valid']
    #tras_scores = scores['trasversale']
    
    file = "DNN_scores.pkl"
    with open(file, 'rb') as f :
        dnn_scores = pickle.load(f)
        
    y_prob_dnn = dnn_scores['longitudinale']
    y_true = dnn_scores['true']
    
    fp_bdt, tp_bdt, _ = roc_curve(y_valid, y_prob_bdt)
    auc_bdt = roc_auc_score(y_valid, y_prob_bdt)
    fp_dnn, tp_dnn, _ = roc_curve(y_true, y_prob_dnn)
    auc_dnn = roc_auc_score(y_true, y_prob_dnn)
    
    file = "dataframes_z_long.pkl"
    with open(file, 'rb') as f:
        df_dict = pickle.load(f)
        
    variable_index = df_dict['names'].index("Pt Z")
    pt_values = df_dict['events'][:, variable_index]    
    
    fpr_dnn, tpr_dnn, _ = roc_curve(y_true, pt_values)
    fpr_bdt, tpr_bdt, _ = roc_curve(y_valid, pt_values)
    #auc_pt = auc(fpr, tpr)

    
    plt.figure()
    sb.set(style='whitegrid')
    plt.plot(fp_bdt, tp_bdt, color='royalblue', label= f"BDT: AUC = {auc_bdt:.3f}", linewidth = 2)
    plt.plot(fp_dnn, tp_dnn, color='firebrick', label= f"DNN: AUC = {auc_dnn:.3f}", linewidth = 2)
    plt.title('ROC curve - BDT vs DNN response')
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.legend()
    #plt.savefig("BDT_repsonse_score.png", dpi = 300)
    plt.show()
    
    plt.figure()
    sb.set(style='whitegrid')
    plt.plot(fp_bdt, tp_bdt, color='royalblue', label= f"BDT: AUC = {auc_bdt:.3f}", linewidth = 2)
    plt.plot(fp_dnn, tp_dnn, color='firebrick', label= f"DNN: AUC = {auc_dnn:.3f}", linewidth = 2)
    plt.title('ROC curve - BDT vs DNN response')
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.legend()
    #plt.savefig("BDT_repsonse_score.png", dpi = 300)
    plt.show()
    
    
    '''
    plt.figure()
    sb.histplot(y_prob_bdt, bins = 50)
    plt.show()
    '''
if __name__ == '__main__':
    main()  
           
