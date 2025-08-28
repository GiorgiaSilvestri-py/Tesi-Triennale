import numpy as np
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import math
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve
import seaborn as sb
import xgboost as xgb
import configparser
from model_MLP_06 import MLP_06
from torch.utils.data import Dataset, DataLoader

from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt

'''
Torch Module documentation:
https://pytorch.org/docs/stable/notes/modules.html


A custom Dataset class must implement three functions: 
  - **__init__**:
    The __init__ function is run once when instantiating the Dataset object
  - **__len__**:
    The __len__ function returns the number of samples in our dataset
  - **__getitem__**:
    The __getitem__ function loads a sample,
    transforms it if needed,
    and returns the tensor image and corresponding label in a tuple
'''

class TabularDataset (Dataset):
    """Torch container for a spiral dataset."""

    def __init__(self, X, y):
        """Initializes instance of class StudentsPerformanceDataset.
        Args:
            input_collection: sample randomly generated with generate_spiral_data
        """
        # Save target and predictors
        self.values = X
        self.labels = y

    def __len__ (self) :
        return len (self.values)

    def __getitem__ (self, idx) :
        return [self.values[idx], self.labels[idx]]


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

#nel training vogliamo applicare la back propagation e il dropout serve solo nel training -> training e test si comportano in modo diverso



def train_loop (dataloader, model, loss_fn, optimizer, device, verbosity = 0):
    ''' 
    the function to be called to train the model for one epoch,
    looping over serval batches
    '''
    size = len (dataloader.dataset)
    # Set the model to training mode - 
    # important for batch normalization and dropout layers
    # Unnecessary in this situation but added for best practices
    model.train ()
    loss = -1.
    current = -1.
    # loop on batches
    for batch, (X, y) in enumerate (dataloader):      
        X, y = X.to (device), y.to (device)

        # Compute prediction error
        pred = model (X)
        loss = loss_fn (pred, y)       #calcola funzione di loss 

        # Backpropagation
        loss.backward ()		#propaga info della loss nel modello
        optimizer.step ()
        optimizer.zero_grad ()

        if verbosity != 0 and batch % 100 == 0:
            loss, current = loss.item (), (batch + 1) * len (X)
            print (f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

#non ha più backprop

def test_loop (dataloader, model, loss_fn, verbosity = 1):
    ''' 
    the function to be called to test the model
    verbosity == 1 : talkative
    verbosity == 0 : silent 
    '''
    # Set the model to evaluation mode - 
    # important for batch normalization and dropout layers
    # Unnecessary in this situation but added for best practices
    # print(f"Is the model in training mode after calling model.eval()? {model.training}")
    model.eval ()
    size = len (dataloader.dataset)
    num_batches = len (dataloader)
    test_loss, correct = 0., 0.
    TP, TN, FP, FN = 0., 0., 0., 0.

    # Evaluating the model with torch.no_grad() ensures that 
    # no gradients are computed during test mode
    # also serves to reduce unnecessary gradient computations 
    # and memory usage for tensors with requires_grad=True
    with torch.no_grad () :
        # in this loop each entry is a batch and X.size (0) its size
        for X, y in dataloader :
            pred = model (X)

            # print (' testing: ', len (y), len (X), len (pred))
            test_loss += loss_fn (pred, y).item () * X.size (0)

            TP += ((pred.argmax (1) == 1) & (y == 1)).sum ().item () #true positive etc
            TN += ((pred.argmax (1) == 0) & (y == 0)).sum ().item ()
            FP += ((pred.argmax (1) == 1) & (y == 0)).sum ().item ()
            FN += ((pred.argmax (1) == 0) & (y == 1)).sum ().item ()

            # using y as a mask (~ negates it) 
            # leaving float instead of int for cases where weights may be used
            # TP += (pred.argmax (1) == y).type (torch.float)[y.bool ()].sum ().item ()		
            # FN += (pred.argmax (1) != y).type (torch.float)[y.bool ()].sum ().item ()
            # TN += (pred.argmax (1) == y).type (torch.float)[~y.bool ()].sum ().item ()
            # FP += (pred.argmax (1) != y).type (torch.float)[~y.bool ()].sum ().item ()
            # correct += (pred.argmax (1) == y).type (torch.float).sum ().item ()

    correct = TP + TN
    test_loss /= size
    accuracy = correct / size				                    #quante volte predizione = atteso
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0          #quante previsioni positive sono corrette
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0	            #quanti veri positivi sono letti
    F1 = 2 * (precision * recall) / (precision + recall) if TP > 0 else 0 #media geom*2

    if verbosity > 0 : 
        # print (f"Test Error: \n Accuracy: {(100*correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")
        print ('Test metrics:')
        print (f' - Accuracy : {(100*accuracy):>0.1f}%')
        print (f' - Avg loss : {test_loss:>8f}')
        print (f' - Precision: {precision:>8f}')
        print (f' - Recall   : {recall:>8f}')
        print (f' - F1       : {F1:>8f}')

    return test_loss, accuracy, precision, recall, F1
    

# ------------------------------------------------------------------------------------------


def compute_variable_importance(model,  dataloader, variables_name, metric = roc_auc_score, N=10, device = 'cpu'):


    print("Variabili: ", variables_name)
    
    X_total = []
    y_total = []
    roc_summary = []
    
    for X, y in dataloader:
        print(dataloader)
        X_total.append(X)
        y_total.append(y)
        
        X = torch.cat(X_total, dim=0).to(device)
        y = torch.cat(y_total, dim=0).cpu().numpy()
        
        with torch.no_grad():
            pred = model(X)
            prob = torch.softmax(pred, dim=1)[:, 1].cpu().numpy()
            in_score = metric(y, prob)                              #initial variable roc score
            
        relevance = []
        X_np = X.cpu().numpy()
        
    for i in range(X_np.shape[1]):

        score_list = []
    
        
        for n in range(N):
            X_perm = X_np.copy()
            np.random.shuffle(X_perm[:, i])
            X_tensor = torch.from_numpy(X_perm).float().to(device)
            
            with torch.no_grad():
                pred = model(X_tensor)
                prob = torch.softmax(pred, dim=1)[:, 1].cpu().numpy()
                perm_score = metric(y, prob)
                score_list.append(perm_score)
        
        name = variables_name[i]
        metric_difference = in_score - np.mean(score_list)
        st_dev = np.std(score_list)
        relevance.append( {
                            "Variabile": name, 
                            "Drop AUC": metric_difference, 
                            "Deviazione Std": st_dev
                            })
    
    
    log_roc = []
    
    
    for i in range(X_np.shape[1]):
        name = variables_name[i]
        X_s = (X_np.copy())
        X_single = X_s[:, i]
        
        fp, tp, _ = roc_curve(y, X_single)
        roc_score = metric(y, X_single)
        
        X_s_tensor = torch.from_numpy(X_s).float().to(device)
        
        with torch.no_grad():
            pred = model(X_s_tensor)
            prob = torch.softmax(pred, dim=1)[:, 1].cpu().numpy()
            fp_dnn, tp_dnn, _ = roc_curve(y, prob)
            auc_dnn = metric(y, prob)
            
        
        X_train, X_valid, y_train, y_valid = train_test_split(X_s, y, test_size=0.2, random_state=11)

        model_xgb = xgb.XGBClassifier(tree_method="hist", enable_categorical=True, device="cpu")
        model_xgb.fit(X_train, y_train)
        y_pred = model_xgb.predict(X_valid)
        y_prob_bdt = model_xgb.predict_proba(X_valid)[:, 1]  
        fp_bdt, tp_bdt, _ = roc_curve(y_valid, y_prob_bdt)
        auc_bdt = metric(y_valid, y_prob_bdt)
        
        
        
        log_roc.append({
                        "variable": name,
                        "feature": {
                            "fp": fp,
                            "tp": tp,
                            "auc": roc_score
                        },
                        "dnn": {
                            "fp": fp_dnn,
                            "tp": tp_dnn,
                            "auc": auc_dnn
                        },
                        "bdt": {
                            "fp": fp_bdt,
                            "tp": tp_bdt,
                            "auc": auc_bdt
                        }
                    })
        '''
        plt.figure()
        sb.set(style='whitegrid')
        plt.plot(fp, tp, color='royalblue', label= f"{name}, AUC = {roc_score:.3f}", linewidth = 2)
        plt.plot(fp_dnn, tp_dnn, color='firebrick', label= f"DNN: {name}, AUC = {auc_dnn:.3f}", linewidth = 2)
        plt.plot(fp_bdt, tp_bdt, color='firebrick', label= f"BDT: {name}, AUC = {auc_bdt:.3f}", linewidth = 2)
        plt.title('ROC curve')
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        plt.legend()
        #plt.savefig("BDT_repsonse_score.png", dpi = 300)
        plt.show()
    
            
        df_relevance = pd.DataFrame(relevance).sort_values(by="Drop AUC", ascending=False)
    
    return df_relevance
    '''        
    df = pd.DataFrame(relevance)
    df = df.sort_values(by = "Drop AUC", ascending = False)
            
    return df, log_roc


      
# ------------------------------------------------------------------------------------------------------------------------------------------     
      
def plot_roc_blocks(results, block_size=9, filename_prefix="ROC_block"):
    n = len(results)
    n_blocks = math.ceil(n / block_size)
    
    for block_idx in range(n_blocks):
        start = block_idx * block_size
        end = min(start + block_size, n)
        subset = results[start:end]
        
        fig, axs = plt.subplots(3, 3, figsize=(18, 15))  # 3x3 grid for 9 plots
        axs = axs.flatten()
        
        for ax_idx, res in enumerate(subset):
            name = res["variable"]
            fp, tp, auc = res["feature"]["fp"], res["feature"]["tp"], res["feature"]["auc"]
            fp_dnn, tp_dnn, auc_dnn = res["dnn"]["fp"], res["dnn"]["tp"], res["dnn"]["auc"]
            fp_bdt, tp_bdt, auc_bdt = res["bdt"]["fp"], res["bdt"]["tp"], res["bdt"]["auc"]
            
            sb.set(style='whitegrid')
            ax = axs[ax_idx]
            
            ax.plot(fp, tp, color='royalblue', label=f"{name}, AUC={auc:.3f}", linewidth=2)
            ax.plot(fp_dnn, tp_dnn, color='firebrick', label=f"DNN, AUC={auc_dnn:.3f}", linewidth=2)
            ax.plot(fp_bdt, tp_bdt, color='darkgreen', label=f"BDT, AUC={auc_bdt:.3f}", linewidth=2)
            
            ax.set_title(name)
            ax.set_xlabel('False positive rate')
            ax.set_ylabel('True positive rate')
            ax.legend(fontsize='small')
            ax.grid(True)
        
        
        
        plt.tight_layout()
        plt.savefig(f"{filename_prefix}_{block_idx+1}.png", dpi=300)
        plt.close()
