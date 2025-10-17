import configparser
import importlib
import sys
import numpy as np
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import seaborn as sb
import pandas as pd
import xgboost as xgb
import random
import pickle
import os

from torch.utils.data import Dataset, DataLoader  
from torch.utils.data import random_split
from torch.optim.lr_scheduler import ExponentialLR, LambdaLR                   #variazioni dei LR
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split



# to use local modules
import sys
sys.path.append ('/home/giorgia/MG/Risultati/ML/')

from learning_rate_annealers import expo_lr_scheduler, EarlyStopping
from TwoD_tabular_dataset import TabularDataset, compute_variable_importance, plot_roc_blocks
from TwoD_tabular_dataset import train_loop, test_loop
#from generator_torch import gen_noisy_khole, gen_noisy_ksphere
from spiral import init_weights
#from libs.utils import setup_logger, log_file_content
#from early_stopping import EarlyStopping                                 #ferma training


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


def read_samples (config, device):

    #legge training fraction e batch size dal file cfg (0.8, 512)
    train_frac   = config.getfloat ('training','train_frac', fallback = 0.5)                     #capitolo training cerca train_frac, altrimenti 0.5
    batch_size   = config.getint ('training','batch_size', fallback = 64)

    sample_files = [x.strip () for x in config.get ('samples', 'files').split (',')]            #legge il nome del file da cui prendere i dati
    sample_names = [x.strip () for x in config.get ('samples', 'names').split (',')]            #legge il nome associato ai due eventi

    if len (sample_files) > 2 :
        print ('expecting two sample files only to be compared')
        sys.exit (1)

    samples = {}
    variables = []
    
    for name, file in zip (sample_names, sample_files):             #apre file longitudinale/trasversale
        with open (file, 'rb') as f : 
            collection = pickle.load (f)
        sam = collection['events']
        nam = collection['names']
        samples[name] = sam             #sample
        variables.append (nam)          #nomi variabili
 
    for nam, sam in samples.items () :
        print('--- READING ---', nam, (sam.shape))

    all_equal = all (x == variables[0] for x in variables)
    if not all_equal :
        print ('not all files contain the same list of variables')
        sys.exit (1)
    
    #lista delle 18 variabili
    var_names = variables[0]                                            #variables è la lista di liste delle variabili in ciascuna tabella, prende [0] perché le variabili sono le stesse
    
    X_one = torch.from_numpy (samples['longitudinale']).float ()               #X = evento; converte a mat pyt
    y_one = torch.ones (X_one.shape[0], dtype=torch.long).to (device)          #y = 1 è etichetta per indicare longitudinale [y contiene solo 1]

    X_two = torch.from_numpy (samples['trasversale']).float ()          
    y_two = torch.zeros (X_two.shape[0], dtype=torch.long).to (device)          #y = 0 è etichetta per indicare trasversale [y contiene solo 0]

    X_total = torch.cat ((X_one, X_two))                                        #tensore con tutti i dati
    y_total = torch.cat ((y_one, y_two))                                       #tensore con tutte le etichette [avrà 0 e 1]

    indices = torch.randperm (X_total.shape [0])                        #rimescolamento degli indici -> crea permutazione di indici e la applica sia a segnale che fondo
    X_total = X_total[indices]
    y_total = y_total[indices]

    input_dataset = TabularDataset (X_total, y_total)                    
    print('Sample with ' + str(len(input_dataset)) + ' elements generated')

    N_test = int (len(input_dataset) * train_frac)                                                  # numero eventi test
    N_val = len (input_dataset) - N_test                                                            # numero eventi validation
    train_set, valid_set = random_split (input_dataset, [N_test, N_val])                          #divisione casuale del dataset in train_set e valid_set
    # train_set, valid_set, test_set = random_split (dataset, [XXX, YYY, ZZZ])
    
    X_np = X_total.cpu().numpy()
    y_np = y_total.cpu().numpy()

    print("---- IMPLEMENTING BDT ----")
    X_train, X_valid, y_train, y_valid = train_test_split(X_np, y_np, test_size=0.2, random_state=11)

    model_xgb = xgb.XGBClassifier(tree_method="hist", enable_categorical=True, device="cpu")
    
    
    model_xgb.fit(X_train, y_train)
    y_pred = model_xgb.predict(X_valid)
    y_prob = model_xgb.predict_proba(X_valid)[:, 1]  #array di probabilità per ogni evento che sia 1 (longitudinale)
    y_prob_t = model_xgb.predict_proba(X_valid)[:, 0]  #array di probabilità per ogni evento che sia 0 (trasversale)
    
    # Confusion matrix
    cm = confusion_matrix(y_valid, y_pred)
    sb.heatmap(cm, annot=True, fmt='d', cmap='Reds')
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.title('Confusion Matrix - XGBoost')
    plt.savefig("BDT_confusion_matrix.png", dpi = 300)
    #plt.show()

    # Distribuzione del BDT response score
    
    #crea distribuzione di di prob (che sia longitudinale) solo per eventi 1
    long_scores = [p for p, y in zip(y_prob, y_valid) if y == 1] 
    #crea distribuzione di prob (che sia longitudinale) sono per eventi 0
    tras_scores = [p for p, y in zip(y_prob, y_valid) if y == 0] 

    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(long_scores, bins=40, color='royalblue', stat='probability', label='Longitudinale', alpha=0.5)
    sb.histplot(tras_scores, bins=40, color='firebrick', stat='probability', label='Trasversale', alpha=0.5)
    plt.title('BDT Response Score')
    plt.xlabel('Score')
    plt.ylabel('Event Fraction')
    plt.legend()
    plt.savefig("BDT_repsonse_score.png", dpi = 300)
    #plt.show()
    
    print("--- CREATING TXT FILE ---")
    print("--- ADDING BDT SCORES TO FILE ---")
    
    #creazione dizionario di probabilità
    BDT_prob_dict = {
               "longitudinale": y_prob,
               "trasversale"  : y_prob_t,
               "valid"        : y_valid
               }
   
    #creazione del dataframe
    df = pd.DataFrame(BDT_prob_dict)
            
    #sovrascrittura file
    with open("BDT_scores.pkl", "wb") as f:
        pickle.dump(df, f)
        
    print("---- file created successfully ----")


    
    print('Training dataset with ' + str(len(train_set)) + ' elements generated')
    print('Validation dataset with ' + str(len(valid_set)) + ' elements generated')

    #creazione dataloader
    train_dataloader = DataLoader (train_set, batch_size = batch_size, shuffle = True)     
    valid_dataloader = DataLoader (valid_set, batch_size = batch_size, shuffle = True)

    return train_dataloader, valid_dataloader, var_names


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


def models_training (config, train_dataloader, valid_dataloader, device, verbosity, model_class) :

    epochs               = config.getint ('training','epochs', fallback = 1000)
    layer_norm           = config.getboolean ('training','layer_norm', fallback = True)
    learning_rate        = config.getfloat ('training','learning_rate', fallback = 1e-3)
    lr_annealing         = config.getfloat ('training','lr_annealing', fallback = 0.)
    lr_minimum           = config.getfloat ('training','lr_minimum', fallback = 1.)
    lr_delay             = config.getfloat ('training','lr_delay', fallback = 1.)
    weight_decay         = config.getfloat ('training','weight_decay', fallback = 1e-5)   #evita overtraining e ridurre importanza di pesi troppo alti
    dropout_prob         = config.getfloat ('training','dropout_prob', fallback = 0.)
    dropout_annealing    = config.getfloat ('training','dropout_annealing', fallback = 0.) #riduce il drouput all'aumentare del tempo
    if dropout_annealing >= 1. : dropout_annealing = -1
    out_dim              = config.getint ('model', 'classes_num')
    hidden_layers        = config.getint ('training','hidden_layers', fallback = 3)
    hidden_units         = config.getint ('training','hidden_units', fallback = 2)
    save_model           = config.getboolean ('saving','save_model', fallback = True)

    apply_early_stopping = config.getboolean ('training','apply_early_stopping', fallback = True)
    es_patience          = config.getint ('training','es_patience', fallback = 10)  #verifica che la loss sia effettivamente al minimo
    es_starting          = config.getint ('training','es_starting', fallback = 10)  #early stopping aspetta un po'
    
    
    log_names = [
    "epochs", "layer_norm", "learning_rate", "lr_annealing", "lr_minimum", "lr_delay",
    "weight_decay", "dropout_prob", "dropout_annealing", "out_dim", "hidden_layers",
    "save_model", "apply_early_stopping", "es_patience", "es_starting"
    ]

    log_values = [
        epochs, layer_norm, learning_rate, lr_annealing, lr_minimum, lr_delay,
        weight_decay, dropout_prob, dropout_annealing, out_dim, hidden_layers,
        save_model, apply_early_stopping, es_patience, es_starting
    ]

    for name, value in zip(log_names, log_values):
        print(f"{name}: {value}")


    

    out_path             = config.get ('saving', 'out_path', fallback = None)
    out_file_base        = config.get ('saving', 'out_file_base', fallback = None)

    loss_fn = nn.CrossEntropyLoss (reduction ='mean')           #metrica per valutare il modello: viene minimizzata dal training [è misura di distanza]
    # reduction ='sum' to get the sum of losses in the whole sample

    # get the number of features per entry in the sample
    in_dim = 0                                              
    for batch in train_dataloader:                      #prende il set, legge il batch e trova la dimensione del batch in input
        inputs, targets = batch
        in_dim = inputs.shape[1]                        #numero di colonne per la prima riga -> dice numero di variabili nel campione
        break                                           # Only need the first batch

    hidden_units = in_dim  # overriding the cfg file here

    print('--- TRAINING ---', hidden_layers, hidden_units)
    # ---------------------------------------

    model = model_class (in_dim, hidden_units, hidden_layers, out_dim, dropout_prob, layer_norm).to (device)  
    print(str(model))       

    optimizer = torch.optim.Adam (model.parameters (), lr = learning_rate, weight_decay = weight_decay)      #dice come fare la minimizzazione della funzione di loss 
    model.apply (init_weights)    #inizializza pesi

    train_loss_acc = []
    valid_loss_acc = []

    # initialize Early Stopping
    early_stopping = EarlyStopping (patience = es_patience, verbose = True)

    # initialize learning rate annealing
    # scheduler = ExponentialLR (optimizer, gamma = lr_annealing)
    lr_annealer = lambda e : expo_lr_scheduler (e, learning_rate, lr_minimum, lr_delay)
    scheduler = LambdaLR (optimizer, lr_lambda = lr_annealer)

    for t in range (epochs):                #loop su epoche
        if t%100 == 0 : print(f' - epoch {t+1}')
        train_loop (train_dataloader, model, loss_fn, optimizer, device)
        train_loss_acc.append (test_loop (train_dataloader, model, loss_fn, 0))  #confronto la loss con validation e con train per capire se siamo in overtraining
        valid_loss_acc.append (test_loop (valid_dataloader, model, loss_fn, 0))

        if dropout_annealing > 0. : model.scale_do (dropout_annealing)
        if lr_annealing > 0       : scheduler.step ()

        if not apply_early_stopping or t < es_starting : continue
        early_stopping (valid_loss_acc[-1][0], model, t)
        if early_stopping.early_stop:
            # print("Early stopping triggered!")
            print(f'EARLY STOPPING triggered at epoch {t} with best epoch {early_stopping.best_epoch}\n')

            # riempi le collezioni con l'ultimo numero inserito fino alla durata epochs
            while len (train_loss_acc) < epochs: train_loss_acc.append (train_loss_acc[-1])
            while len (valid_loss_acc) < epochs: valid_loss_acc.append (valid_loss_acc[-1])
            break  #esce dal ciclo se la funzione di early stopping è soddisfatta dal training

    print('training finished')

    saveObject = (train_loss_acc, valid_loss_acc)
    out_file = out_file_base + '_trends.pickle'
    with open (out_path + out_file,'wb') as f:  #salva i trend
        pickle.dump (saveObject, f)
    print(' loss and accuracy trends saved in ', out_path, out_file)

    if save_model :             #salva modello
        out_file = out_file_base + '_model.pth'
        torch.save (early_stopping.best_model_state_dict, out_path + out_file)
        print(' model saved in ', out_path, out_file)

    return model
    

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


def main () :

    if (len (sys.argv) < 2):
        print ('usage:', sys.argv[0], 'input.cfg')
        sys.exit (1)

    # reading the config file
    # -----------------------

    config = configparser.ConfigParser ()
    config.read ('test_MLP_11.cfg')

    out_path  = config.get ('saving', 'out_path', fallback = None)
    print("out_path:", out_path)
    
    print("Sezioni trovate:", config.sections())
    print("Contenuto saving:", dict(config['saving']) if 'saving' in config else "Non trovata")
    #os.makedirs (out_path, exist_ok = True)     #create directory when non-existing

    with open (out_path + 'config.cfg', 'w') as configfile :
        config.write (configfile)

    out_log   = config.get ('saving', 'out_log', fallback = None)
    verbosity = config.getint ('running','verbosity', fallback = 20)

    #setup_logger (logger, out_path + out_log, verbosity)
    #log_file_content (logger, sys.argv[1])

    print('\n --- PREPARING SAMPLES ---\n')
    # ---------------------------------------

    seed = config.getint ('training','seed', fallback = -1)
    if seed > 0 : 
        random.seed (seed)
        np.random.seed (seed)
        torch.manual_seed (seed)

    device = torch.device ("cuda:0" if torch.cuda.is_available () else "cpu") #tecnica
    #logger.debug (f"Using {device} device")


    #creazione dei dataset di traning e validazione
    train_dataloader, valid_dataloader, var_names = read_samples (config, device)
    
    print('--- IMPORTING THE MODEL ---')
    # ---------------------------------------

    library_name  = config.get ('model', 'library_name')
    model_name    = config.get ('model', 'model_name', fallback=None)

    print('Using library: ', library_name)
    print('Using model: ', model_name)
    model_library = importlib.import_module (library_name)
    model_class = getattr (model_library, model_name)

    model = models_training (config, train_dataloader, valid_dataloader, device, verbosity, model_class)
    
    # --- CONFUSION MATRIX ---
    print("\n--- EVALUATING MODEL ---")
    
    in_dim = 0                                              
    for batch in train_dataloader:                      #prende il set, legge il batch e trova la dimensione del batch in input
        inputs, targets = batch
        in_dim = inputs.shape[1]                        #numero di colonne per la prima riga -> dice numero di variabili nel campione
        break                                           # Only need the first batch


    #legge file salvato come miglior modello (se early stopping)
    out_path = config.get ('saving', 'out_path')
    out_file_base = config.get ('saving', 'out_file_base')
    #model = models_training (config, train_dataloader, valid_dataloader, device, verbosity, model_class)
    model.load_state_dict(torch.load(out_path + out_file_base + '_model.pth'))
    
    y_true, y_pred = [], []
    dnn_values = []

    model.eval()
    with torch.no_grad():
        for inputs, targets in valid_dataloader:
            inputs = inputs.to(device)
            targets = targets.to(device)
            outputs = model(inputs)
            preds = torch.argmax(outputs, dim=1)
            probs = torch.softmax(outputs, dim=1)[:, 1]  #prob per ogni evento che sia 1 (long)
            probs_t = torch.softmax(outputs, dim=1)[:, 0]  #probabilità per ogni evento che sia 0 (tr)
            dnn_values.extend(probs.cpu().numpy())
            y_true.extend(targets.cpu().numpy())
            y_pred.extend(preds.cpu().numpy())

    #matrice di confusione
    class_names = ['trasversale', 'longitudinale']  # 0, trasversale; 1, longitudinale
    
    cm = confusion_matrix(y_true, y_pred, labels = [0, 1])
    #sb.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=class_names, yticklabels=class_names)
    sb.heatmap(cm, annot=True, fmt='d', cmap='Reds', xticklabels=class_names, yticklabels=class_names)
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.title('Confusion Matrix')
    plt.savefig("Confusion_matrix.png", dpi = 300)
    
    print(len(y_true), len(dnn_values), len(probs_t))
   
    #DNN response
    long_scores = [l for l, y in zip(dnn_values, y_true) if y == 1]
    tras_scores = [t for t, y in zip(dnn_values, y_true) if y == 0]
    
    #print("Lunghezza DNN: ", len(bdt_values), len(y_true))
    print("Range DNN: ", min(dnn_values), max(dnn_values))
    #print("Lunghezza scores:", len(long_scores), len(tras_scores))
    
    #costruzione di dizionario di probabilità
    DNN_prob_dict = {
               "longitudinale": dnn_values,
               "true"         : y_true
               }
               
    dnn_df = pd.DataFrame(DNN_prob_dict)
    
    with open("DNN_scores.pkl", "wb") as file:
        pickle.dump(dnn_df, file)

    
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(long_scores, bins = 40, color = 'royalblue', stat='probability', edgecolor = 'steelblue', element='step', linewidth=1.5, alpha = 0.4, label='Polarizzazione longitudinale')
    sb.histplot(tras_scores, bins = 40, color='firebrick', stat='probability', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, label='Polarizzazione trasversale')
    
    plt.title('DNN distribution')
    plt.xlabel('DNN response score')
    plt.ylabel('Event fraction')    
    plt.savefig("DNN_distribution.png", dpi = 300)
    plt.show()
    
    
    log_variables, log_roc = compute_variable_importance(model, valid_dataloader, variables_name = var_names, device=device)
    
    print("\n---- VARIABLES RELEVANCE ----")
    print("\n - Variables score drop and standard deviation -")


    print(log_variables) 
    plot_roc_blocks(log_roc, block_size=9, filename_prefix="ROC_curves")

    '''
    for log in log_roc:
        name = log["variable"]
        
        fp, tp, auc = log["feature"]["fp"], log["feature"]["tp"], log["feature"]["auc"]
        fp_dnn, tp_dnn, auc_dnn = log["dnn"]["fp"], log["dnn"]["tp"], log["dnn"]["auc"]
        fp_bdt, tp_bdt, auc_bdt = log["bdt"]["fp"], log["bdt"]["tp"], log["bdt"]["auc"]
        
        plt.figure()
        sb.set(style='whitegrid')
        
        plt.plot(fp, tp, color='royalblue', label=f"{name}, AUC = {auc:.3f}", linewidth=2)
        plt.plot(fp_dnn, tp_dnn, color='firebrick', label=f"DNN: {name}, AUC = {auc_dnn:.3f}", linewidth=2)
        plt.plot(fp_bdt, tp_bdt, color='darkgreen', label=f"BDT: {name}, AUC = {auc_bdt:.3f}", linewidth=2)
        
        plt.title(f'ROC Curve - {name}')
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        plt.legend()
        plt.show()
    '''
        
        
   


    
    
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


if __name__ == "__main__" : 
    main ()
    




