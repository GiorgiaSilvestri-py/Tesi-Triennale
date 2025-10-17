import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils import *


    
            
def main():
    '''
    Confronto tramite pt Z e W + grafico scatter
    '''
    apply_smearing = True
    print("---- BUILDING DICTIONARIES ----")
    #analisi polarizzazione longitudinale
    LHE_L = "unweighted_events_L.lhe"
    eventsL_dict = read_file(LHE_L)
    
    #costruzione quadrimomenti
    W_fm_L = get_fm_of([24, -24], eventsL_dict, apply_smearing = apply_smearing)                                 #W
    Z_fm_L = get_fm_of([23], eventsL_dict , apply_smearing = apply_smearing)                                      #Z    
    H_fm_L = get_fm_of([25], eventsL_dict, apply_smearing = apply_smearing)                                      #H
    zlep_l, zantilep_l, wlep_l, wantilep_l = build_decay_products(eventsL_dict, apply_smearing = apply_smearing)     #dizionari di leptoni e antileptoni da Z/W
    
    
    #analisi polarizzazione trasversale
    LHE_T = "unweighted_events_T.lhe"
    eventsT_dict = read_file(LHE_T)
    
    #quadrimomenti
    W_fm_T = get_fm_of([24, -24], eventsT_dict , apply_smearing = apply_smearing)               
    Z_fm_T = get_fm_of([23], eventsT_dict, apply_smearing = apply_smearing)               
    zlep_t, zantilep_t, wlep_t, wantilep_t = build_decay_products(eventsT_dict, apply_smearing = apply_smearing)
    
    
    ptb_l, ptb_lw = get_lep_ptbalance(eventsL_dict, apply_smearing = apply_smearing)

    ptb_t, ptb_tw = get_lep_ptbalance(eventsT_dict, apply_smearing = apply_smearing)
       
    print("---- BUILDING HISTOGRAMS ----")
    data = ptb_l
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)

    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(ptb_l, bins=35, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(ptb_t, bins=35, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversale')
    ax[0].set_xlabel("Pt balance")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione del lepton pt balance (Z)")
    ax[0].legend()
    
    sb.histplot(ptb_lw, bins=40, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(ptb_tw, bins=40, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversale')
    ax[1].set_xlabel("Pt balance")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione del lepton pt balance (W)")
    ax[1].legend()

    plt.savefig("Confronto lepton pt balance.png", dpi=300)
    plt.show()
    
        
if __name__ == '__main__':
    main()  
                    
