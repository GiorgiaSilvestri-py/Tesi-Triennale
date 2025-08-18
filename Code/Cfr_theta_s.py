import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils import *

def main () :
    '''
    Confronta le distribuzioni dell'angolo theta*, definito come la direzione di volo di Z (o W) calcolato nel sdr del centro di massa VH
    ''' 
    
    #longitudinale
    LHE_L = "unweighted_events_L.lhe"
    eventsL_dict = read_file(LHE_L)
    
    W_fm_L = get_fm_of(24, eventsL_dict)               #dizionario di quadrimomenti con ID evento per la W
    Z_fm_L = get_fm_of(23, eventsL_dict)               #dizionario di quadrimomenti con ID evento per la Z
    
    H_fm_L = get_fm_of(25, eventsL_dict)               #dizionario di quadrimomenti con ID evento per H
    #compute_tot_fm(Z_fm, H_fm)                      #dizionario di quadrimomenti totali (ZH) con ID evento
    
    theta_star_ZL = get_thetastar_of(Z_fm_L, H_fm_L)
    theta_star_WL = get_thetastar_of(W_fm_L, H_fm_L)
 
    #trasversale
    LHE_T = "unweighted_events_T.lhe"
    eventsT_dict = read_file(LHE_T)

    W_fm_T = get_fm_of(24, eventsT_dict)               #dizionario di quadrimomenti con ID evento per la W
    Z_fm_T = get_fm_of(23, eventsT_dict)               #dizionario di quadrimomenti con ID evento per la Z
    
    H_fm_T = get_fm_of(25, eventsT_dict)               #dizionario di quadrimomenti con ID evento per la Z
    #fm_tot = compute_tot_fm(Z_fm_T, H_fm_T)                    #dizionario di quadrimomenti totali (ZH) con ID evento
    
    
    #calcolo l'angolo
    theta_star_ZT = get_thetastar_of(Z_fm_T, H_fm_T)
    theta_star_WT = get_thetastar_of(W_fm_T, H_fm_T)
    
        
    #visualizzazione distribuzioni
    data = theta_star_ZL
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)
    
    data = theta_star_WL
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_bins = int((max(data) - min(data)) / bin_width)
    
     
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(theta_star_ZL, bins=40, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(theta_star_ZT, bins=40, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversale')
    ax[0].set_xlabel("cosθ* (rad)")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione angolo θ* (Z)")
    ax[0].legend()
    
    sb.histplot(theta_star_WL, bins=40, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(theta_star_WT, bins=40, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversale')
    ax[1].set_xlabel("cosθ* (rad)")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione angolo θ* (W)")
    ax[1].legend()

    #plt.savefig("Confronto massa invariante.png", dpi=300)
    plt.show()
  
if __name__ == '__main__' :
    main() 
