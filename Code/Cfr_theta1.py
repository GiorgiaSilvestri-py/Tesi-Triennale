import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils import *

def main () :
    '''
    Confronta le distribuzioni dell'angolo theta 1, definito come l'angolo tra il leptone e la direzione di volo dell'Higgs
    '''
    '''
    Procedimento:
    - SDR comovente con la Z
	    - trovo il quadrimomento della Z da quello dei suoi prodotti di decadimento
	    - calcolo beta come px/E, py/E etc con px, py, pz momento della Z
	    - fatto ciò, applico il boost (con -beta) al leptone
	    - ho ottenuto ora il quadrimomento dei due leptoni e il quadrimomento somma nel SDR della Z
    - applico il boost anche al bosone di Higgs -> ottengo il quadrimomento di h nel SDR della Z
    - a questo punto, calcolo l'angolo, tramite compute_angle, con p1 = p_leptone e p2 = p_Higgs
    '''
    
    #longitudinale
    LHE_L = "unweighted_events_L.lhe"
    eventsL_dict = read_file(LHE_L)
    
    #costruzione quadrimomenti
    W_fm_L = get_fm_of([24, -24], eventsL_dict)                                 #W
    Z_fm_L = get_fm_of([23], eventsL_dict)                                      #Z    
    H_fm_L = get_fm_of([25], eventsL_dict)                                      #H
    lep_ZL, antilep_ZL, lep_WL, antilep_WL = connect_lep_to_V(eventsL_dict)     #dizionari di leptoni e antileptoni da Z/W

    #trasversale
    LHE_T = "unweighted_events_T.lhe"
    eventsT_dict = read_file(LHE_T)

    W_fm_L = get_fm_of([24, -24], eventsL_dict)                                 #W
    Z_fm_T = get_fm_of([23], eventsT_dict)                                      #Z
    H_fm_T = get_fm_of([25], eventsT_dict)                                      #H
    lep_ZT,  antilep_ZT, lep_WT, antilep_WT = connect_lep_to_V(eventsT_dict)    #lep/antilep
    
    #liste di cos(theta)
    Z_cos_list_L = get_theta1_of(Z_fm_L, lep_ZL, antilep_ZL, H_fm_L)
    Z_cos_list_T = get_theta1_of(Z_fm_T, lep_ZT, antilep_ZT, H_fm_T)
    
       
    
    #visualizzazione distribuzioni
    data = Z_cos_list_L
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)
     
    sb.set(style="whitegrid")
    #fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    #plt.subplots_adjust(wspace=0.5)

    ax = sb.histplot(Z_cos_list_L, bins=nBins, color='royalblue', edgecolor = 'steelblue', stat='density', label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(Z_cos_list_T, bins=30, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, label='Polarizzazione trasversale')
    ax.set_xlabel("cosθ1")
    ax.set_ylabel("dN/N")
    ax.set_title("Distribuzione angolo θ1 (Z)")
    ax.legend()   
    #plt.savefig("Confronto coseno.png", dpi=300)
    plt.show()

  
if __name__ == '__main__' :
    main() 
