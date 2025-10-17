import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils import *

def main () :
    '''
    Confronta le distribuzioni dell'angolo delta phi, definito come la differenza di angolo azimutale tra i leptoni prodotti dal decadimento nel SDR VH
    ''' 
    
    print("---- BUILDING DICTIONARIES ----")
    apply_smearing = True
    #longitudinale
    LHE_L = "unweighted_events_L.lhe"
    eventsL_dict = read_file(LHE_L)
    
    #costruzione quadrimomenti
    W_fm_L = get_fm_of([24, -24], eventsL_dict, apply_smearing = apply_smearing)                                 #W
    Z_fm_L = get_fm_of([23], eventsL_dict, apply_smearing = apply_smearing)                                      #Z    
    H_fm_L = get_fm_of([25], eventsL_dict, apply_smearing = apply_smearing)                                      #H
    lep_ZL, antilep_ZL, lep_WL, antilep_WL = build_decay_products(eventsL_dict, apply_smearing = apply_smearing)     
       
    #trasversale
    LHE_T = "unweighted_events_T.lhe"
    eventsT_dict = read_file(LHE_T)
    
    W_fm_T = get_fm_of([24, -24], eventsT_dict, apply_smearing = apply_smearing)                                 #W
    Z_fm_T = get_fm_of([23], eventsT_dict, apply_smearing = apply_smearing)                                      #Z
    H_fm_T = get_fm_of([25], eventsT_dict, apply_smearing = apply_smearing)                                      #H
    lep_ZT,  antilep_ZT, lep_WT, antilep_WT = build_decay_products(eventsT_dict, apply_smearing = apply_smearing)
    
    #calcolo angoli
    print("---- COMPUTING VARIABLES ----")
    
    delta_phi_ZL = get_deltaphi_of(Z_fm_L, lep_ZL, antilep_ZL, H_fm_L)
    delta_phi_ZT = get_deltaphi_of(Z_fm_T, lep_ZT, antilep_ZT, H_fm_T)
    delta_phi_WL = get_deltaphi_of(W_fm_L, lep_WL, antilep_WL, H_fm_L)
    delta_phi_WT = get_deltaphi_of(W_fm_T, lep_WT, antilep_WT, H_fm_T)
  
    phi1_ZL = get_phi1_of(Z_fm_L, lep_ZL, antilep_ZL, H_fm_L)
    phi1_ZT = get_phi1_of(Z_fm_T, lep_ZT,  antilep_ZT, H_fm_T)
    phi1_WL = get_phi1_of(W_fm_L, lep_WL, antilep_WL, H_fm_L)
    phi1_WT = get_phi1_of(W_fm_T, lep_WT,  antilep_WT, H_fm_T)
       
    print("---- BUILDING HISTOGRAMS ----")
    #visualizzazione distribuzioni
    data = phi1_ZL
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)
    
    data = phi1_WL
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_bins = int((max(data) - min(data)) / bin_width)
    
     
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(phi1_ZL, bins=35, color='royalblue', edgecolor = 'steelblue', stat='probability', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(phi1_ZT, bins=35, color='firebrick', stat='probability', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversale')
    ax[0].set_xlabel("|φ1| (rad)")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione angolo φ1 (Z)")
    ax[0].legend()
    
    sb.histplot(phi1_WL, bins=45, color='royalblue', edgecolor = 'steelblue', stat='probability', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(phi1_WT, bins=45, color='firebrick', stat='probability', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversale')
    ax[1].set_xlabel("|φ1| (rad)")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione angolo φ1 (W)")
    ax[1].legend()

    #plt.savefig("Confronto massa invariante.png", dpi=300)
    plt.show()
    
    
    data = delta_phi_ZL
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)
    
    data = delta_phi_WL
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_bins = int((max(data) - min(data)) / bin_width)
    
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(delta_phi_ZL, bins=35, color='royalblue', edgecolor = 'steelblue', stat='probability', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(delta_phi_ZT, bins=35, color='firebrick', stat='probability', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversale')
    ax[0].set_xlabel("Δφ (rad)")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione angolo Δφ (Z)")
    ax[0].legend()
    
    sb.histplot(delta_phi_WL, bins=45, color='royalblue', edgecolor = 'steelblue', stat='probability', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(delta_phi_WT, bins=45, color='firebrick', stat='probability', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversale')
    ax[1].set_xlabel("Δφ (rad)")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione angolo Δφ (W)")
    ax[1].legend()

    #plt.savefig("Confronto massa invariante.png", dpi=300)
    plt.show()
   
    
if __name__ == '__main__' :
    main() 
