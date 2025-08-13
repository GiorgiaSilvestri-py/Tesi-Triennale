import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils import read_file, contains_particle, build_fm_ZW, boost_to_rf, theta_star

def main () :
    '''
    Confronta le distribuzioni dell'angolo theta*, definito come la direzione di volo di Z (o W) calcolato nel sdr del centro di massa VH
    '''
    
    #calcolo del quadrimomento VH
    H_fm_L, H_fm_T = [], []
    Z_fm_L, Z_fm_T = [], []
    W_fm_L, W_fm_T = [], []    
    
    #longitudinale
    LHE_L = "unweighted_events_L.lhe"
    events_L = read_file(LHE_L)
    
    quadrivector_ZL, quadrivector_WL = build_fm_ZW(events_L)                                #liste di quadrimomenti per Z e W
    quadrivector_HL = []                                                                    #lista di quadrimomenti per l'higgs
    
    for event in events_L :
        for p in event :
            if p["pid"] == 25 :
                h_vec = vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])
                
                quadrivector_HL.append(h_vec)
                break
    '''
    #somma dei due
    fm_ZH_L = [x+y for x, y in zip(quadrivector_ZL, quadrivector_HL)]                       #quadrimomento ZH
    fm_WH_L = [x+y for x, y in zip(quadrivector_WL, quadrivector_HL)]                       #quadrimomento WH
    '''
    
    #trasversale
    LHE_T = "unweighted_events_T.lhe"
    events_T = read_file(LHE_T)

    quadrivector_ZT, quadrivector_WT = build_fm_ZW(events_T)                                #liste di quadrimomenti per Z e W
    quadrivector_HT = []                                                                    #lista di quadrimomenti per l'higgs

    for event in events_T :
        for p in event :
            if p["pid"] == 25 :
                h_vec = vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])
                
                quadrivector_HT.append(h_vec)
                break
    '''
    #somma dei due
    fm_ZH_T = [x+y for x, y in zip(quadrivector_ZT, quadrivector_HT)]                       #quadrimomento ZH
    fm_WH_T = [x+y for x, y in zip(quadrivector_WT, quadrivector_HT)]                       #quadrimomento WH
    
    
    #applico boost di lorentz
    fm_ZH_L_rf = boost_to_rf(quadrivector_ZL, quadrivector_HL)
    fm_WH_L_rf = boost_to_rf(quadrivector_WL, quadrivector_HL)
    fm_ZH_T_rf = boost_to_rf(quadrivector_ZT, quadrivector_HT)
    fm_WH_T_rf = boost_to_rf(quadrivector_WT, quadrivector_HT)
    '''
    
    #calcolo l'angolo
    theta_star_ZL = theta_star(quadrivector_ZL, quadrivector_HL)
    theta_star_ZT = theta_star(quadrivector_ZT, quadrivector_HT)
    theta_star_WL = theta_star(quadrivector_WL, quadrivector_HL)
    theta_star_WT = theta_star(quadrivector_WT, quadrivector_HT)

    #visualizzazione distribuzioni
    data = theta_star_ZL
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)
    print(nBins)
    
    data = theta_star_WL
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_bins = int((max(data) - min(data)) / bin_width)
    
     
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(theta_star_ZL, bins=50, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(theta_star_ZT, bins=50, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversale')
    ax[0].set_xlabel("cosθ*")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione angolo θ* (Z)")
    ax[0].legend()
    
    sb.histplot(theta_star_WL, bins=50, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(theta_star_WT, bins=50, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversale')
    ax[1].set_xlabel("cosθ*")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione angolo θ* (W)")
    ax[1].legend()

    plt.savefig("Confronto coseno.png", dpi=300)
    plt.show()
    
  
if __name__ == '__main__' :
    main() 
