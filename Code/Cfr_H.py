import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils import read_file, contains_particle, build_fm_ZW
import sys

            
def main():
    '''
    Confronto tramite pt_h, eta_h, e massa invariante VH
    '''
    
    #analisi polarizzazione longitudinale
    LHE_L = "unweighted_events_L.lhe"
    events_H_L = read_file(LHE_L)
    #print("N events_L: ", len(events_L))
    
    pt_H_L, eta_H_L = [], []
        
    for event in events_H_L :
        vector_H_L = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for i in event :
            vecH_L = vector.obj(px = i["px"], py = i["py"], pz = i["pz"], E = i["E"])
            if i["pid"] == 25 :
                pt_H_L.append(vecH_L.pt)
                eta_H_L.append(vecH_L.eta)
                    
     
    #analisi polarizzazione trasversale
    LHE_T = "unweighted_events_T.lhe"
    events_H_T = read_file(LHE_T)
    #print("N events_T: ", len(events_T))

    pt_H_T, eta_H_T = [], []
    
    for event in events_H_T :
        vector_H_T = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for i in event :
            vecH_T = vector.obj(px = i["px"], py = i["py"], pz = i["pz"], E = i["E"])
            if i["pid"] == 23 :
                pt_H_T.append(vecH_T.pt)
                eta_H_T.append(vecH_T.eta)
                
 
    #rappresentazione confronto
    data = pt_H_L
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)
    
    data = eta_H_L
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_bins = int((max(data) - min(data)) / bin_width)
     
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(pt_H_L, bins=nBins, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(pt_H_T, bins=nBins, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversale')
    ax[0].set_xlabel("Pt (GeV)")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione momento trasverso (H)")
    ax[0].legend()
    
    sb.histplot(eta_H_L, bins=n_bins, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(eta_H_T, bins=n_bins, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversale')
    ax[1].set_xlabel("η (rad)")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione pseudorapidità (H)")
    ax[1].legend()

    #plt.savefig("Confronto Higgs.png", dpi=300)
    plt.show()
    
    
    #wasserstein distance pt
    distance_z = wasserstein_distance(pt_H_L, pt_H_T)
    k_s_z = ks_2samp(pt_H_L, pt_H_T)

    #bin con massima distanza
    height_l, bin_l = np.histogram(pt_H_L, bins=40, density=True)
    height_t, bin_t = np.histogram(pt_H_T, bins=40, density=True)
    bin_center = [0.5 * (bin_l[i] + bin_l[i+1]) for i in range(len(bin_l) - 1)]
    difference = [abs(x-y) for x,y in zip(height_l, height_t)]
    index_max = np.argmax(difference)
    pt_max = bin_center[index_max]

    print(f"EMD (pt): {distance_z:.2f} Gev")
    print(f"Kolmogorov-Smirnov stat: {k_s_z.statistic:.3f}, p-value: {k_s_z.pvalue:.3f}")
    print(f"Valore di pt per cui si ha massima differenza: {pt_max:0f}")
    
    #wasserstein distanze eta
    distance_w = wasserstein_distance(eta_H_L, eta_H_T)
    k_s_w = ks_2samp(eta_H_L, eta_H_T)
    
    #bin con massima distanza w
    height_l, bin_l = np.histogram(eta_H_L, bins=40, density=True)
    height_t, bin_t = np.histogram(eta_H_T, bins=40, density=True)
    bin_center = [0.5 * (bin_l[i] + bin_l[i+1]) for i in range(len(bin_l) - 1)]
    difference = [abs(x-y) for x,y in zip(height_l, height_t)]
    index_max = np.argmax(difference)
    eta_max = bin_center[index_max]
    
    print(f"EMD (eta): {distance_w:.2f} rad")
    print(f"Kolmogorov-Smirnov stat: {k_s_w.statistic:.3f}, p-value: {k_s_w.pvalue:.3f}")
    print(f"Valore di eta per cui si ha massima differenza: {eta_max:0f}")
    
    
    
    #confronto massa invariante
    H_fm_L, H_fm_T = [], []
    Z_fm_L, Z_fm_T = [], []
    W_fm_L, W_fm_T = [], []    
    
    #longitudinale
    LHE_L = "unweighted_events_L.lhe"
    events_L = read_file(LHE_L)
    
    quadrivector_ZL, quadrivector_WL = build_fm_ZW(events_L)
    quadrivector_HL = []                                                                #lista di quadrimomenti per l'higgs
    
    for event in events_L :
        for p in event :
            if p["pid"] == 25 :
                h_vec = vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])
                
                quadrivector_HL.append(h_vec)
                break
    
    
    #somma dei due e massa invariante
    fm_ZH_L = [x+y for x, y in zip(quadrivector_ZL, quadrivector_HL)]
    fm_WH_L = [x+y for x, y in zip(quadrivector_WL, quadrivector_HL)]
    
    M_ZH_l = [zh.M for zh in fm_ZH_L if zh.M > 0]                                        #lista con massa invariante di ZH
    M_WH_l = [wh.M for wh in fm_WH_L if wh.M > 0]                                        #lista con massa invariante di WH
    
   
    #trasversale
    LHE_T = "unweighted_events_T.lhe"
    events_T = read_file(LHE_T)

    quadrivector_ZT, quadrivector_WT = build_fm_ZW(events_T)
    quadrivector_HT = []                                                                #lista di quadrimomenti per l'higgs

    for event in events_T :
        for p in event :
            if p["pid"] == 25 :
                h_vec = vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])
                
                quadrivector_HT.append(h_vec)
                break


    #somma dei due e massa invariante
    fm_ZH_T = [x+y for x, y in zip(quadrivector_ZT, quadrivector_HT)]
    fm_WH_T = [x+y for x, y in zip(quadrivector_WT, quadrivector_HT)]

    M_ZH_t = [zh.M for zh in fm_ZH_T]                                                   #lista con massa invariante di ZH
    M_WH_t = [wh.M for wh in fm_WH_T]                                                   #lista con massa invariante di WH
    
    
    #rappresentazione sul grafico
    data = M_ZH_l
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)
    
    data = M_WH_l
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_Bins = int((max(data) - min(data)) / bin_width)
    
    
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(M_ZH_l, bins=nBins, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(M_ZH_t, bins=nBins, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversale')
    ax[0].set_xlabel("M (GeV)")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione massa invariante ZH")
    ax[0].legend()
    
    sb.histplot(M_WH_l, bins=n_Bins, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(M_WH_t, bins=n_Bins, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversale')
    ax[1].set_xlabel("M (GeV)")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione massa invariante WH")
    ax[1].legend()

    plt.savefig("Confronto massa invariante.png", dpi=300)
    plt.show()
    
    #wasserstein distance zh
    distance_z = wasserstein_distance(M_ZH_l, M_ZH_t)
    k_s_z = ks_2samp(M_ZH_l, M_ZH_t)

    #bin con massima distanza
    height_l, bin_l = np.histogram(M_ZH_l, bins=40, density=True)
    height_t, bin_t = np.histogram(M_ZH_t, bins=40, density=True)
    bin_center = [0.5 * (bin_l[i] + bin_l[i+1]) for i in range(len(bin_l) - 1)]
    difference = [abs(x-y) for x,y in zip(height_l, height_t)]
    index_max = np.argmax(difference)
    m_max = bin_center[index_max]

    print(f"EMD (M_ZH): {distance_z:.2f} Gev")
    print(f"Kolmogorov-Smirnov stat: {k_s_z.statistic:.3f}, p-value: {k_s_z.pvalue:.3f}")
    print(f"Valore di M per cui si ha massima differenza: {m_max:0f}")
    
    #wasserstein distanze wh
    distance_w = wasserstein_distance(M_WH_l, M_WH_t)
    k_s_w = ks_2samp(M_WH_l, M_WH_t)
    
    #bin con massima distanza w
    height_l, bin_l = np.histogram(M_WH_l, bins=40, density=True)
    height_t, bin_t = np.histogram(M_WH_t, bins=40, density=True)
    bin_center = [0.5 * (bin_l[i] + bin_l[i+1]) for i in range(len(bin_l) - 1)]
    difference = [abs(x-y) for x,y in zip(height_l, height_t)]
    index_max = np.argmax(difference)
    m_max = bin_center[index_max]
    
    print(f"EMD (M_WH): {distance_w:.2f} rad")
    print(f"Kolmogorov-Smirnov stat: {k_s_w.statistic:.3f}, p-value: {k_s_w.pvalue:.3f}")
    print(f"Valore di eta per cui si ha massima differenza: {m_max:0f}")
    
if __name__ == '__main__':
    main()  
                    
