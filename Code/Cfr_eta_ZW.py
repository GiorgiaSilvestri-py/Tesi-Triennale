import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils import read_file
import sys

def contains_particle (particle_list, particle_ID) :
    for p in particle_list :
        if (p["pid"]) == particle_ID :
            return True
            
def main():
    '''
    confronto tramite pseudorapidità di Z e W e (solo grafico) tramite angolo phi
    '''
    
    #analisi polarizzazione longitudinale
    LHE_L = "unweighted_events_L.lhe"
    events_L = read_file(LHE_L)
    
    events_Z_L = [e for e in events_L if contains_particle (e, 23)]
    events_W_L = [e for e in events_L if contains_particle (e, 24) or contains_particle (e, -24)]
    #print(len(events_Z_L), len(events_W_L))
    
    eta_Z_L, phi_Z_L = [], []
    eta_W_L, phi_W_L = [], []
    
    for event in events_Z_L :
        vector_Z_L = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for i in event :
            vecZ_L = vector.obj(px = i["px"], py = i["py"], pz = i["pz"], E = i["E"])
            if i["pid"] == 23 :
                eta_Z_L.append(vecZ_L.eta)
                phi_Z_L.append(vecZ_L.phi)
                
    for event in events_W_L :
        vector_W_L = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for i in event :
            vecW_L = vector.obj(px = i["px"], py = i["py"], pz = i["pz"], E = i["E"])
            if abs(i["pid"]) == 24 :
                eta_W_L.append(vecW_L.eta)
                phi_W_L.append(vecW_L.phi)
      
    #analisi polarizzazione trasversale
    LHE_T = "unweighted_events_T.lhe"
    events_T = read_file(LHE_T)
    #print(len(events_T))
    
    events_Z_T = [e for e in events_T if contains_particle (e, 23)]
    events_W_T = [e for e in events_T if contains_particle (e, 24) or contains_particle (e, -24)]
    #print(len(events_Z_T), len(events_W_T))
    
    eta_Z_T, phi_Z_T = [], []
    eta_W_T, phi_W_T = [], []    
    
    for event in events_Z_T :
        vector_Z_T = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for i in event :
            vecZ_T = vector.obj(px = i["px"], py = i["py"], pz = i["pz"], E = i["E"])
            if i["pid"] == 23 :
                eta_Z_T.append(vecZ_T.eta)
                phi_Z_T.append(vecZ_T.phi)
                
    for event in events_W_T :
        vector_W_T = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for i in event :
            vecW_T = vector.obj(px = i["px"], py = i["py"], pz = i["pz"], E = i["E"])
            if abs(i["pid"]) == 24 :
                eta_W_T.append(vecW_T.eta)
                phi_W_T.append(vecW_T.phi)
    
    #numero totale di eventi per Z
    N_events = len(events_Z_T) + len(events_Z_L)
    #print(N_events) 
    
    #rappresentazione distribuzioni eta
    data = eta_Z_L
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_binsZ = int((max(data) - min(data)) / bin_width)
      
    data = eta_W_L
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_binsW = int((max(data) - min(data)) / bin_width)
    
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(19, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(eta_Z_L, bins=n_binsZ, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(eta_Z_T, bins=n_binsZ, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversa')
    ax[0].set_xlabel("Pseudorapidità (rad)")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione pseudorapidità Z")
    ax[0].legend()
    
    sb.histplot(eta_W_L, bins=n_binsW, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(eta_W_T, bins=n_binsW, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversa')
    ax[1].set_xlabel("Pseudorapidità (rad)")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione pseudorapidità W")
    ax[1].legend()

    #plt.savefig("Confronto pseudorapidità.png", dpi=300)
    plt.show()
    
    
    
    #rappresentazione distribuzioni phi
    data = phi_Z_L
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_binsZ = int((max(data) - min(data)) / bin_width)
      
    data = phi_W_L
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_binsW = int((max(data) - min(data)) / bin_width)
    
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(19, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(phi_Z_L, bins=n_binsZ, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(phi_Z_T, bins=n_binsZ, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversa')
    ax[0].set_xlabel("Pseudorapidità (rad)")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione angolo azimutale Z")
    ax[0].legend()
    
    sb.histplot(phi_W_L, bins=n_binsW, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(phi_W_T, bins=n_binsW, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversa')
    ax[1].set_xlabel("Pseudorapidità (rad)")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione angolo azimutale W")
    ax[1].legend()

    #plt.savefig("Confronto phi.png", dpi=300)
    plt.show()
    
    
    
    
    #wasserstein distance
    distance_z = wasserstein_distance(eta_Z_L, eta_Z_T)
    k_s_z = ks_2samp(eta_Z_L, eta_Z_T)

    #bin con massima distanza
    height_l, bin_l = np.histogram(eta_Z_L, bins=40, density=True)
    height_t, bin_t = np.histogram(eta_Z_T, bins=40, density=True)
    bin_center_Z = [0.5 * (bin_l[i] + bin_l[i+1]) for i in range(len(bin_l) - 1)]
    difference = [abs(x-y) for x,y in zip(height_l, height_t)]
    index_max = np.argmax(difference)
    eta_max = bin_center_Z[index_max]

    print(f"EMD (Z): {distance_z:.2f} rad")
    print(f"Kolmogorov-Smirnov stat: {k_s_z.statistic:.3f}, p-value: {k_s_z.pvalue:.3f}")
    print(f"Valore di eta per cui si ha massima differenza: {eta_max:0f}")
    
    
    distance_w = wasserstein_distance(eta_W_L, eta_W_T)
    k_s_w = ks_2samp(eta_W_L, eta_W_T)
    
    #bin con massima distanza w
    height_l, bin_l = np.histogram(eta_W_L, bins=40, density=True)
    height_t, bin_t = np.histogram(eta_W_T, bins=40, density=True)
    bin_center_W = [0.5 * (bin_l[i] + bin_l[i+1]) for i in range(len(bin_l) - 1)]
    difference = [abs(x-y) for x,y in zip(height_l, height_t)]
    index_max = np.argmax(difference)
    eta_max = bin_center_W[index_max]
    
    print(f"EMD (W): {distance_w:.2f} rad")
    print(f"Kolmogorov-Smirnov stat: {k_s_w.statistic:.3f}, p-value: {k_s_w.pvalue:.3f}")
    print(f"Valore di eta per cui si ha massima differenza: {eta_max:0f}")
 
if __name__ == '__main__':
    main() 
