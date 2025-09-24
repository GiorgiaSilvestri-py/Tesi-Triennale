import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils import *
import sys

def contains_particle (particle_list, particle_ID) :
    for p in particle_list :
        if (p["pid"]) == particle_ID :
            return True
            
def main():
    '''
    Confronto tramite pt Z e W + grafico scatter
    '''
    apply_smearing = True
    print("---- BUILDING DICTIONARIES ----")
    #analisi polarizzazione longitudinale
    LHE_L = "unweighted_events_L.lhe"
    events_l = read_file(LHE_L)
    
    #analisi polarizzazione trasversale
    LHE_T = "unweighted_events_T.lhe"
    events_t = read_file(LHE_T)

    pt_z_l = get_pt_of([23], events_l, apply_smearing = apply_smearing)
    pt_w_l = get_pt_of([24, -24], events_l, apply_smearing = apply_smearing)
    pt_z_t = get_pt_of([23], events_t, apply_smearing = apply_smearing)
    pt_w_t = get_pt_of([24, -24], events_t, apply_smearing = apply_smearing)
    
    eta_z_l = get_eta_of([23], events_l, apply_smearing = apply_smearing)
    #print(len(eta_z_l))
    # print(len(ta_z_l))
    eta_w_l = get_eta_of([24, -24], events_l, apply_smearing = apply_smearing)
    eta_z_t = get_eta_of([23], events_t, apply_smearing = apply_smearing)
    eta_w_t = get_eta_of([24, -24], events_t, apply_smearing = apply_smearing)
    
    print("---- BUILDING HISTOGRAMS ----")
    #rappresentazione
    data = list(pt_z_l.values())
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_binsZ = int((max(data) - min(data)) / bin_width)
      
    data = list(pt_w_l.values())
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_binsW = int((max(data) - min(data)) / bin_width)
    
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(pt_z_l, bins=n_binsZ, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(pt_z_t, bins=n_binsZ, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversale')
    ax[0].set_xlabel("Pt (GeV)")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione momento trasverso Z")
    ax[0].legend()
    
    sb.histplot(pt_w_l, bins=n_binsW, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(pt_w_t, bins=n_binsW, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversale')
    ax[1].set_xlabel("Pt (GeV)")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione momento trasverso W")
    ax[1].legend()

    plt.savefig("Confronto pt V.png", dpi=300)
    plt.show()
    
    print("---- BUILDING HISTOGRAMS ----")
    #rappresentazione distribuzioni eta
    data = list(eta_z_l.values())
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_binsZ = int((max(data) - min(data)) / bin_width)
      
    data = list(eta_w_l.values())
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_binsW = int((max(data) - min(data)) / bin_width)
    
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(19, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(eta_z_l, bins=n_binsZ, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(eta_z_t, bins=n_binsZ, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversale')
    ax[0].set_xlabel("η (rad)")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione pseudorapidità Z")
    ax[0].legend()
    
    sb.histplot(eta_w_l, bins=n_binsW, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(eta_w_t, bins=n_binsW, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversale')
    ax[1].set_xlabel("η (rad)")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione pseudorapidità W")
    ax[1].legend()

    plt.savefig("Confronto eta V.png", dpi=300)
    plt.show()
    
    '''
    #confronto polarizzazioni W
    height_L, b = np.histogram(pt_w_l, 30)         #altezza delle colonne della distribuzione longitudinale
    height_T, b_t = np.histogram(pt_w_t, 30)       #altezza delle colonne della distribuzione trasversale
    #print(len(height_L), len(height_T))
    
    height_tot = np.array([x+y for x, y in zip(height_L, height_T)])
    #print((height_tot))
    
    l_fraction = np.divide(height_L, height_tot, out=np.zeros_like(height_tot, dtype=float), where=height_tot > 0)  #y1
    t_fraction = np.divide(height_T, height_tot, out=np.zeros_like(height_tot, dtype=float), where=height_tot > 0)  #y2
    #print(len(l_fraction), len(t_fraction))
    
    bin_center_W = [0.5 * (b[i] + b[i+1]) for i in range(len(b) - 1)]
    
    #confronto polarizzazioni Z
    height_L_Z, b_Z = np.histogram(pt_z_l, 30)         #altezza delle colonne della distribuzione longitudinale
    height_T_Z, b_t_Z = np.histogram(pt_z_t, 30)       #altezza delle colonne della distribuzione trasversale
    
    
    height_tot_Z = np.array([x+y for x, y in zip(height_L_Z, height_T_Z)])
    
    
    l_fraction_Z = np.divide(height_L_Z, height_tot_Z, out=np.zeros_like(height_tot_Z, dtype=float), where=height_tot_Z > 0)  #y1
    t_fraction_Z = np.divide(height_T_Z, height_tot_Z, out=np.zeros_like(height_tot_Z, dtype=float), where=height_tot_Z > 0)  #y2
    
    bin_center_Z = [0.5 * (b_Z[i] + b_Z[i+1]) for i in range(len(b_Z) - 1)]
    
    
    sb.set(style='whitegrid')
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)

    #polarizzazione W
    sb.scatterplot(x=bin_center_W, y=l_fraction, label='Longitudinale', color='royalblue', ax=axes[0])
    sb.scatterplot(x=bin_center_W, y=t_fraction, label='Trasversale', color='firebrick', ax=axes[0])
    axes[0].set_title('Eventi di polarizzazione (W)')
    axes[0].set_xlabel('Momento trasverso (Gev)')
    axes[0].set_ylabel('dN/N')
    axes[0].legend()

    #polarizzazione Z
    sb.scatterplot(x=bin_center_Z, y=l_fraction_Z, label='Longitudinale', color='royalblue', ax=axes[1])
    sb.scatterplot(x=bin_center_Z, y=t_fraction_Z, label='Trasversale', color='firebrick', ax=axes[1])
    axes[1].set_title('Eventi di polarizzazione (Z)')
    axes[1].set_xlabel('Momento trasverso (Gev)')
    axes[1].set_ylabel('dN/N')
    axes[1].legend()
    
    plt.tight_layout()
    plt.show()
    '''
    
    print("---- COMPUTING METRICS ----")
    #wasserstein distance
    distance_z = wasserstein_distance(list(pt_z_l.values()), list(pt_z_t.values()))
    k_s_z = ks_2samp(list(pt_z_l.values()), list(pt_z_t.values()))

    #bin con massima distanza
    height_l, bin_l = np.histogram(list(pt_z_l.values()), bins=40, density=True)
    height_t, bin_t = np.histogram(list(pt_z_t.values()), bins=40, density=True)
    bin_center_Z = [0.5 * (bin_l[i] + bin_l[i+1]) for i in range(len(bin_l) - 1)]
    difference = [abs(x-y) for x,y in zip(height_l, height_t)]
    index_max = np.argmax(difference)
    pt_max = bin_center_Z[index_max]

    print(f"EMD (Z): {distance_z:.2f} Gev")
    print(f"Kolmogorov-Smirnov stat: {k_s_z.statistic:.3f}, p-value: {k_s_z.pvalue:.3f}")
    print(f"Valore di pt per cui si ha massima differenza: {pt_max:0f}")
    
    
    distance_w = wasserstein_distance(list(pt_w_l.values()), list(pt_w_t.values()))
    k_s_w = ks_2samp(list(pt_w_l.values()), list(pt_w_t.values()))
    
    #bin con massima distanza w
    height_l, bin_l = np.histogram(list(pt_w_l.values()), bins=40, density=True)
    height_t, bin_t = np.histogram(list(pt_w_t.values()), bins=40, density=True)
    bin_center_W = [0.5 * (bin_l[i] + bin_l[i+1]) for i in range(len(bin_l) - 1)]
    difference = [abs(x-y) for x,y in zip(height_l, height_t)]
    index_max = np.argmax(difference)
    pt_max = bin_center_W[index_max]
    
    print(f"EMD (W): {distance_w:.2f} Gev")
    print(f"Kolmogorov-Smirnov stat: {k_s_w.statistic:.3f}, p-value: {k_s_w.pvalue:.3f}")
    print(f"Valore di pt per cui si ha massima differenza: {pt_max:0f}")
  
    
if __name__ == '__main__':
    main()  
                    
