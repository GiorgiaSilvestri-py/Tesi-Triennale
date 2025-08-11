import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from utils import read_file
import sys

def contains_particle (particle_list, particle_ID) :
    for p in particle_list :
        if (p["pid"]) == particle_ID :
            return True
            
def main():
    
    #analisi polarizzazione longitudinale
    LHE_L = "unweighted_events_L.lhe"
    events_L = read_file(LHE_L)
    print(len(events_L))
    
    events_Z_L = [e for e in events_L if contains_particle (e, 23)]
    events_W_L = [e for e in events_L if contains_particle (e, 24) or contains_particle (e, -24)]
    print(len(events_Z_L), len(events_W_L))
    
    pt_Z_L = []
    pt_W_L = []    
    
    for event in events_Z_L :
        vector_Z_L = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for i in event :
            vecZ_L = vector.obj(px = i["px"], py = i["py"], pz = i["pz"], E = i["E"])
            if i["pid"] == 23 :
                pt_Z_L.append(vecZ_L.pt)
                
    for event in events_W_L :
        vector_W_L = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for i in event :
            vecW_L = vector.obj(px = i["px"], py = i["py"], pz = i["pz"], E = i["E"])
            if abs(i["pid"]) == 24 :
                pt_W_L.append(vecW_L.pt)
      
    #analisi polarizzazione trasversale
    LHE_T = "unweighted_events_T.lhe"
    events_T = read_file(LHE_T)
    #print(len(events_T))
    
    events_Z_T = [e for e in events_T if contains_particle (e, 23)]
    events_W_T = [e for e in events_T if contains_particle (e, 24) or contains_particle (e, -24)]
    #print(len(events_Z_T), len(events_W_T))
    
    pt_Z_T = []
    pt_W_T = []    
    
    for event in events_Z_T :
        vector_Z_T = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for i in event :
            vecZ_T = vector.obj(px = i["px"], py = i["py"], pz = i["pz"], E = i["E"])
            if i["pid"] == 23 :
                pt_Z_T.append(vecZ_T.pt)
                
    for event in events_W_T :
        vector_W_T = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for i in event :
            vecW_T = vector.obj(px = i["px"], py = i["py"], pz = i["pz"], E = i["E"])
            if abs(i["pid"]) == 24 :
                pt_W_T.append(vecW_T.pt)  
    
    #numero totale di eventi per Z
    N_events = len(events_Z_T) + len(events_Z_L)
    #print(N_events) 
    
    
    data = pt_Z_L
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_binsZ = int((max(data) - min(data)) / bin_width)
      
    data = pt_W_L
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_binsW = int((max(data) - min(data)) / bin_width)
    
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(pt_Z_L, bins=n_binsZ, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(pt_Z_T, bins=n_binsZ, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversa')
    ax[0].set_xlabel("Momento trasverso (GeV)")
    ax[0].set_ylabel("Densità")
    ax[0].set_title("Distribuzione momento Z")
    ax[0].legend()
    
    sb.histplot(pt_W_L, bins=n_binsW, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(pt_W_T, bins=n_binsW, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversa')
    ax[1].set_xlabel("Momento trasverso (GeV)")
    ax[1].set_ylabel("Densità")
    ax[1].set_title("Distribuzione momento W")
    ax[1].legend()

    #plt.savefig("Confronto momento trasverso.png", dpi=300)
    plt.show()
    
    
    #confronto polarizzazioni W
    height_L, b = np.histogram(pt_W_L, 30)         #altezza delle colonne della distribuzione longitudinale
    height_T, b_t = np.histogram(pt_W_T, 30)       #altezza delle colonne della distribuzione trasversale
    #print(len(height_L), len(height_T))
    
    height_tot = np.array([x+y for x, y in zip(height_L, height_T)])
    #print((height_tot))
    
    l_fraction = np.divide(height_L, height_tot, out=np.zeros_like(height_tot, dtype=float), where=height_tot > 0)  #y1
    t_fraction = np.divide(height_T, height_tot, out=np.zeros_like(height_tot, dtype=float), where=height_tot > 0)  #y2
    #print(len(l_fraction), len(t_fraction))
    
    bin_center_W = [0.5 * (b[i] + b[i+1]) for i in range(len(b) - 1)]
    
    #confronto polarizzazioni Z
    height_L_Z, b_Z = np.histogram(pt_Z_L, 30)         #altezza delle colonne della distribuzione longitudinale
    height_T_Z, b_t_Z = np.histogram(pt_Z_T, 30)       #altezza delle colonne della distribuzione trasversale
    
    
    height_tot_Z = np.array([x+y for x, y in zip(height_L_Z, height_T_Z)])
    
    
    l_fraction_Z = np.divide(height_L_Z, height_tot_Z, out=np.zeros_like(height_tot_Z, dtype=float), where=height_tot_Z > 0)  #y1
    t_fraction_Z = np.divide(height_T_Z, height_tot_Z, out=np.zeros_like(height_tot_Z, dtype=float), where=height_tot_Z > 0)  #y2
    
    bin_center_Z = [0.5 * (b_Z[i] + b_Z[i+1]) for i in range(len(b_Z) - 1)]
    
    
    sb.set(style='whitegrid')
    plt.figure(figsize=(8,5))
    
    sb.scatterplot(x=bin_center_W, y=l_fraction, label='Longitudinale', color='royalblue')
    sb.scatterplot(x=bin_center_W, y=t_fraction, label='Trasversale', color='firebrick')
    
    plt.title('Eventi di polarizzazione (W)')
    plt.xlabel('Momento trasverso (Gev)')
    plt.ylabel('dN/N')
    plt.legend()
    plt.tight_layout()
    #plt.savefig("Confronto polarizzazioni.png", dpi=300)
    plt.show()
    
    
    sb.set(style='whitegrid')
    plt.figure(figsize=(8,5))
    
    sb.scatterplot(x=bin_center_Z, y=l_fraction_Z, label='Longitudinale', color='royalblue')
    sb.scatterplot(x=bin_center_Z, y=t_fraction_Z, label='Trasversale', color='firebrick')
    
    plt.title('Eventi di polarizzazione (Z)')
    plt.xlabel('Momento trasverso (Gev)')
    plt.ylabel('dN/N')
    plt.legend()
    plt.tight_layout()
    plt.show()
    
if __name__ == '__main__':
    main()  
                    
