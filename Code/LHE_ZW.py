import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sns
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
    
    sns.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sns.histplot(pt_Z_L, bins=n_binsZ, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sns.histplot(pt_Z_T, bins=n_binsZ, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversa')
    ax[0].set_xlabel("Momento trasverso (GeV)")
    ax[0].set_ylabel("Densità")
    ax[0].set_title("Distribuzione momento Z")
    ax[0].legend()
    
    sns.histplot(pt_W_L, bins=n_binsW, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sns.histplot(pt_W_T, bins=n_binsW, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversa')
    ax[1].set_xlabel("Momento trasverso (GeV)")
    ax[1].set_ylabel("Densità")
    ax[1].set_title("Distribuzione momento W")
    ax[1].legend()

    plt.savefig("Confronto momento trasverso.png", dpi=300)
    plt.show()
    
    '''
    ax[0].hist(pt_Z_L, bins=n_binsZ, color = 'blue', density = True,  alpha = 0.5)
    ax[0].hist(pt_Z_T, bins=n_binsZ, color = 'red', edgecolor = 'red', histtype='stepfilled', density = True, hatch='/', alpha = 0.3)
    ax[0].set_xlabel("Momento trasverso (Gev)")
    ax[0].set_ylabel("Conteggi")
    ax[0].set_title("Distribuzione momento Z (LT)")
    
    ax[1].hist(pt_W_L, bins=n_binsW, color = 'blue', density = True, alpha = 0.5)
    ax[1].hist(pt_W_T, bins=n_binsW, color = 'red', edgecolor = 'red', histtype='stepfilled', hatch='/', density = True, alpha = 0.3)
    ax[1].set_xlabel("Momento trasverso (Gev)")
    ax[1].set_ylabel("Conteggi")
    ax[1].set_title("Distribuzione momento W (LT)")
    
    #plt.savefig("Pt polarizzazione longitudinale.png", dpi = 300)
    plt.show()

    '''
  
if __name__ == '__main__':
    main()  
                    
