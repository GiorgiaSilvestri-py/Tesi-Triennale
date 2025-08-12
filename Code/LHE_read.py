import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from utils import read_file
import sys

def sturges (N) :
    return int(np.ceil(1+np.log2(N)))
    
def contains_particle (particle_list, particle_ID) :
    for p in particle_list :
        if (p["pid"]) == particle_ID :
            return True
    
def main():
    
    file_lhe = "unweighted_events_50000.lhe"
    events = read_file(file_lhe)
    #print(len(events))
    
    events_Z = [e for e in events if contains_particle (e, 23)]
    events_W = [e for e in events if contains_particle (e, 24) or contains_particle (e, -24)]
    #print(len(events_Z), len(events_W))
    M_z = []
    pt_lep, eta_lep, phi_lep = [], [], []
    
    for event in events_Z :
        #print(len(event))
        v_Z = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for p in event :
            #print(p)
            my_vec = vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])
            if p["status"] == 1 and abs(p["pid"]) in [11, 13, 15] :
                pt_lep.append(my_vec.pt)        #momento pt leptone
                eta_lep.append(my_vec.eta)      #eta leptone
                phi_lep.append(my_vec.phi)      #phi leptone
                v_Z += my_vec                   #quadrimomento somma
                
        if v_Z.M > 0 : 
            M_z.append(v_Z.M)                   #massa invariante Z
        
    data = pt_lep  
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)
            
    fig, ax = plt.subplots(2, 2, figsize = (16, 5))
    
    sb.histplot(pt_lep, bins = nBins, color = 'firebrick', edgecolor = 'firebrick', element='step', stat = 'density', ax = ax[0,0], linewidth=1.5, alpha = 0.4, label = 'Polarizzazione trasversale')
    ax[0, 0].set_xlabel("P_t (GeV)")
    ax[0, 0].set_ylabel("dN/N")
    ax[0, 0].set_title("Distribuzione momento trasverso leptone")
    ax[0, 0].legend()

    sb.histplot(eta_lep, bins = nBins, color = 'firebrick', edgecolor = 'firebrick', element='step',stat = 'density', ax = ax[0,1], linewidth=1.5, alpha = 0.4, label = 'Polarizzazione trasversale')
    ax[0, 1].set_xlabel("η (rad)")
    ax[0, 1].set_ylabel("dN/N")
    ax[0, 1].set_title("Distribuzione pseudorapidità leptone")
    ax[0, 1].legend()
    
    sb.histplot(phi_lep, bins = nBins, color = 'firebrick', edgecolor = 'firebrick', element='step',stat = 'density', ax = ax[1,0], linewidth=1.5, alpha = 0.4, label = 'Polarizzazione trasversale')
    ax[1, 0].set_xlabel("φ (rad)")
    ax[1, 0].set_ylabel("dN/N")
    ax[1, 0].set_title("Distribuzione angolo azimutale leptone")
    ax[1, 0].legend()
    
    sb.histplot(M_z, bins = nBins, color = 'firebrick', edgecolor = 'firebrick', element='step',stat = 'density', ax = ax[1,1], linewidth=1.5, alpha = 0.4, label = 'Polarizzazione trasversale')
    ax[1, 1].set_xlabel("M (Gev)")
    ax[1, 1].set_ylabel("dN/N")
    ax[1, 1].set_title("Massa invariante")
    ax[1, 1].legend()
    
    plt.tight_layout()
    plt.savefig("Grafici_Z.png", dpi=300)
    plt.show()
    
    
    #uguale procedimento per W
    M_W = []
    pt_lep, eta_lep, phi_lep = [], [], []
    
    for event in events_W :
        v_W = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for p in event :
            #print(p)
            my_vec = vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])
            if p["status"] == 1 and abs(p["pid"]) in [11, 13, 15] :
                pt_lep.append(my_vec.pt)        #momento pt leptone
                eta_lep.append(my_vec.eta)      #eta leptone
                phi_lep.append(my_vec.phi)      #phi leptone
                v_W += my_vec                   #quadrimomento somma
                
        if v_W.M > 0 : 
            M_W.append(v_W.M)                   #massa invariante W
        
    data = pt_lep  
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)
            
    fig, ax = plt.subplots(2, 2, figsize = (16, 5))
    
    sb.histplot(pt_lep, bins = nBins, color = 'steelblue', edgecolor = 'steelblue', element='step', stat = 'density', ax = ax[0,0], linewidth=1.5, alpha = 0.6, label = 'Polarizzazione longitudinale')
    ax[0, 0].set_xlabel("P_t (GeV)")
    ax[0, 0].set_ylabel("dN/N")
    ax[0, 0].set_title("Distribuzione momento trasverso leptone")
    ax[0, 0].legend()

    sb.histplot(eta_lep, bins = nBins, color = 'steelblue', edgecolor = 'steelblue', element='step',stat = 'density', ax = ax[0,1], linewidth=1.5, alpha = 0.6, label = 'Polarizzazione longitudinale')
    ax[0, 1].set_xlabel("η (rad)")
    ax[0, 1].set_ylabel("dN/N")
    ax[0, 1].set_title("Distribuzione pseudorapidità leptone")
    ax[0, 1].legend()
    
    sb.histplot(phi_lep, bins = nBins, color = 'steelblue', edgecolor = 'steelblue', element='step',stat = 'density', ax = ax[1,0], linewidth=1.5, alpha = 0.6, label = 'Polarizzazione longitudinale')
    ax[1, 0].set_xlabel("φ (rad)")
    ax[1, 0].set_ylabel("dN/N")
    ax[1, 0].set_title("Distribuzione angolo azimutale leptone")
    ax[1, 0].legend()
    
    sb.histplot(M_z, bins = nBins, color = 'steelblue', edgecolor = 'steelblue', element='step',stat = 'density', ax = ax[1,1], linewidth=1.5, alpha = 0.6, label = 'Polarizzazione longitudinale')
    ax[1, 1].set_xlabel("M (Gev)")
    ax[1, 1].set_ylabel("dN/N")
    ax[1, 1].set_title("Massa invariante")
    ax[1, 1].legend()
    
    plt.tight_layout()
    plt.savefig("Grafici_W.png", dpi=300)
    plt.show()
    
if __name__ == '__main__':
    main()
    


