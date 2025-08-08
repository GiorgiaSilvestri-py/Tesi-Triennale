import matplotlib.pyplot as plt
import math
import numpy as np
import vector
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
    print(len(events))
    
    events_Z = [e for e in events if contains_particle (e, 23)]
    events_W = [e for e in events if contains_particle (e, 24) or contains_particle (e, -24)]
    #print(len(events_Z), len(events_W))
    M_z = []
    
    for event in events_Z :
        #print(len(event))
        v_Z = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        for p in event :
            #print(p)
            my_vec = vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])
            if p["status"] == 1 and abs(p["pid"]) in [11, 13, 15] :
                v_Z += my_vec
                
        if v_Z.M > 0 : 
            M_z.append(v_Z.M)
        
    plt.hist(M_z, bins = 300)
    plt.xlabel ("M (Gev)")
    plt.ylabel ("Conteggi")
    plt.title ("Distribuzione momento trasverso da Z")
    plt.show()
        
        
    
    
    return 0
        

    #id_number = int(input("ID particella: "))
    
    #momento trasverso
    pt_z, theta_list_z, eta_list_z, phi_list_z, pt_w, theta_list_w, eta_list_w, phi_list_w, z_decay, w_decay, n_events = pt_da_lhe(file_lhe)
    
    print("Numero totale di eventi: ", n_events), 
    print("Numero produzioni ZH: ", len(z_decay))
    print("Numero produzioni WH: ", len(w_decay))
    
    #pseudorapidità e angolo polare
    #print("Angolo polare: ", theta_list)
    print("Max Pseudorapidity Z: ", max(eta_list_z))
    print("Max Pseudorapidity W: ", max(eta_list_w))
    
    #angolo phi
    #print("Angolo azimutale phi: ", phi_list)
    
    data = pt_z  
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_bins = int((max(data) - min(data)) / bin_width)
    
    '''
    n_bins = int(input("Numero bin: "))
    
    
    '''
    plt.hist(pt_z, bins = n_bins)
    plt.xlabel ("Momento trasverso (Gev)")
    plt.ylabel ("Conteggi")
    plt.title ("Distribuzione momento trasverso da Z")
    plt.show()
    
    plt.hist(pt_w, bins = n_bins)
    plt.xlabel ("Momento trasverso (Gev)")
    plt.ylabel ("Conteggi")
    plt.title ("Distribuzione momento trasverso da W")
    plt.show()
    
    '''
    sturges(len(eta_list_z))
    xMin = min(eta_list_z)
    xMax = max(eta_list_z)
    bin_edges = np.linspace (xMin, xMax, n_bins)
    '''
    
    data = eta_list_z  
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_bins = int((max(data) - min(data)) / bin_width)
    
    plt.hist(eta_list_z, bins = n_bins)
    plt.xlabel ("Psuedorapidity (rad)")
    plt.ylabel ("Conteggi")
    plt.title ("Distribuzione pseudorapidità (Z)")
    plt.show()
    
    '''
    sturges(len(eta_list_w))
    xMin = min(eta_list_w)
    xMax = max(eta_list_w)
    bin_edges = np.linspace (xMin, xMax, n_bins)
    '''
    
    data = eta_list_w  
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_bins = int((max(data) - min(data)) / bin_width)
    
    plt.hist(eta_list_w, bins = n_bins)
    plt.xlabel ("Psuedorapidity (rad)")
    plt.ylabel ("Conteggi")
    plt.title ("Distribuzione pseudorapidità (W)")
    plt.show()
    '''
    sturges(len(phi_list_z))
    xMin = min(phi_list_z)
    xMax = max(phi_list_z)
    bin_edges = np.linspace (xMin, xMax, n_bins)
    '''
    
    n = len(phi_list_z)
    #nbins = int(np.ceil(2 * n ** (1/3)))
    std_dev = np.std(phi_list_z)
    bin_width = 3.5 * std_dev / n ** (1/3)

    # Calcolo del numero di bin
    phi_list_z_range = np.max(phi_list_z) - np.min(phi_list_z)
    nbins = int(np.ceil(phi_list_z_range / bin_width))
    
    plt.hist(phi_list_z, bins = 'auto')
    plt.xlabel ("Angolo azimutale (rad)")
    plt.ylabel ("Conteggi")
    plt.title ("Distribuzione angolo azimutale (Z)")
    plt.show()
    '''
    sturges(len(phi_list_w))
    xMin = min(phi_list_w)
    xMax = max(phi_list_w)
    bin_edges = np.linspace (xMin, xMax, n_bins)
    '''
    
    n = len(phi_list_w)
    #nbins = int(np.ceil(2 * n ** (1/3)))
    std_dev = np.std(phi_list_w)
    bin_width = 3.5 * std_dev / n ** (1/3)

    # Calcolo del numero di bin
    phi_list_w_range = np.max(phi_list_w) - np.min(phi_list_w)
    nbins = int(np.ceil(phi_list_w_range / bin_width))
    
    plt.hist(phi_list_w, bins = nbins)
    plt.xlabel ("Angolo azimutale (rad)")
    plt.ylabel ("Conteggi")
    plt.title ("Distribuzione angolo azimutale (W)")
    plt.show()
    
if __name__ == '__main__':
    main()
    


