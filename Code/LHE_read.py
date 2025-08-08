import matplotlib.pyplot as plt
import math
import numpy as np
import sys

def sturges (N) :
    return int(np.ceil(1+np.log2(N)))

def pt_da_lhe(nome_file):
    pt_list_z, eta_list_z, theta_list_z, phi_list_z = [], [], [], []
    pt_list_w, eta_list_w, theta_list_w, phi_list_w = [], [], [], []
    event__z = set()
    event__w = set()
    event_counter = 0  

    with open(nome_file, 'r') as f:
        in_event = False
        event_lines = []

        for line in f:
            line = line.strip()
            if "<event>" in line:
                in_event = True
                event_lines = []
                continue
            elif "</event>" in line:
                in_event = False
                event_counter +=1  
                particles = []

                for data in event_lines:
                    parameter = data.strip().split()
                    if len(parameter) > 9:
                        pid = int(parameter[0])
                        status = int(parameter[1])
                        mother1 = int(parameter[2])
                        mother2 = int(parameter[3])
                        px = float(parameter[6])
                        py = float(parameter[7])
                        pz = float(parameter[8])
                        particles.append({
                            "pid": pid,
                            "status": status,
                            "mother1": mother1,
                            "mother2": mother2,
                            "px": px,
                            "py": py,
                            "pz": pz
                        })
                found_z = False
                found_w = False 
                
                for p in particles:
                    if p["status"] == 1 and abs(p["pid"]) in [11, 13, 15] :
                        mother_index = set([p["mother1"], p["mother2"]])
                        for mother_idx in mother_index :
                            if 1 <= mother_idx <= len(particles):
                                mother_pid = particles[mother_idx - 1]["pid"]
                                pt = math.sqrt(p["px"]**2 + p["py"]**2)
                                p_mod = math.sqrt(p["px"]**2 + p["py"]**2 + p["pz"]**2)
                                theta = np.arccos(p["pz"] / p_mod)
                                eta = -np.log(np.tan(theta / 2))
                                phi = math.atan2(p["py"], p["px"])

                                if mother_pid == 23:
                                    pt_list_z.append(pt)
                                    theta_list_z.append(theta)
                                    eta_list_z.append(eta)
                                    phi_list_z.append(phi)
                                    found_z = True
                                elif abs(mother_pid) == 24:
                                    pt_list_w.append(pt)
                                    theta_list_w.append(theta)
                                    eta_list_w.append(eta)
                                    phi_list_w.append(phi)
                                    found_w = True
                
                if found_z:
                    event__z.add(event_counter)
                if found_w:
                    event__w.add(event_counter)
                
                
                continue

            if in_event:
                event_lines.append(line)

    return pt_list_z, theta_list_z, eta_list_z, phi_list_z, pt_list_w, theta_list_w, eta_list_w, phi_list_w, event__z, event__w, event_counter
    
def main():
    
    file_lhe = "unweighted_events_50000.lhe"

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
    


