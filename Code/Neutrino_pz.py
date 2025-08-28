import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils import *

def solve_equation(a, b, c):
    
    delta = b**2 - 4*a*c
    
    #if delta < b :
        #print(f"{delta} < {b}")
    
    if delta < 0:
        delta = 0
        #print("Invalid sqrt value")
        
        
       
    s1 = ( -b + np.sqrt(delta) ) / (2*a)
    s2 = ( -b - np.sqrt(delta) ) / (2*a)
    
    return s1, s2
    
def big_solve_equation(a, b, c):
    
    delta = b**2 - 4*a*c
    
    if delta < b :
        print(f"{delta} < {b}")
    
    if delta < 0:
        delta = 0
        print("Invalid sqrt value")
    
    return -b/(2*a), delta/(4*a**2)    
    
    
def main():

    print("---- BUILDING DICTIONARIES ----")
    #analisi polarizzazione longitudinale
    LHE_L = "unweighted_events_L.lhe"
    eventsL_dict = read_file(LHE_L)
    
    counter_plus = 0
    counter_minus = 0
    counter_evt = 0
    counter_neg = 0
    
    pz_1_n_list = []
    pz_2_n_list = []
    vec2_pz_list = []
    abs_diff_plus_list = []
    abs_diff_minus_list = []
    fm_w_built = {}
    theta_p_list, theta_m_list = [], []
    
    pz_neutrino_exp = {}
    pz_neutrino_built = {}
    pz_elected = {}
    
    
    fm_exp = {}
    
    
    for e_id, p_dict in eventsL_dict.items():
        
        vec_1 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        vec_2 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        
        for pid, cinematic in p_dict.items():
            if abs(pid) in [11, 13, 15]:
                vec_1 = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
            if abs(pid) in [12, 14, 16]:
                vec_2 = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
    
        if not is_nonzero(vec_1) or not is_nonzero(vec_2):
            continue                                        
        fm_exp[e_id] = vec_1 + vec_2
        pz_neutrino_exp[e_id] = vec_2.pz
        counter_evt += 1
        
        vec_1_2d = vector.obj(px = vec_1.px, py = vec_1.py, pz = 0)
        
        vec_2_2d = vector.obj(px = vec_2.px, py = vec_2.py, pz = 0)
        
        M = 80
        #scritto a mano con massa leptone nulla
        '''
        a_z = vec_1.pz / vec_1.E
        d = 1 / vec_1.E * ( M**2/2 +  vec_1_2d.dot(vec_2_2d))
        a=1
        b = 2*d*a_z / (a_z**2-1)
        c = (d**2 - (vec_2_2d.px**2 + vec_2_2d.py**2) )/ (a_z**2-1)
        '''
        #soluzione a mano (me)
        #alpha = M**2 + 2 * vec_1_2d.dot(vec_2_2d)
        alpha = M**2 + 2 * (vec_1.px*vec_2.px + vec_1.py*vec_2.py)
        a = 4 * (vec_1.pz ** 2) - 4 * vec_1.E ** 2
        b = 4 * alpha * vec_1.pz
        c = alpha**2 - 4 * vec_1.E**2 * (vec_2.px**2 + vec_2.py**2) 
        
        
        
        #scritto a mano
        #A = M**2 + 2*(vec_1_2d.dot(vec_2_2d))
        
        
        #a = 4*(vec_1.E**2 - vec_1.pz**2)
        #b = -4 * A * vec_1.pz
        #c = 4*vec_1.E**2 * (vec_2_2d.px**2 + vec_2_2d.py**2) - A**2
        
        #TEORIco
        #a = (vec_1.pz)**2 - (vec_1.E)**2
        #b = M**2 * (vec_1.pz) + 2 * (vec_1.pz) * (vec_1_2d.dot(vec_2_2d))
        #c = M**4 / 4 + (vec_1_2d.dot(vec_2_2d))**2 + M**2 * (vec_1_2d.dot(vec_2_2d)) - (vec_1.E)**2 * (vec_2_2d.px**2 + vec_2_2d.py**2) 
        
        #risoluzione eq di secondo grado
        pz_delta_plus, pz_delta_minus = solve_equation(a, b, c)
        #first, second = big_solve_equation(a, b, c)
        '''
        d = -b / (2*a)
    
        if pz_delta_plus > d:
            pz_elected[e_id] = pz_delta_plus
        if pz_delta_minus > d:
            pz_elected[e_id] = pz_delta_minus
        if pz_delta_minus > d and pz_delta_plus > d:
            #print("Both values are acceptable")
            pz_elected[e_id] = random.choice([pz_delta_plus, pz_delta_minus])
        '''  
        
        if pz_delta_plus is None and pz_delta_minus is None:
            counter_neg += 1
            continue
        
        #calcolo pz neutrino per i due casi
        pz_1_n = pz_delta_plus 
        pz_2_n = pz_delta_minus

        
        neutrino_E_1 = np.sqrt(vec_2.px**2 + vec_2.py**2 + pz_1_n**2)
        new_vec_neutrino_plus = vector.obj(px = vec_2.px, py = vec_2.py, pz = pz_1_n, E = neutrino_E_1) 
        vec_2_p = new_vec_neutrino_plus.to_3D()
        vec_1_p = vec_1.to_3D()
        theta_p = compute_angle(vec_1_p, vec_2_p)
        theta_p_list.append(theta_p)
        
        neutrino_E_2 = np.sqrt(vec_2.px**2 + vec_2.py**2 + pz_2_n**2)
        new_vec_neutrino_minus = vector.obj(px = vec_2.px, py = vec_2.py, pz = pz_2_n, E = neutrino_E_2) 
        vec_2_m = new_vec_neutrino_minus.to_3D()
        vec_1_m = vec_1.to_3D()
        theta_m = compute_angle(vec_1_m, vec_2_m)
        theta_m_list.append(theta_m)
        
        if abs(theta_p) < abs(theta_m):
            vec_2_rec = new_vec_neutrino_plus
            pz_neutrino_built[e_id] = vec_2_rec.pz
            
        else :
            vec_2_rec = new_vec_neutrino_minus
            pz_neutrino_built[e_id] = vec_2_rec.pz
            
        
            
        fm_w_built[e_id] = vec_1 + vec_2_rec
        #print(f"Evento: {e_id}; {fm_w} \n")
    
        
        
        
    R = [] #confronto energie w
    R_pz = [] #confronto pz neutrino con angoli
    R_pz_2 = [] #confronto pz con minimo parabola
    
    for event_id in fm_w_built.keys():
        w_built = fm_w_built[event_id]
        w_exp = fm_exp[event_id]
        #pz_2 = pz_elected[event_id] 
        pz_neutrino = pz_neutrino_built[event_id]
        pz_exp = pz_neutrino_exp[event_id]      
        
        R.append(w_built.E / w_exp.E)
        R_pz.append(pz_neutrino / pz_exp)
        #R_pz_2.append(pz_2 - pz_exp)
        
    
    print(min(R), max(R))
    
    plt.figure()
    sb.set(style='whitegrid')
    #sb.histplot(R, bins=60, binrange = (0, 2), color = 'firebrick', edgecolor = 'firebrick', stat = 'probability', alpha = 0.4, label = 'E W')
    sb.histplot(R_pz, bins = 60, binrange = (0, 2), color = 'steelblue', edgecolor = 'steelblue', stat = 'probability', alpha = 0.6, label = 'Pz neutrino')
    #sb.histplot(R_pz_2, bins = 50, binrange = (np.min(R), 10), color = 'green', edgecolor = 'brown', stat = 'density', alpha = 0.3, label = 'Pz neutrino - 2')
    #plt.xticks(np.linspace(0, 10, 11))
    #plt.yticks(np.linspace(0, 1, 11)) 
    plt.xlabel("Computed / Expected")
    plt.ylabel("Event fraction")        
    plt.title("WL polarization")        
    plt.show()

    

  
        
if __name__ == '__main__':
    main()  
                    
