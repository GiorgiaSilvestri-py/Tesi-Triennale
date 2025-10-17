import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
import random
from scipy.stats import wasserstein_distance, ks_2samp
from utils import read_file, f_event_type, is_nonzero, get_thetastar_of

def solve_equation(a, b, c):
    
    delta = b**2 - 4*a*c
    
    if delta < 0:
        delta = 0
        print("Negative delta, setting delta = 0")
        
    s1 = ( -b + np.sqrt(delta) ) / (2*a)
    s2 = ( -b - np.sqrt(delta) ) / (2*a)
    
    return s1, s2

def apply_smearing_to(fm_vector, is_neutrino = False):
    '''
    Apply smearing to a four-momentum vector with the given resolution
    '''
    #neutrinos
    if is_neutrino:
        resolution = 0.20
        rand_smearing = random.gauss(mu = 1.0, sigma = resolution)
        px_new = fm_vector.px * rand_smearing
        py_new = fm_vector.py * rand_smearing
        pz_new = fm_vector.pz * rand_smearing
        E_new = fm_vector.E * rand_smearing
        smeared_vector = vector.obj( px = px_new, py = py_new, pz = py_new, E = E_new )
    
    else:
        #leptons
        resolution = 0.01
        rand_smearing = random.gauss(mu = 1.0, sigma = resolution)
        smeared_vector = vector.obj( px = (fm_vector.px * rand_smearing), py = (fm_vector.py * rand_smearing), pz = (fm_vector.pz * rand_smearing), E =  (fm_vector.E * rand_smearing) ) 
    
    
    return smeared_vector    
    
def compute_neutrino_momentum(events_dict, e_id):
    '''
    Ricostruisce il pz del neutrino a partire dalla massa invariante della W
    '''
    '''
    Ricostruisce il pz del neutrino a partire dalla massa invariante della W
    '''
   
    vec_lepton            = vector.obj(px = 0, py = 0, pz = 0, E = 0)
    vec_neutrino_expected = vector.obj(px = 0, py = 0, pz = 0, E = 0)
    
    p_dict = events_dict[e_id]
    for pid, cinematic in p_dict.items():
        if abs(pid) in [11, 13, 15]:
            vec_lepton = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
            
        if abs(pid) in [12, 14, 16]:
            vec_neutrino_expected = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)

    #skip event if either one is None
    if not is_nonzero(vec_lepton) or not is_nonzero(vec_neutrino_expected):
        return
        
    
    return compute_neutrino_momentum_from_particles(vec_lepton, vec_neutrino_expected)
    
    
    
def compute_neutrino_momentum_from_particles(vec_lepton, vec_neutrino_expected):
    '''
    Ricostruisce il pz del neutrino a partire dalla massa invariante della W -> angolo maggiore
    '''
   
    expected_pz_neutrino = vec_neutrino_expected.pz
    
    #solvin second grade equation
    M = 80
    pt_lepton_dot_pt_neutrino = vec_lepton.px * vec_neutrino_expected.px + vec_lepton.py * vec_neutrino_expected.py
    pt_neutrino = vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2
    
    alpha = M**2 + 2 * pt_lepton_dot_pt_neutrino
    a = 4 * (vec_lepton.pz ** 2) - 4 * vec_lepton.E ** 2
    b = 4 * alpha * vec_lepton.pz
    c = alpha**2 - 4 * vec_lepton.E**2 * pt_neutrino
    
    #second grade equation solutions
    pz_neutrino_delta_plus, pz_neutrino_delta_minus = solve_equation(a, b, c)
    
    
    #computing angles between lepton and neutrino
    neutrino_energy_deltaplus = np.sqrt( vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2 + pz_neutrino_delta_plus**2 )
    built_neutrino_vector_plus = vector.obj(
                                            px = vec_neutrino_expected.px,
                                            py = vec_neutrino_expected.py, 
                                            pz = pz_neutrino_delta_plus,
                                            E  = neutrino_energy_deltaplus
                                            )
                                            
 
    neutrino_direction_plus = built_neutrino_vector_plus.to_3D()
    lepton_direction        = vec_lepton.to_3D()
    
    #first solution
    theta_plus = compute_angle_between(lepton_direction, neutrino_direction_plus)

    neutrino_energy_deltaminus = np.sqrt( vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2 + pz_neutrino_delta_minus**2 )
    built_neutrino_vector_minus = vector.obj(
                                            px = vec_neutrino_expected.px,
                                            py = vec_neutrino_expected.py, 
                                            pz = pz_neutrino_delta_minus,
                                            E  = neutrino_energy_deltaminus
                                            )
                                            
    neutrino_direction_minus = built_neutrino_vector_minus.to_3D()
    
    #second solution
    theta_minus = compute_angle_between(lepton_direction, neutrino_direction_minus)
    
    #selecting correct solution
    
    if abs(theta_plus) < abs(theta_minus):
        built_neutrino_vector = built_neutrino_vector_plus
        #pz_neutrino_built = built_neutrino_vector_plus.pz
        
    else:
        built_neutrino_vector = built_neutrino_vector_minus
        #pz_neutrino_built = built_neutrino_vector_minus.pz
   
    return built_neutrino_vector
    

def compute_angle_between(v1, v2):
    """
    Calcola il coseno dell'angolo tra due vettori tridimensionali
    """
    
    p1 = np.array([v1.x, v1.y, v1.z])
    p2 = np.array([v2.x, v2.y, v2.z])

    cos_theta = np.dot(p1, p2) / (np.linalg.norm(p1) * np.linalg.norm(p2))
    
    return cos_theta
    
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    

    
def get_fm_of(p_id, events_dict, apply_smearing = False):

    fm_dict = {}
    
    if apply_smearing:                
        if set(p_id) == set([24, -24]):
            for event_id, particles_dict in events_dict.items():
                event_type = f_event_type(particles_dict)
                
                #24 = W; 23 = Z
                if event_type == 23:
                    continue
                    
                fm_vec_1 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                fm_vec_2 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                
                for pid, cinematic in particles_dict.items():
                    if abs(pid) in [11, 13, 15]:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum)
                        fm_vec_1 = smeared_momentum
                        
                for pid, cinematic in particles_dict.items():
                    if abs(pid) in [12, 14, 16]:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, is_neutrino = True)
                        fm_vec_2 = compute_neutrino_momentum_from_particles(fm_vec_1, smeared_momentum)
                        
                if not is_nonzero(fm_vec_1) or not is_nonzero(fm_vec_2):
                    continue                       # salta l’evento se manca uno dei due
                fm_dict[event_id] = fm_vec_1 + fm_vec_2
        
        if set(p_id) == set([25]):
            for event_id, particles_dict in events_dict.items():
                for pid, cinematic in particles_dict.items():
                    if pid in p_id:
                        fm_dict[event_id] = cinematic.build_fm()
        
    
    
def main():


    #analisi polarizzazione longitudinale
    LHE_L = "unweighted_events_L.lhe"
    eventsL_dict = read_file(LHE_L)
    
    vec_exp = {}
    vec_built = {}
    w_exp = {}
    pz_exp = []
    pz_built = []
    R_l = []
    
    for evt_id, p_dict in eventsL_dict.items():
        e_type = f_event_type(p_dict)
        
        if e_type == 23:
            continue
        for pid, cinematic in p_dict.items():
            if abs(pid) in [12, 14, 16]:
                vec_exp[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                vec_built[evt_id] = compute_neutrino_momentum(eventsL_dict, evt_id)
                
            if abs(pid) in [24]:
                w_exp[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
    
    #print(vec_built)
    computed_list, exp_list_t, w_expected_l = [], [], []
    
    for e_id in vec_exp.keys():
        fm_exp = vec_exp[e_id]
        fm_built = vec_built[e_id]
        w_expected_l.append(w_exp[e_id].pz)
        
        pz_exp.append(fm_exp.pz)
        pz_built.append(fm_built.pz)
        R_l.append(fm_built.pz/fm_exp.pz)
        
        computed_list.append(fm_built.pz)
        exp_list_t.append(fm_exp.pz)
    
    
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(computed_list, bins = 80, binrange= (-500, 500), color = 'steelblue', edgecolor = 'steelblue', stat = 'probability', alpha = 0.6, label = 'comp')
    sb.histplot(exp_list_t, bins = 80,  binrange= (-500, 500), color = 'firebrick', edgecolor = 'steelblue', stat = 'probability', alpha = 0.4, label = 'exp')
    plt.xlabel("pz (GeV)")
    plt.ylabel("Event fraction")
    plt.title("Cfr computed expected")
    plt.legend()
    plt.show()
    
    
    LHE_T = "unweighted_events_T.lhe"
    eventsT_dict = read_file(LHE_T)
    
    vec_exp = {}
    vec_built = {}
    pz_exp = []
    pz_built = []
    R_t = []
    
    for evt_id, p_dict in eventsT_dict.items():
        e_type = f_event_type(p_dict)
        
        if e_type == 23:
            continue
        for pid, cinematic in p_dict.items():
            if abs(pid) in [12, 14, 16]:
                vec_exp[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                vec_built[evt_id] = compute_neutrino_momentum(eventsT_dict, evt_id)
    
    #print(vec_built)
    
    computed_list, exp_list_l = [], []
    
    for e_id in vec_exp.keys():
        fm_exp = vec_exp[e_id]
        fm_built = vec_built[e_id]
        
        pz_exp.append(fm_exp.pz)
        pz_built.append(fm_built.pz)
        R_t.append(fm_built.pz/fm_exp.pz)
            
        computed_list.append(fm_built.pz)
        exp_list_l.append(fm_exp.pz)
    
    print("min: ", min(R_l), min(R_t))
    print("max: ", max(R_l), max(R_t))
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(computed_list, bins = 80, binrange= (-500, 500), color = 'steelblue', edgecolor = 'steelblue', stat = 'probability', alpha = 0.6, label = 'computed')
    sb.histplot(exp_list_l, bins = 80,  binrange= (-500, 500), color = 'firebrick', edgecolor = 'steelblue', stat = 'probability', alpha = 0.4, label = 'expected')
    plt.xlabel("pz (GeV)")
    plt.ylabel("Event fraction")
    plt.title("Cfr computed expected")
    plt.legend()
    plt.show()
    
    sb.set(style='whitegrid')
    sb.histplot(R_l, bins = 80, binrange= (0, 7), color = 'steelblue', edgecolor = 'steelblue', stat = 'probability', alpha = 0.6, label = 'long')
    sb.histplot(R_t, bins = 80,  binrange= (0, 7), color = 'firebrick', edgecolor = 'steelblue', stat = 'probability', alpha = 0.4, label = 'trans')
    plt.xlabel("Computed pz / expected pz")
    plt.ylabel("Event fraction")
    plt.legend()
    plt.show()
    
    #compute RMS
    

    
    #costruzione theta star
    fm_w_dict = get_fm_of([24, -24], eventsL_dict, apply_smearing = True)     
    #build_decay_products(eventsL_dict, apply_smearing = True)
    fm_h_all_dict = get_fm_of([25], eventsL_dict, apply_smearing = True)                               
    
    
    #filtro H
    fm_h_dict = {}   
    fm_w_l = []
    
    for event_id in fm_w_dict.keys() & fm_h_all_dict.keys():
        fm_w     = fm_w_dict[event_id]
        fm_higgs = fm_h_all_dict[event_id]        
        
        fm_h_dict[event_id] = fm_higgs
        fm_w_l.append(fm_w.pz)
        
    cos_l = get_thetastar_of(fm_w_dict, fm_h_dict)
   
    
    
    fm_w_dict = get_fm_of([24, -24], eventsT_dict, apply_smearing = True)     
    fm_h_all_dict = get_fm_of([25], eventsT_dict, apply_smearing = True)                               
    
    
    #filtro H
    fm_h_dict = {}           
    fm_w_t = []
    
    for event_id in fm_w_dict.keys() & fm_h_all_dict.keys():
        fm_w     = fm_w_dict[event_id]
        fm_higgs = fm_h_all_dict[event_id]        
        
        fm_h_dict[event_id] = fm_higgs
        fm_w_t.append(fm_w.pz)
        
    cos_t = get_thetastar_of(fm_w_dict, fm_h_dict)
    

    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(cos_l, bins = 70, binrange= (-1, 1), color = 'steelblue', edgecolor = 'steelblue', stat = 'probability', alpha = 0.6, label = 'longitudinal polarisation')
    sb.histplot(cos_t, bins = 70,  binrange= (-1, 1), color = 'firebrick', edgecolor = 'steelblue', stat = 'probability', alpha = 0.4, label = 'transverse polarisation')
    plt.xlabel("cos(θ*)")
    plt.ylabel("Event fraction")         
    plt.legend()
    plt.show()
    
    '''
    #analisi polarizzazione longitudinale
    LHE_L = "unweighted_events_L.lhe"
    eventsL_dict = read_file(LHE_L)
    
    vec_exp = {}
    vec_built = {}
    pz_exp = []
    pz_built = []
    R_l = []
    
    for evt_id, p_dict in eventsL_dict.items():
        e_type = f_event_type(p_dict)
        
        if e_type == 23:
            continue
        for pid, cinematic in p_dict.items():
            if abs(pid) in [12, 14, 16]:
                vec_exp[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                vec_built[evt_id] = compute_neutrino_momentum(eventsL_dict, evt_id)
    
    #print(vec_built)
    
    for e_id in vec_exp.keys():
        fm_exp = vec_exp[e_id]
        fm_built = vec_built[e_id]
        
        pz_exp.append(fm_exp.pz)
        pz_built.append(fm_built.pz)
        R_l.append(fm_built.pz/fm_exp.pz)
        
    LHE_T = "unweighted_events_T.lhe"
    eventsT_dict = read_file(LHE_T)
    
    vec_exp = {}
    vec_built = {}
    pz_exp = []
    pz_built = []
    R_t = []
    
    for evt_id, p_dict in eventsT_dict.items():
        e_type = f_event_type(p_dict)
        
        if e_type == 23:
            continue
        for pid, cinematic in p_dict.items():
            if abs(pid) in [12, 14, 16]:
                vec_exp[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                vec_built[evt_id] = compute_neutrino_momentum(eventsT_dict, evt_id)
    
    #print(vec_built)
    
    computed_list, exp_list = [], []
    
    for e_id in vec_exp.keys():
        fm_exp = vec_exp[e_id]
        fm_built = vec_built[e_id]
        
        pz_exp.append(fm_exp.pz)
        pz_built.append(fm_built.pz)
        R_t.append(fm_built.pz/fm_exp.pz)
            
        computed_list.append(fm_built.pz)
        exp_list.append(fm_exp.pz)
    
    print("min: ", min(R_l), min(R_t))
    print("max: ", max(R_l), max(R_t))
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(computed_list, bins = 80, binrange= (-500, 500), color = 'steelblue', edgecolor = 'steelblue', stat = 'probability', alpha = 0.6, label = 'longitudinal polarisation')
    sb.histplot(exp_list, bins = 80,  binrange= (-500, 500), color = 'firebrick', edgecolor = 'steelblue', stat = 'probability', alpha = 0.4, label = 'transverse polarisation')
    plt.xlabel("Computed pz / expected pz")
    plt.ylabel("Event fraction")         
    plt.legend()
    plt.show()
    
    #compute RMS
    '''
    '''
    N = len(R_l)
    r_ridotto = [x for x in R_l if x < 2]
    print(r_ridotto)
    R_sq = [x*x for x in R_l]
    somma = sum(R_sq)
    RMS = np.sqrt(np.mean(r_ridotto))
    print(RMS)
    '''
    '''
    r_ridotto = [x for x in R_l if x < 2 and x > 0]
    R_sq = [x**2 for x in r_ridotto]
    RMS = np.sqrt(np.mean(R_sq))
    print("longitudinale:" ,RMS)
    
    r_ridotto = [x for x in R_t if x < 2 and x > 0]
    R_sq = [x**2 for x in r_ridotto]
    RMS = np.sqrt(np.mean(R_sq))
    print("transv:", RMS)

    '''
    
    
    '''
    pz_1_n_list = []
    pz_2_n_list = []
    vec2_pz_list = []

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
        
        a_z = vec_1.pz / vec_1.E
        d = 1 / vec_1.E * ( M**2/2 +  vec_1_2d.dot(vec_2_2d))
        a=1
        b = 2*d*a_z / (a_z**2-1)
        c = (d**2 - (vec_2_2d.px**2 + vec_2_2d.py**2) )/ (a_z**2-1)
        
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
        
        d = -b / (2*a)
    
        if pz_delta_plus > d:
            pz_elected[e_id] = pz_delta_plus
        if pz_delta_minus > d:
            pz_elected[e_id] = pz_delta_minus
        if pz_delta_minus > d and pz_delta_plus > d:
            #print("Both values are acceptable")
            pz_elected[e_id] = random.choice([pz_delta_plus, pz_delta_minus])
        
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
        theta_p = compute_angle_between(vec_1_p, vec_2_p)
        theta_p_list.append(theta_p)
        
        neutrino_E_2 = np.sqrt(vec_2.px**2 + vec_2.py**2 + pz_2_n**2)
        new_vec_neutrino_minus = vector.obj(px = vec_2.px, py = vec_2.py, pz = pz_2_n, E = neutrino_E_2) 
        vec_2_m = new_vec_neutrino_minus.to_3D()
        vec_1_m = vec_1.to_3D()
        theta_m = compute_angle_between(vec_1_m, vec_2_m)
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
    '''   
    
    #print(min(R), max(R))
    
   

    

  
        
if __name__ == '__main__':
    main()  
                    
