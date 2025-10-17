import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import random
import seaborn as sb
from utils_nonlhe import solve_equation, read_file, f_event_type, get_fm_of, get_thetastar_of


def compute_neutrino_momentum(events_dict, e_id):
   
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
    
def compute_neutrino_momentum_diff(events_dict, e_id):
   
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
        
  
    return compute_neutrino_momentum_from_particles_diff(vec_lepton, vec_neutrino_expected)
    

def compute_neutrino_momentum_from_particles(vec_lepton, vec_neutrino_expected):
    
    #compute pz with minimum angle 3d
    
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
    cos_theta_plus = compute_angle_between(lepton_direction, neutrino_direction_plus)
    
    neutrino_energy_deltaminus = np.sqrt( vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2 + pz_neutrino_delta_minus**2 )
    built_neutrino_vector_minus = vector.obj(
                                            px = vec_neutrino_expected.px,
                                            py = vec_neutrino_expected.py, 
                                            pz = pz_neutrino_delta_minus,
                                            E  = neutrino_energy_deltaminus
                                            )
            
    neutrino_direction_minus = built_neutrino_vector_minus.to_3D()
    #second solution
    cos_theta_minus = compute_angle_between(lepton_direction, neutrino_direction_minus)         #coseno
    
    
    #selecting correct solution
    if (cos_theta_plus) > (cos_theta_minus):
        built_neutrino_vector = built_neutrino_vector_plus
        #pz_neutrino_built = built_neutrino_vector_plus.pz
        cos = cos_theta_plus
        
    else:
        built_neutrino_vector = built_neutrino_vector_minus
        #pz_neutrino_built = built_neutrino_vector_minus.pz
        cos = cos_theta_minus
   
        
    return built_neutrino_vector, cos


def compute_neutrino_momentum_from_particles_diff(vec_lepton, vec_neutrino_expected):
   
    #compute pz cosidering lepton pz   
   
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
    
    neutrino_energy_deltaplus = np.sqrt( vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2 + pz_neutrino_delta_plus**2 )
    built_neutrino_vector_plus = vector.obj(
                                            px = vec_neutrino_expected.px,
                                            py = vec_neutrino_expected.py, 
                                            pz = pz_neutrino_delta_plus,
                                            E  = neutrino_energy_deltaplus
                                            )
  
    neutrino_direction_plus = built_neutrino_vector_plus.to_3D()
    lepton_direction        = vec_lepton.to_3D()
    
    
  
    neutrino_energy_deltaminus = np.sqrt( vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2 + pz_neutrino_delta_minus**2 )
    built_neutrino_vector_minus = vector.obj(
                                            px = vec_neutrino_expected.px,
                                            py = vec_neutrino_expected.py, 
                                            pz = pz_neutrino_delta_minus,
                                            E  = neutrino_energy_deltaminus
                                            )
            
    
    #computing differences with lepton pz
    diff_plus = abs(pz_neutrino_delta_plus - vec_lepton.pz)
    diff_minus = abs(pz_neutrino_delta_minus - vec_lepton.pz)
   
    #selecting correct solution

    if diff_plus < diff_minus:
        built_neutrino_vector = built_neutrino_vector_plus
        
    else:
        built_neutrino_vector = built_neutrino_vector_minus
        
    return built_neutrino_vector

def compute_angle_between(v1, v2):
    """
    Calcola il coseno dell'angolo tra due vettori tridimensionali
    """
    
    p1 = np.array([v1.x, v1.y, v1.z])
    p2 = np.array([v2.x, v2.y, v2.z])

    cos_theta = np.dot(p1, p2) / (np.linalg.norm(p1) * np.linalg.norm(p2))
    
    return cos_theta
    
def main():

   
    
    LHE_L = "unweighted_events_L.lhe"
    eventsL_dict = read_file(LHE_L)
    

    #I method - minimum angle
    vec_lepton = {}
    vec_neutrino = {}
    w_exp = {}
    fm_w_dict = {}
    vec_built = {}
    cos_list = []
    
    for evt_id, p_dict in eventsL_dict.items():
        e_type = f_event_type(p_dict)
        
        if e_type == 23:
            continue
        for pid, cinematic in p_dict.items():
            if abs(pid) in [12, 14, 16]:
                vec_neutrino[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
        for pid, cinematic in p_dict.items():     
            if abs(pid) in [11, 13, 15]:
                vec_lepton[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
                vec_built[evt_id], cos = compute_neutrino_momentum_from_particles(vec_lepton[evt_id], vec_neutrino[evt_id])
                cos_list.append(cos)
                
            if abs(pid) in [24]:
                w_exp[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
        fm_w_dict[evt_id] = vec_lepton[evt_id] + vec_built[evt_id]
        
    fm_h_all_dict = get_fm_of([25], eventsL_dict, apply_smearing = False)                               
    
    
    fm_h_dict = {} 
    fm_pz_list_l = []
    
    for event_id in fm_w_dict.keys() & fm_h_all_dict.keys() & vec_lepton.keys():
        fm_lep = vec_lepton[event_id]
        fm_w     = fm_w_dict[event_id]
        fm_higgs = fm_h_all_dict[event_id]        
        
        fm_h_dict[event_id] = fm_higgs
        fm_pz_list_l.append(abs(fm_lep.pz))
        
        
    cos_l_angle = get_thetastar_of(fm_w_dict, fm_h_dict)
    
    decay_angle_list, decay_angle_list_built = [], []
    decay_cos_list, decay_cos_list_built = [], []
    pz_diff, pz_diff_built_angle_l = [], []
    
    R_l_angle = []
    
    for e_id in w_exp.keys():
        fm_w = w_exp[e_id]
        fm_n = vec_neutrino[e_id]
        fm_n_built = vec_built[e_id]
        fm_l = vec_lepton[e_id]
        
        neutrino_direction = fm_n.to_3D()
        n_built_direction = fm_n_built.to_3D()
        lepton_direction   = fm_l.to_3D()
        
        cos_theta = compute_angle_between(neutrino_direction, lepton_direction)
        cos_theta_b = compute_angle_between(n_built_direction, lepton_direction)
        
        decay_cos_list.append(cos_theta)
        decay_cos_list_built.append(cos_theta_b)
        
        decay_angle_list.append( np.arccos(cos_theta) )
        decay_angle_list_built.append( np.arccos(cos_theta_b) )
        
        pz_diff.append(abs(fm_n.pz - fm_l.pz))
        pz_diff_built_angle_l.append(abs(fm_n_built.pz - fm_l.pz))
        
        R_l_angle.append(fm_n_built.pz / fm_n.pz)
        
    
    #------------------------------------------------------------------------------------------------------------
    
    #CFR histograms - minimum angle
    
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(decay_cos_list, bins = 80, binrange= (-1, 1), color = 'steelblue',  element = 'step', label = 'expected', edgecolor = 'steelblue', stat = 'probability', alpha = 0.6)
    sb.histplot(decay_cos_list_built, bins = 80, binrange= (-1, 1), color = 'firebrick', label = 'computed', element = 'step', edgecolor = 'steelblue', stat = 'probability', alpha = 0.4)
    plt.xlabel("coseno")
    plt.ylabel("Event fraction")  
    plt.title("Cos distribution - minimum angle")
    plt.legend()
    plt.show()
    
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(decay_angle_list, bins = 80, color = 'steelblue', element = 'step', edgecolor = 'steelblue', label = 'expected', stat = 'probability', alpha = 0.6)
    sb.histplot(decay_angle_list_built, bins = 80, color = 'firebrick', element = 'step', label = 'computed',  edgecolor = 'steelblue', stat = 'probability', alpha = 0.4)
    plt.xlabel("Angle (rad)")
    plt.ylabel("Event fraction")
    plt.title("Angle distribution - minimum angle")
    plt.legend()
    plt.show()
    
    
    
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(pz_diff, bins = 80, color = 'steelblue',  element = 'step', edgecolor = 'steelblue', label = 'expected', stat = 'probability', alpha = 0.6)
    sb.histplot(pz_diff_built, bins = 80, color = 'firebrick', element = 'step', label = 'computed', edgecolor = 'steelblue', stat = 'probability', alpha = 0.4)
    plt.xlabel("Pz Difference")
    plt.ylabel("Event fraction")  
    plt.title("Pz lepton-neutrino distribution - minimum angle")   
    plt.legend()
    plt.show()


    #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    
    #ii method - PZ difference
    
    vec_lepton = {}
    vec_neutrino = {}
    w_exp = {}
    vec_built = {}
    cos_list = []
    
    for evt_id, p_dict in eventsL_dict.items():
        e_type = f_event_type(p_dict)
        
        if e_type == 23:
            continue
        for pid, cinematic in p_dict.items():
            if abs(pid) in [12, 14, 16]:
                vec_neutrino[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
        for pid, cinematic in p_dict.items():     
            if abs(pid) in [11, 13, 15]:
                vec_lepton[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
                vec_built[evt_id] = compute_neutrino_momentum_from_particles_diff(vec_lepton[evt_id], vec_neutrino[evt_id])
                
            if abs(pid) in [24]:
                w_exp[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
        fm_w_dict[evt_id] = vec_lepton[evt_id] + vec_built[evt_id]
                
                
            
    fm_h_all_dict = get_fm_of([25], eventsL_dict, apply_smearing = False)                               
    
    
    fm_h_dict = {}   
    
    for event_id in fm_w_dict.keys() & fm_h_all_dict.keys():
        fm_w     = fm_w_dict[event_id]
        fm_higgs = fm_h_all_dict[event_id]        
        
        fm_h_dict[event_id] = fm_higgs
        
        
    cos_l_diff = get_thetastar_of(fm_w_dict, fm_h_dict)
    
    decay_angle_list, decay_angle_list_built = [], []
    decay_cos_list, decay_cos_list_built = [], []
    pz_diff, pz_diff_built = [], []
    
    R_l_diff = []
    
    for e_id in w_exp.keys():
        fm_w = w_exp[e_id]
        fm_n = vec_neutrino[e_id]
        fm_n_built = vec_built[e_id]
        fm_l = vec_lepton[e_id]
        
        neutrino_direction = fm_n.to_3D()
        n_built_direction = fm_n_built.to_3D()
        lepton_direction   = fm_l.to_3D()
        
        cos_theta = compute_angle_between(neutrino_direction, lepton_direction)
        cos_theta_b = compute_angle_between(n_built_direction, lepton_direction)
        
        decay_cos_list.append(cos_theta)
        decay_cos_list_built.append(cos_theta_b)
        
        decay_angle_list.append( np.arccos(cos_theta) )
        decay_angle_list_built.append( np.arccos(cos_theta_b) )
        
        pz_diff.append(abs(fm_n.pz - fm_l.pz))
        pz_diff_built.append(abs(fm_n_built.pz - fm_l.pz))
        
        R_l_diff.append(fm_n_built.pz / fm_n.pz)
        
    
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    

    #CFR histograms - pz difference
    
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(decay_cos_list, bins = 80, binrange= (-1, 1), color = 'steelblue',  element = 'step', label = 'expected', edgecolor = 'steelblue', stat = 'probability', alpha = 0.6)
    sb.histplot(decay_cos_list_built, bins = 80, binrange= (-1, 1), color = 'firebrick', label = 'computed', element = 'step', edgecolor = 'steelblue', stat = 'probability', alpha = 0.4)
    plt.xlabel("coseno")
    plt.ylabel("Event fraction")  
    plt.title("Cos distribution - pz difference")
    plt.legend()
    plt.show()
    
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(decay_angle_list, bins = 80, color = 'steelblue', element = 'step', edgecolor = 'steelblue', label = 'expected', stat = 'probability', alpha = 0.6)
    sb.histplot(decay_angle_list_built, bins = 80, color = 'firebrick', element = 'step', label = 'computed',  edgecolor = 'steelblue', stat = 'probability', alpha = 0.4)
    plt.xlabel("Angle (rad)")
    plt.ylabel("Event fraction")
    plt.title("Angle distribution - pz difference")
    plt.legend()
    plt.show()
    
    
    
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(pz_diff, bins = 80, color = 'steelblue',  element = 'step', edgecolor = 'steelblue', label = 'expected', stat = 'probability', alpha = 0.6)
    sb.histplot(pz_diff_built, bins = 80, color = 'firebrick', element = 'step', label = 'computed', edgecolor = 'steelblue', stat = 'probability', alpha = 0.4)
    plt.xlabel("Pz Difference")
    plt.ylabel("Event fraction")  
    plt.title("Pz lepton-neutrino distribution - pz difference")   
    plt.legend()
    plt.show()
    
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    LHE_T = "unweighted_events_T.lhe"
    eventsT_dict = read_file(LHE_T)
    
    vec_lepton = {}
    vec_neutrino = {}
    w_exp = {}
    vec_built = {}
    cos_list = []
    
    #I method - minimum angle
    for evt_id, p_dict in eventsT_dict.items():
        e_type = f_event_type(p_dict)
        
        if e_type == 23:
            continue
        for pid, cinematic in p_dict.items():
            if abs(pid) in [12, 14, 16]:
                vec_neutrino[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
        for pid, cinematic in p_dict.items():     
            if abs(pid) in [11, 13, 15]:
                vec_lepton[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
                vec_built[evt_id], cos = compute_neutrino_momentum_from_particles(vec_lepton[evt_id], vec_neutrino[evt_id])
                cos_list.append(cos)
                
            if abs(pid) in [24]:
                    w_exp[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
        fm_w_dict[evt_id] = vec_lepton[evt_id] + vec_built[evt_id]
     
    fm_h_all_dict = get_fm_of([25], eventsL_dict, apply_smearing = False)                               
    
    
    fm_h_dict = {}   
    
    for event_id in fm_w_dict.keys() & fm_h_all_dict.keys():
        fm_w     = fm_w_dict[event_id]
        fm_higgs = fm_h_all_dict[event_id]        
        
        fm_h_dict[event_id] = fm_higgs
        
        
    cos_t_angle = get_thetastar_of(fm_w_dict, fm_h_dict)
            
    decay_angle_list, decay_angle_list_built = [], []
    decay_cos_list, decay_cos_list_built = [], []
    pz_diff, pz_diff_built = [], []
    
    R_t_angle = []
    
    for e_id in w_exp.keys():
        fm_w = w_exp[e_id]
        fm_n = vec_neutrino[e_id]
        fm_n_built = vec_built[e_id]
        fm_l = vec_lepton[e_id]
        
        neutrino_direction = fm_n.to_3D()
        n_built_direction = fm_n_built.to_3D()
        lepton_direction   = fm_l.to_3D()
        
        cos_theta = compute_angle_between(neutrino_direction, lepton_direction)
        cos_theta_b = compute_angle_between(n_built_direction, lepton_direction)
        
        decay_cos_list.append(cos_theta)
        decay_cos_list_built.append(cos_theta_b)
        
        decay_angle_list.append( np.arccos(cos_theta) )
        decay_angle_list_built.append( np.arccos(cos_theta_b) )
        
        pz_diff.append(abs(fm_n.pz - fm_l.pz))
        pz_diff_built.append(abs(fm_n_built.pz - fm_l.pz))
        
        R_t_angle.append(fm_n_built.pz / fm_n.pz)
        
        
        
    #second method - pz difference
    vec_lepton = {}
    vec_neutrino = {}
    w_exp = {}
    vec_built = {}
    cos_list = []
    
    for evt_id, p_dict in eventsT_dict.items():
        e_type = f_event_type(p_dict)
        
        if e_type == 23:
            continue
        for pid, cinematic in p_dict.items():
            if abs(pid) in [12, 14, 16]:
                vec_neutrino[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
        for pid, cinematic in p_dict.items():     
            if abs(pid) in [11, 13, 15]:
                vec_lepton[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
                vec_built[evt_id] = compute_neutrino_momentum_from_particles_diff(vec_lepton[evt_id], vec_neutrino[evt_id])
                
            if abs(pid) in [24]:
                w_exp[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
            
        fm_w_dict[evt_id] = vec_lepton[evt_id] + vec_built[evt_id] 
    
    fm_h_all_dict = get_fm_of([25], eventsL_dict, apply_smearing = False)                               
    
    
    fm_h_dict = {}   
    
    for event_id in fm_w_dict.keys() & fm_h_all_dict.keys():
        fm_w     = fm_w_dict[event_id]
        fm_higgs = fm_h_all_dict[event_id]        
        
        fm_h_dict[event_id] = fm_higgs
        
        
    cos_t_diff = get_thetastar_of(fm_w_dict, fm_h_dict)
   
           
    decay_angle_list, decay_angle_list_built = [], []
    decay_cos_list, decay_cos_list_built = [], []
    pz_diff, pz_diff_built = [], []
    
    R_t_diff = []
    
    for e_id in w_exp.keys():
        fm_w = w_exp[e_id]
        fm_n = vec_neutrino[e_id]
        fm_n_built = vec_built[e_id]
        fm_l = vec_lepton[e_id]
        
        neutrino_direction = fm_n.to_3D()
        n_built_direction = fm_n_built.to_3D()
        lepton_direction   = fm_l.to_3D()
        
        cos_theta = compute_angle_between(neutrino_direction, lepton_direction)
        cos_theta_b = compute_angle_between(n_built_direction, lepton_direction)
        
        decay_cos_list.append(cos_theta)
        decay_cos_list_built.append(cos_theta_b)
        
        decay_angle_list.append( np.arccos(cos_theta) )
        decay_angle_list_built.append( np.arccos(cos_theta_b) )
        
        pz_diff.append(abs(fm_n.pz - fm_l.pz))
        pz_diff_built.append(abs(fm_n_built.pz - fm_l.pz))
        
        R_t_diff.append(fm_n_built.pz / fm_n.pz)
    
    #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    #plots
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(decay_cos_list, bins = 80, binrange= (-1, 1), color = 'steelblue',  element = 'step', label = 'expected', edgecolor = 'steelblue', stat = 'probability', alpha = 0.6)
    sb.histplot(decay_cos_list_built, bins = 80, binrange= (-1, 1), color = 'firebrick', element = 'step', label = 'computed', edgecolor = 'steelblue', stat = 'probability', alpha = 0.4)
    plt.xlabel("coseno")
    plt.ylabel("Event fraction") 
    plt.title("Cos distribution")  
    plt.legend()
    plt.show()
    
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(decay_angle_list, bins = 80, color = 'steelblue', element = 'step', edgecolor = 'steelblue', label = 'expected', stat = 'probability', alpha = 0.6)
    sb.histplot(decay_angle_list_built, bins = 80, color = 'firebrick', element = 'step', edgecolor = 'steelblue', label = 'computed', stat = 'probability', alpha = 0.4)
    plt.xlabel("Angle")
    plt.ylabel("Event fraction")
    plt.title("Angle distribution") 
    plt.legend()
    plt.show()
    
    
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(pz_diff, bins = 80, binrange = (0, 500),color = 'steelblue',  element = 'step', edgecolor = 'steelblue', stat = 'probability', label = 'expected', alpha = 0.6)
    sb.histplot(pz_diff_built, bins = 80, binrange = (0, 500),color = 'firebrick', element = 'step', edgecolor = 'steelblue', label = 'computed', stat = 'probability', alpha = 0.4)
    plt.xlabel("Pz Difference")
    plt.ylabel("Event fraction")  
    plt.title("Pz lepton-neutrino distribution")
    plt.legend()
    plt.show()
    
    plt.figure()
    sb.set(style = 'whitegrid')
    sb.histplot(R_t_angle, bins = 80, binrange = (0, 10),color = 'royalblue',  element = 'step', fill = False,  stat = 'probability', label = 'minimum angle', alpha = 0.6, linewidth = 2.5)
    sb.histplot(R_t_diff, bins = 80, binrange = (0, 10), color = 'firebrick', element = 'step', fill = False,  stat = 'probability',  label = 'pz difference', alpha = 0.4, linewidth = 2.5)
    plt.xlabel("Ratio")
    plt.ylabel("Event fraction")  
    plt.title("pz computed / pz expected")
    plt.legend()
    plt.show()
    

    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(cos_l_diff, bins = 70, binrange= (-1, 1), color = 'steelblue', edgecolor = 'steelblue', stat = 'probability', alpha = 0.6, label = 'longitudinal polarisation')
    sb.histplot(cos_t_diff, bins = 70,  binrange= (-1, 1), color = 'firebrick', edgecolor = 'steelblue', stat = 'probability', alpha = 0.4, label = 'transverse polarisation')
    plt.xlabel("cos(θ*)")
    plt.ylabel("Event fraction")
    plt.title("Angle")
    plt.legend()
    plt.show()
    
    plt.figure()
    sb.set(style='whitegrid')
    sb.histplot(cos_l_angle, bins = 70, binrange= (-1, 1), color = 'steelblue', edgecolor = 'steelblue', stat = 'probability', alpha = 0.6, label = 'longitudinal polarisation')
    sb.histplot(cos_t_angle, bins = 70,  binrange= (-1, 1), color = 'firebrick', edgecolor = 'steelblue', stat = 'probability', alpha = 0.4, label = 'transverse polarisation')
    plt.xlabel("cos(θ*)")
    plt.ylabel("Event fraction")   
    plt.title("Pz difference")
    plt.legend()
    plt.show()

    
if __name__ == '__main__':
    main()
