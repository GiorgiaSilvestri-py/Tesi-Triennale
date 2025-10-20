import vector
import pandas as pd
import random
import numpy as np

#classe Cinematic
class Cinematic:
    def __init__(self, pid, status, mother1, mother2, px, py, pz, E):
        self.pid = pid
        self.status = status
        self.mother1 = mother1
        self.mother2 = mother2
        self.px = px
        self.py = py
        self.pz = pz
        self.E = E

    def __str__(self):
        return (
            f"Particle Information:\n"
            #f"  pid     : {self.pid}\n"
            f"  status  : {self.status}\n"
            f"  mother1 : {self.mother1}\n"
            f"  mother2 : {self.mother2}\n"
            f"  px      : {self.px}\n"
            f"  py      : {self.py}\n"
            f"  pz      : {self.pz}\n"
            f"  E       : {self.E}"
        )
    
    def build_fm(self):
        return vector.obj(px = self.px, py = self.py, pz = self.pz, E = self.E)
    
    def compute_pt(self):
        return self.build_fm().pt
        
    def compute_eta(self):
        return self.build_fm().eta

    def compute_phi(self):
        return self.build_fm().phi

#dizionario che associa pid => Cinematic


def add_event(event_id, particles_dict, events_dict) :
    events_dict[event_id]= particles_dict



def add_particle_to_dict(pid, cinematic, particles_dict) :
    particles_dict[pid] = cinematic
    

def print_dict(events_dict):
    print(f" Events Dictionary: {len(events_dict)} events\n")
    
    for event_id, particles in events_dict.items():
        if event_id > 10:
            break
        
        print(f"Event ID: {event_id}")
        print("Particles:")
        
        for pid, cinematic in particles.items():
            print(f"\n  PID: {pid}")
            print(cinematic)
        
        print("\n" + "-"*40 + "\n")


def is_nonzero(v):
    return v.px != 0 or v.py != 0 or v.pz != 0 or v.E != 0

def compute_jet_resolution(energy):
    a0 = 5.43376
    a1 = 0.0593128
    
    return a0 / energy + a1

def apply_smearing_to(fm_vector, p_type):
    '''
    Apply smearing to a four-momentum vector according to particle type
    '''
    #neutrinos
    is_lepton = False
    is_jet = False
    is_neutrino = False
    
    if p_type == "lepton":
        is_lepton = True
        
    if p_type == "neutrino":
        is_lepton = True
        
    if p_type == "jet":
        is_lepton = True
        
    if is_neutrino:
        resolution = 0.20
        rand_smearing = random.gauss(mu = 1.0, sigma = resolution)
        px_new = fm_vector.px * rand_smearing
        py_new = fm_vector.py * rand_smearing
        pz_new = fm_vector.pz * rand_smearing
        E_new = fm_vector.E * rand_smearing
        smeared_vector = vector.obj( px = px_new, py = py_new, pz = py_new, E = E_new )
    
    if is_lepton:
        #leptons
        resolution = 0.01
        rand_smearing = random.gauss(mu = 1.0, sigma = resolution)
        smeared_vector = vector.obj( px = (fm_vector.px * rand_smearing), py = (fm_vector.py * rand_smearing), pz = (fm_vector.pz * rand_smearing), E =  (fm_vector.E * rand_smearing) ) 
    
    if is_jet:
        energy = fm_vector.E
        resolution = compute_jet_resolution(energy)
        rand_smearing = random.gauss(mu = 1.0, sigma = resolution)
        px_new = fm_vector.px * rand_smearing
        py_new = fm_vector.py * rand_smearing
        pz_new = fm_vector.pz * rand_smearing
        E_new = fm_vector.E * rand_smearing
        smeared_vector = vector.obj( px = px_new, py = py_new, pz = py_new, E = E_new )
        
    return smeared_vector
    
def solve_equation(a, b, c):
    
    delta = b**2 - 4*a*c
    
    if delta < 0:
        delta = 0
        #print("Negative delta, setting delta = 0")
        
    s1 = ( -b + np.sqrt(delta) ) / (2*a)
    s2 = ( -b - np.sqrt(delta) ) / (2*a)
    
    return s1, s2


def compute_neutrino_momentum(events_dict, e_id):
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
    '''
    #skip event if either one is None
    if not is_nonzero(vec_lepton) or not is_nonzero(vec_neutrino_expected):
        return vector.obj(px = 0, py = 0, pz = 0, E = 0)
    ''' 
    return compute_neutrino_momentum_from_particles(vec_lepton, vec_neutrino_expected) 



def compute_neutrino_momentum_from_particles(vec_lepton, vec_neutrino_expected):

    expected_pz_neutrino = vec_neutrino_expected.pz
    
    #solving second grade equation
    M = 80
    pt_lepton_dot_pt_neutrino = vec_lepton.px * vec_neutrino_expected.px + vec_lepton.py * vec_neutrino_expected.py
    pt_neutrino_squared = vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2
    
    alpha = M**2 + 2 * pt_lepton_dot_pt_neutrino
    a = 4 * (vec_lepton.pz ** 2) - 4 * vec_lepton.E ** 2
    
    b = 4 * alpha * vec_lepton.pz
    c = alpha**2 - 4 * vec_lepton.E**2 * pt_neutrino_squared
    
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


def get_fm_of(p_id, events_dict, apply_smearing = False):
    """
    Ricostruisce il quadrivettore di una particella -> fm_dict.
    - se p_id è Higgs (25) o leptoni (11, 13, 15), li legge direttamente
    - se p_id è Z (23) o W (±24), somma i quadrivettori dei figli leptoni
    - se neutrini pone pz=0
    """
    fm_dict = {}
    
    if apply_smearing:
        if set(p_id) == set([11, 13, 15]):
            for event_id, particles_dict in events_dict.items():
                for pid, cinematic in particles_dict.items():
                    if pid in p_id:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                        fm_dict[event_id] = smeared_momentum
                        
        if set(p_id) == set([-11, -13, -15]):
            for event_id, particles_dict in events_dict.items():
                for pid, cinematic in particles_dict.items():
                    if pid in p_id:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                        fm_dict[event_id] = smeared_momentum
        
        
        if set(p_id) == set([11, 12, 13, 14, 15, 16]):
            for event_id, particles_dict in events_dict.items():
                for pid, cinematic in particles_dict.items():
                    if pid in [11, 13, 15]:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                        fm_dict[event_id] = smeared_momentum
                    if pid in [12, 14, 16]:
                        non_smeared_momentum = compute_neutrino_momentum(events_dict, event_id)
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "neutrino")
                        fm_dict[event_id] = smeared_momentum
                        
        if set(p_id) == set([-11, -12, -13, -14, -15, -16]):
            for event_id, particles_dict in events_dict.items():
                for pid, cinematic in particles_dict.items():
                    if pid in [-11, -13, -15]:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                        fm_dict[event_id] = smeared_momentum
                    if pid in [-12, -14, -16]:
                        non_smeared_momentum = compute_neutrino_momentum(events_dict, event_id)
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "neutrino")
                        fm_dict[event_id] = smeared_momentum
                   
        if p_id == [23]:
        
            for event_id, particles_dict in events_dict.items():
                fm_vec_1 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                fm_vec_2 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                
                event_type = f_event_type(particles_dict)
                
                #24 = W; 23 = Z
                if abs(event_type) == 24 :
                    continue
                
                for pid, cinematic in particles_dict.items():
                    if (pid) in [11, 13]:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                        fm_vec_1 = smeared_momentum
                        
                    if (pid) in [-11, -13]:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                        fm_vec_2 = smeared_momentum
                        
                        
                if not is_nonzero(fm_vec_1) or not is_nonzero(fm_vec_2):
                    continue                                                    # salta l’evento se manca uno dei due
                    
               
                fm_dict[event_id] = fm_vec_1 + fm_vec_2
                

        if set(p_id) == set([24, -24]):
        
            for event_id, particles_dict in events_dict.items():
                event_type = f_event_type(particles_dict)
                
                #24 = W; 23 = Z
                if event_type == 23:
                    continue
                    
                fm_vec_1 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                fm_vec_2 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                
                for pid, cinematic in particles_dict.items():
                    if abs(pid) in [11, 13]:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                        fm_vec_1 = smeared_momentum
                        
                for pid, cinematic in particles_dict.items():
                    if abs(pid) in [12, 14]:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "neutrino")
                        fm_vec_2 = compute_neutrino_momentum_from_particles(fm_vec_1, smeared_momentum)
                        
                if not is_nonzero(fm_vec_1) or not is_nonzero(fm_vec_2):
                    continue                       # salta l’evento se manca uno dei due
                fm_dict[event_id] = fm_vec_1 + fm_vec_2
            '''
            for event_id, particles_dict in events_dict.items():
                event_type = f_event_type(particles_dict)
                
                #24 = W; 23 = Z
                if event_type == 23:
                    continue
                    
                vec_neutrino = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                vec_lepton = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                vec_built = vector.obj(px = 0, py = 0, pz = 0, E = 0)

                e_type = f_event_type(particles_dict)
        
                if e_type == 23:
                    continue
                    
                for pid, cinematic in particles_dict.items():
                    if abs(pid) in [12, 14, 16]:
                        vec_neutrino = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                        
                for pid, cinematic in particles_dict.items():     
                    if abs(pid) in [11, 13, 15]:
                        vec_lepton_lhe = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                        smeared_momentum = apply_smearing_to(vec_lepton_lhe, "lepton")
                        vec_lepton = smeared_momentum

                        
                eta_lepton = vec_lepton.eta
                
                if abs(eta_lepton) <= 2.5:
                    #usa minimum angle
                    vec_built, _ = compute_neutrino_momentum_from_particles(vec_lepton[evt_id], vec_neutrino[evt_id])
                        
                else:
                    #usa pz difference
                    vec_built = compute_neutrino_momentum_from_particles_diff(vec_lepton[evt_id], vec_neutrino[evt_id])
                
                if not is_nonzero(vec_lepton) or not is_nonzero(vec_built):
                    continue                       

                fm_w_dict[evt_id] = vec_lepton + vec_built
                '''
        
        if set(p_id) == set([25]):
            for event_id, particles_dict in events_dict.items():
                fm_vec_1 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                fm_vec_2 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                
                tau_list, antitau_list = [], []

                for pid, cinematic in particles_dict.items():
                    if pid == 15:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                        fm_vec_1 = smeared_momentum
                        tau_list.append(fm_vec_1)

                    if pid == -15:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                        fm_vec_2 = smeared_momentum
                        antitau_list.append(fm_vec_2)
                
                for fm_tau in tau_list:

                    for fm_anti_tau in antitau_list:
                        fm_sum = fm_tau + fm_anti_tau
                        tollerance = 10

                        if abs(fm_sum.M - 125) < tollerance:
                            fm_dict[event_id] = fm_sum
                
                
                '''
                for pid, cinematic in particles_dict.items():
                    if pid == 15 and cinematic.mother1 == 3:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                        fm_vec_1 = smeared_momentum
                        
                for pid, cinematic in particles_dict.items():
                    if pid == -15 and cinematic.mother1 == 3:
                        non_smeared_momentum = cinematic.build_fm()
                        smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                        fm_vec_2 = smeared_momentum
                
                fm_dict[event_id] = fm_vec_1 + fm_vec_2
                '''        
    else:
        if set(p_id) == set([11, 13, 15]):
            for event_id, particles_dict in events_dict.items():
                for pid, cinematic in particles_dict.items():
                    if pid in p_id:
                        fm_dict[event_id] = cinematic.build_fm()
                        
        if set(p_id) == set([-11, -13, -15]):
            for event_id, particles_dict in events_dict.items():
                for pid, cinematic in particles_dict.items():
                    if pid in p_id:
                        fm_dict[event_id] = cinematic.build_fm()
                        
        if set(p_id) == set([11, 12, 13, 14, 15, 16]):
            for event_id, particles_dict in events_dict.items():
                for pid, cinematic in particles_dict.items():
                    if pid in [11, 13, 15]:
                        fm_dict[event_id] = cinematic.build_fm()
                    if pid in [12, 14, 16]:
                        fm_dict[event_id] = compute_neutrino_momentum(events_dict, event_id)
                        
        if set(p_id) == set([-11, -12, -13, -14, -15, -16]):
            for event_id, particles_dict in events_dict.items():
                for pid, cinematic in particles_dict.items():
                    if pid in [-11, -13, -15]:
                        fm_dict[event_id] = cinematic.build_fm()
                    if pid in [-12, -14, -16]:
                        fm_dict[event_id] = compute_neutrino_momentum(events_dict, event_id)
         
        if p_id == [23]:
            for event_id, particles_dict in events_dict.items():
                event_type = f_event_type(particles_dict)
                
                #24 = W; 23 = Z
                if event_type == 24 or event_type == -24:
                    continue
                    
                fm_vec_1 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                fm_vec_2 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                
                for pid, cinematic in particles_dict.items():
                    if (pid) in [11, 13]:
                        fm_vec_1 = cinematic.build_fm()
                    if (pid) in [-11, -13]:
                        fm_vec_2 = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                           
                if not is_nonzero(fm_vec_1) or not is_nonzero(fm_vec_2):
                    continue                                                    # salta l’evento se manca uno dei due
                fm_dict[event_id] = fm_vec_1 + fm_vec_2
                
        if set(p_id) == set([24, -24]):
            for event_id, particles_dict in events_dict.items():
                event_type = f_event_type(particles_dict)
                
                #24 = W; 23 = Z
                if event_type == 23:
                    continue
                fm_vec_1 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                fm_vec_2 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
                
                for pid, cinematic in particles_dict.items():
                    if abs(pid) in [11, 13]:
                        fm_vec_1 = cinematic.build_fm()
                       
                    if abs(pid) in [12, 14]:
                        fm_vec_2 = compute_neutrino_momentum(events_dict, event_id)
                        
                if not is_nonzero(fm_vec_1) or not is_nonzero(fm_vec_2):
                    continue                                                    # salta l’evento se manca uno dei due
                fm_dict[event_id] = fm_vec_1 + fm_vec_2
                
                
        if set(p_id) == set([25]):
            for event_id, particles_dict in events_dict.items():
                
                tau_list, antitau_list = [], []

                for pid, cinematic in particles_dict.items():
                    if pid == 15:
                        tau_list.append(cinematic.build_fm())

                    if pid == -15:
                        antitau_list.append(cinematic.build_fm())
                
                for fm_tau in tau_list:

                    for fm_anti_tau in antitau_list:
                        fm_sum = fm_tau + fm_anti_tau
                        tollerance = 10

                        if abs(fm_sum.M - 125) < tollerance:
                            fm_dict[event_id] = fm_sum

    return fm_dict


def apply_selections(events_dict):

    selected_events = [] #list of event ids that pass selections
    
    for e_id, part_dict in events_dict.items():

        tau_list, antitau_list = [], []
        for pid, cinematic in part_dict .items():
            if pid == 15:
                non_smeared_momentum = cinematic.build_fm()
                smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                fm_vec_1 = smeared_momentum
                tau_list.append(fm_vec_1)

            if pid == -15:
                non_smeared_momentum = cinematic.build_fm()
                smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                fm_vec_2 = smeared_momentum
                antitau_list.append(fm_vec_2)
        
        for fm_tau in tau_list:

            for fm_anti_tau in antitau_list:
                fm_sum_0 = fm_tau + fm_anti_tau
                tollerance = 6

                if (fm_sum_0.M - 125) < tollerance:
                    fm_sum = fm_sum_0
                

        #selections
        pt_tau_tau = np.sqrt(fm_sum.px**2 + fm_sum.py**2)
        eta = fm_sum.eta
        
        if pt_tau_tau > 40 and abs(eta) < 2.1:
            selected_events.append(e_id)
                
    return selected_events
            
            

    
def f_event_type(particles_dict):
    
    event_type = 0  #24 = W; 23 = Z
    for pid, cinematic in particles_dict.items():
        if pid == 24:
            event_type = 24
        if pid == -24:
            event_type = -24
        elif pid == 23:
            event_type = 23
    return event_type

#-------------------------------------------------------------------------------------------------------

        
def get_pt_of(p_id, events_dict, apply_smearing = False):
    
    particle_fm = get_fm_of(p_id, events_dict, apply_smearing = apply_smearing)
    pt_dict = {}
    
    for event_id, p_fm in particle_fm.items():
        pt_dict[event_id] = p_fm.pt
    '''        
    for event_id, pt in pt_dict.items() :
        if event_id > 10:
            break
        print("Event ID: ", event_id, "\n\t Momento Trasverso: ", pt)
    '''   
    return pt_dict


    
def get_eta_of(p_id, events_dict, apply_smearing = False):
   
    particle_fm = get_fm_of(p_id, events_dict, apply_smearing = apply_smearing)
    eta_dict = {}
    
    for event_id, p_fm in particle_fm.items():
        eta_dict[event_id] = p_fm.eta
    '''
    for event_id, eta in eta_dict.items() :
        if event_id > 10:
            break
        print("Event ID: ", event_id, "\n\t Momento Trasverso: ", eta)
    ''' 
    return eta_dict
    
    
def get_phi_of(p_id, events_dict, apply_smearing = False):
    '''
    Crea la lista di eta di una particella (p_id)
    - scorre gli eventi
    - per ciascun evento cerca se c'è la particella
    - se esiste, calcola eta
    '''
    particle_fm = get_fm_of(p_id, events_dict, apply_smearing = apply_smearing)
    phi_dict = {}
    
    for event_id, p_fm in particle_fm.items():
        phi_dict[event_id] = p_fm.phi

    return phi_dict
    
    
def get_thetastar_of(V_fm_dict, H_fm_dict):
    '''
    Restituisce una lista di cos(theta*) per V con 
    - quad_V = dizionario di quadrimomenti VB; 
    - quad_H = dizionario di quadrimomenti H
    
    Procedimento:
    - calcola il quadrimomento somma VH
    - trova il vettore di boost
    - applica il boost a V e ad H
    - controlla che il SDR sia corretto (momento nullo)
    - identifica la direzione del fascio e della V nel SDR VH e calcola l'angolo
    '''

    cos_list = []
 
    
    common_ids = set.intersection(set(V_fm_dict), set(H_fm_dict))

    
    for event_id in common_ids:
        fm_v = V_fm_dict[event_id]
        fm_h = H_fm_dict[event_id]
        fm_tot = fm_v + fm_h

        #costruzione vettore di boost del singolo quadrimomento fm_tot
        #boost_vec = fm_v.boostCM_of_p4(fm_fh)
        #vector.obj(x = -fm_tot.px/fm_tot.E, y = -fm_tot.py/fm_tot.E, z = -fm_tot.pz/fm_tot.E)  
    
        #boost di V e H nel rest frame VH
        fm_v_rf = fm_v.boostCM_of_p4(fm_tot)
        fm_h_rf = fm_h.boostCM_of_p4(fm_tot)
        
        #controllo del boost
        fm_tot_rf = fm_v_rf + fm_h_rf
        if abs(fm_tot_rf.px) > 1e-5 or abs(fm_tot_rf.py) > 1e-5 or abs(fm_tot_rf.pz) > 1e-5 :
            print("Evento N°:", event_id, ". Trimomento non nullo: ", fm_tot_rf)
        
        #calcolo delle direzioni e dell'angolo
        z_axis = vector.obj(px = 0, py = 0, pz = 1, E=1)                                #asse z
        z_axis_rf = z_axis.boostCM_of_p4(fm_tot)
        #print(z_axis_rf)
        p_vec = vector.obj(x = fm_v_rf.px, y = fm_v_rf.py, z = fm_v_rf.pz)      #direzione di volo di V
        
        cos_theta = compute_angle(p_vec, z_axis_rf)
        cos_list.append((cos_theta))
        
    return cos_list
    
def get_deltaphi_of(V_dict, lep_dict, antilep_dict, H_dict):
    '''
    Calcola la separazione azimutale fra i leptoni prodotti dal decadimento nel SDR VH
    Procedimento:
    - ci si sposta nel sdr giusto tramite boost di lorentz
    - si applica il vettore di boost alla coppia leptone, antileptone
    - si calcola phi per leptone e per antileptone
    - si calcola la differenza
    '''
    deltaphi_list = []
    
    fm_VH = compute_tot_fm(V_dict, H_dict)                #dizionario di quadrimomenti totali

    for event_id in V_dict.keys() & H_dict.keys() & fm_VH.keys() & lep_dict.keys() & antilep_dict.keys():
        fm_v = V_dict[event_id]
        fm_h = H_dict[event_id]
        fm_lep = lep_dict[event_id]
        fm_antilep = antilep_dict[event_id]
        
        fm_tot = fm_v + fm_h

        #costruzione vettore di boost del singolo quadrimomento fm_tot
        boost_vec = vector.obj(x = -fm_tot.px/fm_tot.E, y = -fm_tot.py/fm_tot.E, z = -fm_tot.pz/fm_tot.E)  
        
        fm_v_rf = fm_v.boost(boost_vec)
        fm_lep_rf = fm_lep.boost(boost_vec)
        fm_antilep_rf = fm_antilep.boost(boost_vec)
        fm_h_rf = fm_h.boost(boost_vec)
        
        fm_VH_tot = fm_v_rf + fm_h_rf
        
        #controllo boost
        
        if abs(fm_VH_tot.px) > 1e-5 or abs(fm_VH_tot.py) > 1e-5 or abs(fm_VH_tot.pz) > 1e-5 :
            print("Evento N°:", event_id, ". Trimomento non nullo: ", fm_VH_tot)
        
        deltaphi_list.append(abs(fm_lep_rf.deltaphi(fm_antilep_rf)))
        
    return deltaphi_list
    
        
def get_theta1_of(V_dict, lep_dict, antilep_dict, H_dict):
    '''
    Calcola l'angolo theta [restituisce lista di cos] fra la direzione di volo del leptone (dal decadimento della V) e quella di H, nel SDR della V.
    Procedimento:
    - impone identità su id_evento
    - calcola il vettore di boost dal quadrimomento del VB
    - boosta leptone e higgs nel sdr del VB
    - controlla che il sdr sia corretto (il momento somma deve essere nullo)
    - calcola l'angolo
    '''

    cos_list = []
    
    common_ids = set.intersection(set(V_dict), set(lep_dict), set(antilep_dict), set(H_dict))

    
    for event_id in common_ids:
        fm_lep = lep_dict[event_id]
        fm_antilep = antilep_dict[event_id]
        fm_h = H_dict[event_id]
        
        fm_v = fm_lep + fm_antilep
        
        #boost_vec = vector.obj(x = -fm_v.px/fm_v.E, y = -fm_v.py/fm_v.E, z = -fm_v.pz/fm_v.E)
        
        fm_lep_rf = fm_lep.boostCM_of_p4(fm_v)
        fm_antilep_rf = fm_antilep.boostCM_of_p4(fm_v)
        fm_h_rf = fm_h.boostCM_of_p4(fm_v)
        
        #controllo boost
        fm_tot = fm_lep_rf + fm_antilep_rf
        if abs(fm_tot.px) > 1e-5 or abs(fm_tot.py) > 1e-5 or abs(fm_tot.pz) > 1e-5 :
            print("Evento N°:", event_id, ". Trimomento non nullo: ", fm_tot)
        
        lep_direction = vector.obj(x = fm_lep_rf.px, y = fm_lep_rf.py, z = fm_lep_rf.pz)
        h_direction = vector.obj(x = fm_h_rf.px, y = fm_h_rf.py, z = fm_h_rf.pz)
        
        theta_1 = compute_angle(lep_direction, h_direction)
        cos_list.append(abs(theta_1))
        
    return cos_list
    
def get_phi1_of(V_dict, lep_dict, antilep_dict, H_dict):
    '''
    Calcola l'angolo φ1 [restituisce lista di phi1] tra il piano del decadimento del bosone e il piano di scattering partonico nel sdr VH.
    Procedimento:
    - ci si sposta nel sdr VH
    - calcola il boost e boosta leptone, antileptone, V e H
    - calcola i versori e l'angolo
    '''
    
    phi1_list= []
    
    fm_tot_in = compute_tot_fm(V_dict, H_dict)
    #print(fm_tot_in)

    for event_id in V_dict.keys() & lep_dict.keys() & antilep_dict.keys() & H_dict.keys():
        fm_v = V_dict[event_id]
        fm_lep = lep_dict[event_id]
        fm_antilep = antilep_dict[event_id]
        fm_h = H_dict[event_id]
        z_axis = vector.obj(px = 0, py = 0, pz = 1, E = 1)
        
        fm_tot_i = fm_v + fm_h #posso sommarli perché impongo già che l'id evento sia lo stesso
        
        #boost_vec = vector.obj(x = -fm_tot_i.px/fm_tot_i.E, y = -fm_tot_i.py/fm_tot_i.E, z = -fm_tot_i.pz/fm_tot_i.E)
        
        fm_lep_rf = fm_lep.boostCM_of_p4(fm_tot_i)
        fm_v_rf = fm_v.boostCM_of_p4(fm_tot_i)
        fm_antilep_rf = fm_antilep.boostCM_of_p4(fm_tot_i)
        fm_h_rf = fm_h.boostCM_of_p4(fm_tot_i)
        z_ax_r = z_axis.boostCM_of_p4(fm_tot_i)
        z_ax_rf = z_ax_r.to_3D()
        
        #controllo boost
        fm_tot = fm_v_rf + fm_h_rf
        if abs(fm_tot.px) > 1e-5 or abs(fm_tot.py) > 1e-5 or abs(fm_tot.pz) > 1e-5 :
            print("Evento N°:", event_id, ". Trimomento non nullo: ", fm_tot)
        
        #print(fm_tot)
        
        pL_vec = fm_lep_rf.to_3D()
        p_aL_vec = fm_antilep_rf.to_3D()
        p_v_vec = fm_v_rf.to_3D()
        
        pl = np.array([pL_vec.x, pL_vec.y, pL_vec.z])
        plbar = np.array([p_aL_vec.x, p_aL_vec.y, p_aL_vec.z])
        pv = np.array([p_v_vec.x, p_v_vec.y, p_v_vec.z])
        nz = np.array([z_ax_rf.x, z_ax_rf.y, z_ax_rf.z])
        
        cross = np.cross(pl, plbar)
        nv = cross / np.linalg.norm(cross)
        #nV = pL_vec.cross(p_aL_vec).unit()
        
        cross = np.cross(nz, pv)
        nsc = cross / np.linalg.norm(cross)
        #n_sc = z_ax_rf.cross(p_z_vec).unit()
        
        num = np.dot(pv, np.cross(nv, nsc))
        phi1 = np.sign(num) * np.arccos( np.dot(nv, nsc) )
         
        #phi1 = np.sign(p_z_vec @ nV.cross(n_sc)) * np.arccos( (nV @ n_sc) )
        phi1_list.append(abs(phi1))
   
    return phi1_list
    
    
def get_lep_ptbalance(events_dict, apply_smearing = False):
    '''
    calcola il lepton pt balance: 
    Argomenti: dizionario di eventi totali
    Output: liste di pt_balance per z e per W
    '''
    
    pt_balance_list_z, pt_balance_list_w = [], []
    
    #dizionari necessari
    z_fm = get_fm_of([23], events_dict, apply_smearing = apply_smearing)
    zlep, zantilep, _, _ = build_decay_products(events_dict, apply_smearing = apply_smearing)
    
    for event_id in z_fm.keys() & zlep.keys() & zantilep.keys() :
        p4z = z_fm[event_id]
        p4lep = zlep[event_id]
        p4antilep = zantilep[event_id]
        
        ptlep = p4lep.pt
        ptantilep = p4antilep.pt
        
        pt_sum = ptlep + ptantilep
        
        if pt_sum == 0:
            continue
            
        pt_balance = abs(ptlep - ptantilep) / pt_sum
    
        pt_balance_list_z.append(pt_balance)
        
    w_fm = get_fm_of([24, -24], events_dict)
    _, _, wlep, wantilep = build_decay_products(events_dict, apply_smearing = apply_smearing)
    
    
    for event_id in w_fm.keys() & wlep.keys() & wantilep.keys():
        p4w = w_fm[event_id]
        p4lep = wlep[event_id]
        p4antilep = wantilep[event_id]
        
        ptlep = p4lep.pt
        ptantilep = p4antilep.pt
        #print(ptlep, ptantilep)
        
        pt_sum = ptlep + ptantilep
        
        if pt_sum == 0:
            continue
            
        pt_balance = abs(ptlep - ptantilep) / pt_sum
    
        pt_balance_list_w.append(pt_balance)

    return pt_balance_list_z, pt_balance_list_w
    


def compute_tot_fm(fm1_dict, fm2_dict):
    '''
    Calcola il quadrimomento somma
    '''
    fm_tot_dict = {}
    
    for event_id in fm1_dict.keys() & fm2_dict.keys():
        fm_1 = fm1_dict[event_id]
        fm_2 = fm2_dict[event_id]

        fm_tot_dict[event_id] = fm_1 + fm_2

    return fm_tot_dict



def build_decay_products(events_dict, apply_smearing = False):

    
    lepton_dict_z = {}
    antilepton_dict_z = {}
    
    antilepton_dict_w_plus = {}
    lepton_dict_w_plus = {}
    
    antilepton_dict_w_minus = {}
    lepton_dict_w_minus = {}         
    
    if apply_smearing:
        for event_id, p_dict in events_dict.items():
            #building Z dictionary
            event_type = f_event_type(p_dict)    #24 = W; 23 = Z
            if event_type == 23:
              
                for pid, cinematic in p_dict.items():
                    if pid in [11, 13]:
                        non_smeared_fm = cinematic.build_fm()
                        lepton_dict_z[event_id] = apply_smearing_to(non_smeared_fm, "lepton")
                    if pid in [-11, -13]:
                        non_smeared_fm = cinematic.build_fm()
                        antilepton_dict_z[event_id] = apply_smearing_to(non_smeared_fm, "lepton")
            
            #building W+ dictionary --> l+n
            if event_type == 24:
                for pid, cinematic in p_dict.items():
                    if pid in [-11, -13]:
                        antilepton_vec_non_smeared = cinematic.build_fm()
                        antilepton_vec_smeared = apply_smearing_to(antilepton_vec_non_smeared, "lepton")
                        antilepton_dict_w_plus[event_id] = antilepton_vec_smeared
                for pid, cinematic in p_dict.items():
                    if pid in [12, 14]:
                        vec_neutrino_expected = cinematic.build_fm()
                        vec_neutrino_smeared = apply_smearing_to(vec_neutrino_expected, "neutrino")
                        lepton_dict_w_plus[event_id] = compute_neutrino_momentum_from_particles(antilepton_vec_smeared, vec_neutrino_smeared)
            
            #building W- dictionary --> l-n_bar
            if event_type == -24:
                for pid, cinematic in p_dict.items():
                    if pid in [11, 13]:
                        lepton_vec = cinematic.build_fm()
                        lepton_vec_smeared = apply_smearing_to(lepton_vec, "lepton")
                        lepton_dict_w_minus[event_id] = lepton_vec_smeared
                for pid, cinematic in p_dict.items():
                    if pid in [-12, -14]:
                        antineutrino_vec = cinematic.build_fm()
                        antineutrino_smeared = apply_smearing_to(antineutrino_vec, "neutrino")
                        antilepton_dict_w_minus[event_id] = compute_neutrino_momentum_from_particles(lepton_vec_smeared, antineutrino_smeared)
        
        
    else:
        for event_id, p_dict in events_dict.items():
            #building Z dictionary
            event_type = f_event_type(p_dict)    #24 = W; 23 = Z
            if event_type == 23:
              
                for pid, cinematic in p_dict.items():
                    if pid in [11, 13]:
                        lepton_dict_z[event_id] = cinematic.build_fm()
                    if pid in [-11, -13]:
                        antilepton_dict_z[event_id] = cinematic.build_fm()
            
            #building W+ dictionary --> l+n
            if event_type == 24:
                for pid, cinematic in p_dict.items():
                    if pid in [-11, -13]:
                        antilepton_vec = cinematic.build_fm()
                        antilepton_dict_w_plus[event_id] = antilepton_vec
                for pid, cinematic in p_dict.items():
                    if pid in [12, 14]:
                        vec_neutrino_expected = cinematic.build_fm()
                        lepton_dict_w_plus[event_id] = compute_neutrino_momentum_from_particles(antilepton_vec, vec_neutrino_expected)
            
            #building W- dictionary --> l-n_bar
            if event_type == -24:
                for pid, cinematic in p_dict.items():
                    if pid in [11, 13]:
                        lepton_vec = cinematic.build_fm()
                        lepton_dict_w_minus[event_id] = lepton_vec
                for pid, cinematic in p_dict.items():
                    if pid in [-12, -14]:
                        antineutrino_vec = cinematic.build_fm()
                        antilepton_dict_w_minus[event_id] = compute_neutrino_momentum_from_particles(lepton_vec, antineutrino_vec)
                    
    lepton_dict_w = {}
    antilepton_dict_w = {}
    
    lepton_dict_w.update(lepton_dict_w_plus)
    lepton_dict_w.update(lepton_dict_w_minus)
    
    antilepton_dict_w.update(antilepton_dict_w_plus)
    antilepton_dict_w.update(antilepton_dict_w_minus)
    
    return lepton_dict_z, antilepton_dict_z, lepton_dict_w, antilepton_dict_w


   
def connect_lep_to_V(events_dict):
  
    
    lep_from_Z, lep_from_W = {}, {}
    antilep_from_Z, antilep_from_W = {}, {}
    
    z_fm = get_fm_of([23], events_dict)
    w_fm = get_fm_of([24, -24], events_dict)
    
    lep_dict = get_fm_of([11,  13, 15], events_dict)
    anti_lep_dict = get_fm_of([-11, -13, -15], events_dict)
    
    for event_id, lep_fm in lep_dict.items():
        if event_id in z_fm.keys():
            lep_from_Z[event_id] = lep_fm
        if event_id in w_fm.keys():
            lep_from_W[event_id] = lep_fm
        
            
    for event_id, antilep_fm in anti_lep_dict.items():
        if event_id in z_fm.keys():
            antilep_from_Z[event_id] = antilep_fm
        if event_id in w_fm.keys():
            antilep_from_W[event_id] = antilep_fm
    
    return lep_from_Z, antilep_from_Z, lep_from_W, antilep_from_W
    

def connect_alllep_to_V(events_dict):
 
    
    lep_from_Z, lep_from_W = {}, {}
    antilep_from_Z, antilep_from_W = {}, {}
    
    z_fm = get_fm_of([23], events_dict)
    w_fm = get_fm_of([24, -24], events_dict)
    
    lep_dict = get_fm_of([11, 12, 13, 14, 15, 16], events_dict)
    anti_lep_dict = get_fm_of([-11, -12, -13, -14,  -15, -16], events_dict)
    
    
    for event_id, lep_fm in lep_dict.items():
        if event_id in z_fm.keys():
            lep_from_Z[event_id] = lep_fm
        if event_id in w_fm.keys():
            lep_from_W[event_id] = lep_fm
        
            
    for event_id, antilep_fm in anti_lep_dict.items():
        if event_id in z_fm.keys():
            antilep_from_Z[event_id] = antilep_fm
        if event_id in w_fm.keys():
            antilep_from_W[event_id] = antilep_fm
        
    
    return lep_from_Z, antilep_from_Z, lep_from_W, antilep_from_W


def read_file(nome_file):

    with open(nome_file, 'r') as f:
        in_event = False
        event_lines = []
        events_dict = {}
        event_id = 0

        for line in f:
            line = line.strip()
            
            if "<event>" in line:
                in_event = True
                event_lines = []
            elif "</event>" in line:
                event_id += 1
                particles_dict = {}
                for data in event_lines:
                    parameter = data.strip().split()
                    if len(parameter) > 9:
                        pid = int(parameter[0])
                        cinematic = Cinematic(
                                            pid,
                                            int(parameter[1]),
                                            int(parameter[2]),
                                            int(parameter[3]),
                                            float(parameter[6]),
                                            float(parameter[7]),
                                            float(parameter[8]),
                                            float(parameter[9])
                                            )
                        add_particle_to_dict(pid, cinematic, particles_dict)
                        
                in_event = False
                add_event(event_id, particles_dict, events_dict)
                
            else :
                if in_event:
                    event_lines.append(line)
                        
        #print_dict(events_dict)
        
        
    return events_dict
 
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def contains_particle (particle_list, particle_ID) :
    for p in particle_list :
        if (p["pid"]) == particle_ID :
            return True

#------------------------------------------------------------------------------------------------------------------------------

 
def compute_angle(v1, v2):
    """
    Calcola il coseno dell'angolo tra due vettori tridimensionali
    """
    p1 = np.array([v1.x, v1.y, v1.z])
    p2 = np.array([v2.x, v2.y, v2.z])

    cos_theta = np.dot(p1, p2) / (np.linalg.norm(p1) * np.linalg.norm(p2))
    
    return np.clip(cos_theta, -1.0, 1.0)
    
    
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def compute_angle_between(v1, v2):
    """
    Calcola il coseno dell'angolo tra due vettori tridimensionali
    """
    p1 = np.array([v1.x, v1.y, v1.z])
    p2 = np.array([v2.x, v2.y, v2.z])

    cos_theta = np.dot(p1, p2) / (np.linalg.norm(p1) * np.linalg.norm(p2))
    
    return cos_theta
  
    
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


def build_df_dict_Z(LHE_file, apply_smearing = False):
    '''
    Costruisce un dataframe per la Z a partire dal file inserito come argomento
    '''

    pt_l1, pt_l2, eta_l1, eta_l2, phi_l1, phi_l2, pt_b_z = [], [], [], [], [], [], []
    pt_z, pt_h, eta_z, eta_h, phi_z, phi_h, M_zh = [], [], [], [], [], [], []
    cos_star, cos_1, phi_1, delta_phi = [], [], [], []
    
    #costruzione dizionari
    events_complete = read_file(LHE_file) #all events
    id_events = apply_selections(events_complete) #event ids that pass selections
    print(len(id_events))
    events = {}
    
    for i in id_events:
        events[i] = events_complete[i]
    
    
    fm_z_dict_all = get_fm_of([23], events, apply_smearing = apply_smearing)                                                    #costruisce il dizionario di quadrimomenti della Z
    lep_z_dict, antilep_z_dict, _, _ = build_decay_products(events, apply_smearing = apply_smearing)                                 #dizionario di quadrimomenti di leptone e antileptone dalla Z
    fm_h_all_dict = get_fm_of([25], events, apply_smearing = apply_smearing)                                                    #dizionario con tutti gli eventi, da filtrare su Z
    
    #filtro Z
    fm_z_dict = {}
    for event_id in lep_z_dict.keys() & fm_z_dict_all.keys():
        fm_z     = fm_z_dict_all[event_id]
        fm_lep = lep_z_dict[event_id]        
        
        fm_z_dict[event_id] = fm_z
    
    #filtro H
    fm_h_dict = {}                                                                           #dizionario di fm H solo in presenza di z
    
    for event_id in fm_z_dict.keys() & fm_h_all_dict.keys():
        fm_z     = fm_z_dict[event_id]
        fm_higgs = fm_h_all_dict[event_id]        
        
        fm_h_dict[event_id] = fm_higgs
    
    
    #costruzione liste di variabili
    print("Building lists")
    
    for ((e_id1, fm_lep), (e_id2, fm_antilep)) in zip(lep_z_dict.items(), antilep_z_dict.items()):
        pt_l1.append(fm_lep.pt)
        pt_l2.append(fm_antilep.pt)
        
        eta_l1.append(fm_lep.eta)
        eta_l2.append(fm_antilep.eta)
        
        phi_l1.append(fm_lep.phi)
        phi_l2.append(fm_antilep.phi)
        
    pt_b_z, _ = get_lep_ptbalance(events)
    
    pt_z  = [fm_z.pt for fm_z in list(fm_z_dict.values())]
    #pt_h  = [fm_h.pt for fm_h in list(fm_h_dict.values())]
    eta_z = [fm_z.eta for fm_z in list(fm_z_dict.values())]
    #eta_h = [fm_h.eta for fm_h in list(fm_h_dict.values())]
    phi_z = [fm_z.phi for fm_z in list(fm_z_dict.values())]
    #phi_h = [fm_h.phi for fm_h in list(fm_h_dict.values())]
    M_zh  = [fm_tot.M for fm_tot in list(compute_tot_fm(fm_z_dict, fm_h_dict).values())]
    
    cos_star  = get_thetastar_of(fm_z_dict, fm_h_dict)
    cos_1     = get_theta1_of(fm_z_dict, lep_z_dict, antilep_z_dict, fm_h_dict)
    phi_1     = get_phi1_of(fm_z_dict, lep_z_dict, antilep_z_dict, fm_h_dict)
    delta_phi = get_deltaphi_of(fm_z_dict, lep_z_dict, antilep_z_dict, fm_h_dict)

    #creazione dataframe
    print("Building dataframe")
    log = {
            "Pt lepton"         : pt_l1,
            "Pt antilepton"     : pt_l2,
            "Eta lepton"        : eta_l1,
            "Eta antilepton"    : eta_l2,
            "Phi lepton"        : phi_l1,
            "Phi antilepton"    : phi_l2,
            "Pt Z"               : pt_z,
            #"Pt H"               : pt_h,
            "Eta Z"              : eta_z,
            #"Eta H"              : eta_h,
            "Phi Z"              : phi_z,
            #"Phi H"              : phi_h,
            "ZH invariant mass"  :  M_zh,
            "Pt balance"         : pt_b_z,
            "cos(θ*)"            : cos_star,
            "cos(θ1)"            : cos_1,
            "φ1"                 : phi_1,
            "Δφ"                 : delta_phi
            }

    
    for (name, list_v) in log.items():
        print(name, len(list_v))

    df = pd.DataFrame(log)
    
    df_dict = {
               'events': df.values,
               'names': df.columns.tolist()
              }
    print(df.shape)
    print(df.head())
    
    
    return(df_dict)

def build_df_dict_W(LHE_file, apply_smearing = False):
    '''
    Costruisce un dataframe per la W a partire dal file inserito come argomento
    '''

    pt_l1, pt_l2, eta_l1, eta_l2, phi_l1, phi_l2, pt_b_w = [], [], [], [], [], [], []
    pt_w, pt_h, eta_w, eta_h, phi_w, phi_h, M_wh = [], [], [], [], [], [], []
    cos_star, cos_1, phi_1, delta_phi = [], [], [], []
    
    #costruzione dizionari
    events_complete = read_file(LHE_file)           #all events
    id_events_tot = apply_selections(events_complete)   #event ids that pass selections
    print(len(id_events_tot))
    events = {}
    
    #selected events
    for i in id_events_tot:
        events[i] = events_complete[i]
    
    fm_w_dict = get_fm_of([24, -24], events, apply_smearing = apply_smearing)
    print("len W", len(fm_w_dict.values()))  
    _, _, lep_w_dict, antilep_w_dict = build_decay_products(events, apply_smearing = apply_smearing) 
    fm_h_all_dict = get_fm_of([25], events, apply_smearing = apply_smearing)                               
    
    print("len lep, antilep", len(lep_w_dict.values()), len(antilep_w_dict.values()))
    
    #filtro H
    fm_h_dict = {}           
    
    for event_id in fm_w_dict.keys() & fm_h_all_dict.keys():
        fm_w     = fm_w_dict[event_id]
        fm_higgs = fm_h_all_dict[event_id]        
        
        fm_h_dict[event_id] = fm_higgs

    print("len H",len(fm_h_all_dict.values()))
    
    
    #costruzione liste di variabili
    print("Building lists")
    common_ids = set.intersection(set(fm_w_dict), set(lep_w_dict), set(antilep_w_dict), set(fm_h_all_dict))
    
    for event_id in common_ids:
        fm_w = fm_w_dict[event_id]
        fm_lep = lep_w_dict[event_id]
        fm_antilep = antilep_w_dict[event_id]
        fm_h = fm_h_all_dict[event_id]
        fm_vh = fm_w + fm_h
    
        pt_l1.append(fm_lep.pt)
        pt_l2.append(fm_antilep.pt)
        
        eta_l1.append(fm_lep.eta)
        eta_l2.append(fm_antilep.eta)
        
        phi_l1.append(fm_lep.phi)
        phi_l2.append(fm_antilep.phi)
        
        pt_w.append(fm_w.pt)
        #pt_h.append(fm_h.pt)
        
        eta_w.append(fm_w.eta)
        #eta_h.append(fm_h.eta)
        
        phi_w.append(fm_w.phi)
        #phi_h.append(fm_h.phi)
        
        #M_wh.append(fm_vh.M)
        
    _, pt_b_w = get_lep_ptbalance(events)
    
    
    #pt_w  = [fm_w.pt for fm_w in list(fm_w_dict.values())]
    #pt_h  = [fm_h.pt for fm_h in list(fm_h_dict.values())]
    #eta_w = [fm_w.eta for fm_w in list(fm_w_dict.values())]
    #eta_h = [fm_h.eta for fm_h in list(fm_h_dict.values())]
    #phi_w = [fm_w.phi for fm_w in list(fm_w_dict.values())]
    #phi_h = [fm_h.phi for fm_h in list(fm_h_dict.values())]
    M_wh  = [fm_tot.M for fm_tot in list(compute_tot_fm(fm_w_dict, fm_h_dict).values())]
    
    cos_star  = get_thetastar_of(fm_w_dict, fm_h_dict)
    cos_1     = get_theta1_of(fm_w_dict, lep_w_dict, antilep_w_dict, fm_h_dict)
    phi_1     = get_phi1_of(fm_w_dict, lep_w_dict, antilep_w_dict, fm_h_dict)
    delta_phi = get_deltaphi_of(fm_w_dict, lep_w_dict, antilep_w_dict, fm_h_dict)

    #creazione dataframe
    print("Building dataframe")
    log = {
            "Pt lepton"         : pt_l1,
            "Pt antilepton"     : pt_l2,
            "Eta lepton"        : eta_l1,
            "Eta antilepton"    : eta_l2,
            "Phi lepton"        : phi_l1,
            "Phi antilepton"    : phi_l2,
            "Pt W"               : pt_w,
            #"Pt H"               : pt_h,
            "Eta W"              : eta_w,
            #"Eta H"              : eta_h,
            "Phi W"              : phi_w,
            #"Phi H"              : phi_h,
            "WH invariant mass": M_wh,
            "Pt balance"         : pt_b_w,
            "cos(θ*)"            : cos_star,
            "cos(θ1)"            : cos_1,
            "φ1"                 : phi_1,
            "Δφ"                 : delta_phi
            }

    for (name, list_v) in log.items():
        print(name, len(list_v))

    df = pd.DataFrame(log)
    
    print(df.shape)
    print(df.head())
    
    df_dict = {
               'events': df.values,
               'names': df.columns.tolist()
              }
    
    
    return(df_dict)

#--------------------------------------------------------------------------------------------------------------------------------
from math import ceil
def sturges (N_events) :
    return ceil (1 + 3.322 * np.log (N_events))

import ROOT
def ROOT_hist1d(list1, list2, var, axis_name, nbins = None, xmin = None, xmax = None, ylim = None, xlim = None):

    if type(list1) == dict:
        list1 = list(list1.values())
    
    if type(list2) == dict:
        list2 = list(list2.values())

    if xmin is None:
        xmin1 = min(list1)
        xmin2 = min(list2)
    else:
        xmin1 = xmin
        xmin2 = xmin
        
    if xmax is None:
        xmax1 = max(list1)
        xmax2 = max(list2)
    else:
        xmax1 = xmax
        xmax2 = xmax

    if nbins is None:
        nbins1 =  sturges (len (list1))
        nbins2 =  sturges (len (list2))
    else:
        nbins1, nbins2 = nbins

    hist1 = ROOT.TH1F("hist1", f"{var} distribution", nbins1, xmin1, xmax1)
    hist2 = ROOT.TH1F("hist2", f"{var} distribution", nbins2, xmin2, xmax2)
    
    for value in list1:
        hist1.Fill(value)
    for value in list2:
        hist2.Fill(value)
        
    #normalisation
    if hist1.Integral() > 0:
        hist1.Scale(1.0/hist1.Integral())
    
    if hist2.Integral() > 0:
        hist2.Scale(1.0/hist2.Integral())
        
    steelblue  = ROOT.TColor.GetColor("#4682B4")  #steelblue
    firebrick  = ROOT.TColor.GetColor("#B22222")  #firebrick

    hist1.SetLineColor(steelblue)
    hist2.SetLineColor(firebrick)
    hist1.SetLineWidth(2)
    hist2.SetLineWidth(2)
    
    hist1.SetLineColor(steelblue)
    hist1.SetFillColor(steelblue)
    hist1.SetFillStyle(3002)
    
    hist2.SetLineColor(firebrick)
    hist2.SetFillColor(firebrick)
    hist2.SetFillStyle(3002)

    if ylim is not None:
        hist1.GetYaxis().SetRangeUser(0, ylim)
        hist2.GetYaxis().SetRangeUser(0, ylim)

    if xlim is not None:
        if isinstance(xlim, tuple) and len(xlim) == 2:
            xlim_m, xlim_M = xlim
            hist1.GetXaxis().SetRangeUser(xlim_m, xlim_M)
            hist2.GetXaxis().SetRangeUser(xlim_m, xlim_M)
        else:
            raise ValueError("xlim must be a two element tuple: (min, max)")


    hist1.GetXaxis().SetTitle(f"{axis_name}")
    hist1.GetYaxis().SetTitle("Event Fraction")
    
    ROOT.gROOT.SetBatch(True)
    
    canvas = ROOT.TCanvas("canvas", f"{var} distribution", 1600, 1200)
    hist1.Draw("HIST")
    hist2.Draw("HIST SAME")
    
    legend = ROOT.TLegend(0.6, 0.8, 0.89, 0.89)
    legend.SetFillStyle(0)                
    legend.SetTextFont(42)                
    legend.SetTextSize(0.03)              
    legend.AddEntry(hist1, "longitudinal polarisation", "f")
    legend.AddEntry(hist2, "transverse polarisation", "f")
    legend.Draw()

    canvas.SetFillColor(0)                
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.08)
    
    ROOT.gStyle.SetOptStat(0)
    
    canvas.SaveAs(f"{var} distribution.png")
