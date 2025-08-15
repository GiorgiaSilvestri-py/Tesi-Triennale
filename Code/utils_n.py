import vector
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


def get_fm_of(p_id, events_dict):
    '''
    Crea la lista di quadrimomenti di una particella (p_id)
    - scorre gli eventi
    - per ciascun evento cerca se c'è la particella
    - se esiste, costruisce il quadrimomento
    '''
    fm_dict = {}
    
    for event_id, particles_dict in events_dict.items():
        for pid, cinematic in particles_dict.items():
            if pid == p_id :
                fm_dict[event_id] = cinematic.build_fm()
            
    for event_id, fm in fm_dict.items() :
        if event_id > 10:
            break
        print("Event ID: ", event_id, "\n\t Quadrimomento: ", fm)
    
    return fm_dict
        
def get_pt_of(p_id, events_dict):
    '''
    Crea la lista di pt di una particella (p_id)
    - scorre gli eventi
    - per ciascun evento cerca se c'è la particella
    - se esiste, calcola il pt
    '''
    pt_dict = {}
    
    for event_id, particles_dict in events_dict.items():
        for pid, cinematic in particles_dict.items():
            if pid == p_id :
                pt_dict[event_id] = cinematic.compute_pt()
            
    for event_id, pt in pt_dict.items() :
        if event_id > 10:
            break
        print("Event ID: ", event_id, "\n\t Momento Trasverso: ", pt)
        
    return pt_dict
    
def compute_tot_fm(fm1_dict, fm2_dict):
    '''
    Calcola il quadrimomento somma
    '''
    fm_tot_dict = {}
    
    for event_id1, fm_1 in fm1_dict.items() :
        for event_id2, fm_2 in fm2_dict.items() :
            if event_id1 == event_id2 :
                if event_id1 > 10:
                    break
                fm_tot_dict[event_id1] = fm_1 + fm_2
            
    for event_id, fm_tot in fm_tot_dict.items() :
        if event_id > 10:
            break
        print("Event ID: ", event_id, "\n\t Quadrimomento Somma: ", fm_tot)
        
    return fm_tot_dict
    

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
                        
        print_dict(events_dict)
        
        W_fm = get_fm_of(24, events_dict)
        W_pt = get_pt_of(24, events_dict) 
        
        H_fm = get_fm_of(25, events_dict)
        compute_tot_fm(W_fm, H_fm)
        
        
        
    return events_dict
 
#--------------------------------------------------------------------------------------------------------------------------------

def contains_particle (particle_list, particle_ID) :
    for p in particle_list :
        if (p["pid"]) == particle_ID :
            return True
            
#--------------------------------------------------------------------------------------------------------------------------------
            
def build_fm_ZW (events_list) :
    '''
    restituisce liste di (vector.obj) quadrivettori della particella (W/Z) sommando quelli dei suoi prodotti di decadimento
    '''
    qvectorZ_list, qvectorW_list = [], []
    
    events_Z = [e for e in events_list if contains_particle (e, 23)]
    events_W = [e for e in events_list if contains_particle (e, 24) or contains_particle (e, -24)]
    
    for event in events_Z:
        fm_vec = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        
        for p in event :
            if p["status"] == 1 and abs(p["pid"]) in [11, 12, 14, 16, 13, 15] :
                fm_vec += vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])     #quadrivettore somma
        
        qvectorZ_list.append(fm_vec)
        
    
    for event in events_W:
        fm_vec = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        
        for p in event :
            if p["status"] == 1 and abs(p["pid"]) in [11, 13, 15] :
                fm_vec += vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])      #quadrivettore somma
        
        qvectorW_list.append(fm_vec)
        
    return qvectorZ_list, qvectorW_list

#------------------------------------------------------------------------------------------------------------------------------

def build_decay_lists (events_list) :
    '''
    restituisce liste di (vector.obj) quadrivettori della particella (Z/W) e liste di quadrivettori dei leptoni/antileptoni prodotti dal decadimento della Z
    '''
    
    qvectorZ_list, qvectorW_list = [], []
    fm_lepZ_list, fm_lepW_list = [], []
    fm_antilepZ_list, fm_antilepW_list = [], []
    
    events_Z = [e for e in events_list if contains_particle (e, 23)]
    events_W = [e for e in events_list if contains_particle (e, 24) or contains_particle (e, -24)]
    
    for event in events_Z:
        fm_vec = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        fm_lep_Z = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        fm_antilep_Z = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        
        for p in event :
            if p["status"] == 1 and abs(p["pid"]) in [11, 12, 14, 16, 13, 15] :
                fm_vec += vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])     #quadrivettore somma
            if p["status"] == 1 and p["pid"] in [11, 12, 14, 16, 13, 15] :
                fm_lep_Z = vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])     #quadrivettore leptone
            if p["status"] == 1 and p["pid"] in [-11, -12, -14, -16, -13, -15] :
                fm_antilep_Z = vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])     #quadrivettore antileptone
                            
        qvectorZ_list.append(fm_vec)
        fm_lepZ_list.append(fm_lep_Z)
        fm_antilepZ_list.append(fm_antilep_Z)

    for event in events_W:
        fm_vec = vector.obj(px = 0, py = 0, pz = 0, E = 0)
        
        for p in event :
            if p["status"] == 1 and abs(p["pid"]) in [11, 13, 15] :
                fm_vec += vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])      #quadrivettore somma
            
        qvectorW_list.append(fm_vec)
        
        
    return qvectorZ_list, qvectorW_list, fm_lepZ_list, fm_antilepZ_list
    
#------------------------------------------------------------------------------------------------------------------------------

def boost_to_rf(quad_vec_1, quad_vec_2) :
        '''
        restituisce le liste dei quadrimomenti nel sistema di riferimento del cm:
        quad_vec_1 : lista di quadrimomenti della prima particella
        quad_vec_2 : lista di quadrimomenti della seconda particella
        '''
        list_fm_12_rf = []
        list_1_rf, list_2_rf = [], []
        
        for (fm_1, fm_2) in zip(quad_vec_1, quad_vec_2) :
            fm_12 = fm_1 + fm_2
            
            bx = fm_12.px / fm_12.E
            by = fm_12.py / fm_12.E
            bz = fm_12.pz / fm_12.E
            boost_vec = vector.obj(x=-bx, y=-by, z=-bz)
        
            fm_1_rf = fm_1.boost(boost_vec)
            fm_2_rf = fm_2.boost(boost_vec)
            fm_12_rf = fm_1_rf + fm_2_rf
            
            list_1_rf.append(fm_1_rf)
            list_2_rf.append(fm_2_rf)
            list_fm_12_rf.append(fm_12_rf)
            
                               
            #controllo: momento nullo
            if abs(fm_12_rf.px) > 1e-9 or abs(fm_12_rf.py) > 1e-9 or abs(fm_12_rf.pz) > 1e-9 :
                print("Trimomento non nullo.")
        
        return list_1_rf, list_2_rf
        
#------------------------------------------------------------------------------------------------------------------------------

def compute_angle(v1, v2):
    """
    Calcola il coseno dell'angolo tra due vettori tridimensionali
    """
    p1 = np.array([v1.x, v1.y, v1.z])
    p2 = np.array([v2.x, v2.y, v2.z])

    cos_theta = np.dot(p1, p2) / (np.linalg.norm(p1) * np.linalg.norm(p2))
    return np.clip(cos_theta, -1.0, 1.0)

#------------------------------------------------------------------------------------------------------------------------------

def theta_star(quad_V, quad_H):
    '''
    Calcola l'angolo theta* per V con 
    - quad_V = lista di quadrimomenti VB; 
    - quad_H = lista di quadrimomenti H
    '''
    
    cos_list = []

    for fm_v, fm_h in zip(quad_V, quad_H):
        z_axis = vector.obj(x = 0, y = 0, z = 1)                                            #asse z nel laboratorio
        
        
        fm_VH = fm_v + fm_h
        
        boost_vec = vector.obj(x = -fm_VH.px / fm_VH.E, y = -fm_VH.py / fm_VH.E, z= -fm_VH.pz / fm_VH.E) 
        
        fm_v_rf = fm_v.boost(boost_vec)             #quadrimomento di V nel sdr VH
        fm_h_rf = fm_h.boost(boost_vec)
        
        p_vec = vector.obj(x = fm_v_rf.px, y = fm_v_rf.py, z = fm_v_rf.pz)                  #vettore momento della V nel SDR VH
        z_ax_rf = z_axis.boost(boost_vec)
        
        #angolo tra direzione di volo della V e la direzione del fascio
        cos_theta = compute_angle(p_vec, z_axis_rf)
        cos_list.append(cos_theta)

    return cos_list







'''
def theta_star(quad_vec_1, quad_vec_2) :

    Calcola l'angolo theta star per 1 con
    quad_vec_1 : lista di quadrimomenti della prima particella
    quad_vec_2 : lista di quadrimomenti della seconda particella

    
    list_fm_12_rf, theta_star_list = [], []
        
    for (fm_1, fm_2) in zip(quad_vec_1, quad_vec_2) :
        fm_12 = fm_1 + fm_2
        
        bx = fm_12.px / fm_12.E
        by = fm_12.py / fm_12.E
        bz = fm_12.pz / fm_12.E
        boost_vec = vector.obj(x=-bx, y=-by, z=-bz)
    
        fm_1_rf = fm_1.boost(boost_vec)
        fm_2_rf = fm_2.boost(boost_vec)
        fm_12_rf = fm_1_rf + fm_2_rf
        list_fm_12_rf.append(fm_12_rf)

        #calcolo theta star
        z_axis = vector.obj(x = 0, y = 0, z = 1)      
        p_vec = vector.obj(x=fm_1_rf.px, y=fm_1_rf.py, z=fm_1_rf.pz)
        theta_star = compute_angle(p_vec, z_axis)
        theta_star_list.append(theta_star)

                
        #controllo: momento nullo
        if abs(fm_12_rf.px) > 1e-9 or abs(fm_12_rf.py) > 1e-9 or abs(fm_12_rf.pz) > 1e-9 :
            print("Trimomento non nullo.")
    
    return theta_star_list
'''
