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
            if pid in p_id :
                fm_dict[event_id] = cinematic.build_fm()
    '''        
    for event_id, fm in fm_dict.items() :
        if event_id > 10:
            break
        print("Event ID: ", event_id, "\n\t Quadrimomento: ", fm)
    '''
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
            if abs(pid) == p_id :
                pt_dict[event_id] = cinematic.compute_pt()
    '''        
    for event_id, pt in pt_dict.items() :
        if event_id > 10:
            break
        print("Event ID: ", event_id, "\n\t Momento Trasverso: ", pt)
    '''   
    return pt_dict
    
def get_eta_of(p_id, events_dict):
    '''
    Crea la lista di eta di una particella (p_id)
    - scorre gli eventi
    - per ciascun evento cerca se c'è la particella
    - se esiste, calcola eta
    '''
    eta_dict = {}
    
    for event_id, particles_dict in events_dict.items():
        for pid, cinematic in particles_dict.items():
            if abs(pid) == p_id :
                eta_dict[event_id] = cinematic.compute_eta()
    '''        
    for event_id, pt in pt_dict.items() :
        if event_id > 10:
            break
        print("Event ID: ", event_id, "\n\t Momento Trasverso: ", pt)
    '''   
    return eta_dict
    
    
def get_phi_of(p_id, events_dict):
    '''
    Crea la lista di eta di una particella (p_id)
    - scorre gli eventi
    - per ciascun evento cerca se c'è la particella
    - se esiste, calcola eta
    '''
    phi_dict = {}
    
    for event_id, particles_dict in events_dict.items():
        for pid, cinematic in particles_dict.items():
            if abs(pid) == p_id :
                phi_dict[event_id] = cinematic.compute_phi()
    '''        
    for event_id, pt in pt_dict.items() :
        if event_id > 10:
            break
        print("Event ID: ", event_id, "\n\t Momento Trasverso: ", pt)
    '''   
    return phi_dict
    
    
def get_thetastar_of(V_fm_dict, H_fm_dict):
    '''
    Calcola l'angolo theta* per V con 
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
    
    fm_VH = compute_tot_fm(V_fm_dict, H_fm_dict)                #dizionario di quadrimomenti totali

    for event_id in V_fm_dict.keys() & H_fm_dict.keys() & fm_VH.keys():
        fm_v = V_fm_dict[event_id]
        fm_h = H_fm_dict[event_id]
        fm_tot = fm_VH[event_id]

        #costruzione vettore di boost del singolo quadrimomento fm_tot
        boost_vec = vector.obj(x = -fm_tot.px/fm_tot.E, y = -fm_tot.py/fm_tot.E, z = -fm_tot.pz/fm_tot.E)  
    
        #boost di V e H nel rest frame VH
        fm_v_rf = fm_v.boost(boost_vec)
        fm_h_rf = fm_h.boost(boost_vec)
        
        #controllo del boost
        fm_tot_rf = fm_v_rf + fm_h_rf
        if abs(fm_tot_rf.px) > 1e-5 or abs(fm_tot_rf.py) > 1e-5 or abs(fm_tot_rf.pz) > 1e-5 :
            print("Evento N°:", event_id, ". Trimomento non nullo: ", fm_tot_rf)
        
        #calcolo delle direzioni e dell'angolo
        z_axis = vector.obj(px = 0, py = 0, pz = 1, E=1)                                #asse z
        z_axis_rf = z_axis.boost(boost_vec)
        p_vec = vector.obj(x = fm_v_rf.px, y = fm_v_rf.py, z = fm_v_rf.pz)      #direzione di volo di V
        
        cos_theta = compute_angle(p_vec, z_axis_rf)
        cos_list.append(cos_theta)
        
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
        fm_tot = fm_VH[event_id]

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
    Calcola l'angolo theta fra la direzione di volo del leptone (dal decadimento della V) e quella di H, nel SDR della V.
    Procedimento:
    - impone identità su id_evento
    - calcola il vettore di boost dal quadrimomento del VB
    - boosta leptone e higgs nel sdr del VB
    - controlla che il sdr sia corretto (il momento somma deve essere nullo)
    - calcola l'angolo
    '''

    cos_list = []

    for event_id in V_dict.keys() & lep_dict.keys() & antilep_dict.keys() & H_dict.keys():
        fm_v = V_dict[event_id]
        fm_lep = lep_dict[event_id]
        fm_antilep = antilep_dict[event_id]
        fm_h = H_dict[event_id]
        
        boost_vec = vector.obj(x = -fm_v.px/fm_v.E, y = -fm_v.py/fm_v.E, z = -fm_v.pz/fm_v.E)
        
        fm_lep_rf = fm_lep.boost(boost_vec)
        fm_antilep_rf = fm_antilep.boost(boost_vec)
        fm_h_rf = fm_h.boost(boost_vec)
        
        #controllo boost
        fm_tot = fm_lep_rf + fm_antilep_rf
        if abs(fm_tot.px) > 1e-5 or abs(fm_tot.py) > 1e-5 or abs(fm_tot.pz) > 1e-5 :
            print("Evento N°:", event_id, ". Trimomento non nullo: ", fm_tot)
        
        lep_direction = vector.obj(x = fm_lep_rf.px, y = fm_lep_rf.py, z = fm_lep_rf.pz)
        h_direction = vector.obj(x = fm_h_rf.px, y = fm_h_rf.py, z = fm_h_rf.pz)
        
        theta_1 = compute_angle(lep_direction, h_direction)
        cos_list.append(theta_1)
        
    return cos_list
    
def phi1(V_dict, lep_dict, antilep_dict, H_dict):
    '''
    Calcola l'angolo φ1 tra il piano del decadimento del bosone e il piano di scattering partonico nel sdr VH.
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
        
        boost_vec = vector.obj(x = -fm_tot_i.px/fm_tot_i.E, y = -fm_tot_i.py/fm_tot_i.E, z = -fm_tot_i.pz/fm_tot_i.E)
        
        fm_lep_rf = fm_lep.boost(boost_vec)
        fm_v_rf = fm_v.boost(boost_vec)
        fm_antilep_rf = fm_antilep.boost(boost_vec)
        fm_h_rf = fm_h.boost(boost_vec)
        z_ax_r = z_axis.boost(boost_vec)
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


def compute_tot_fm(fm1_dict, fm2_dict):
    '''
    Calcola il quadrimomento somma
    '''
    fm_tot_dict = {}
    
    for event_id in fm1_dict.keys() & fm2_dict.keys():
        fm_1 = fm1_dict[event_id]
        fm_2 = fm2_dict[event_id]

        fm_tot_dict[event_id] = fm_1 + fm_2
    '''       
    for event_id, fm_tot in fm_tot_dict.items() :
        if event_id > 10:
            break
        print("Event ID: ", event_id, "\n\t Quadrimomento Somma: ", fm_tot)
    '''   
    return fm_tot_dict
    
    
def connect_lep_to_V(events_dict):
    '''
    Ricostruisce le liste di leptoni decaduti dalla Z o dalla W
    '''
    
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
    '''
    Ricostruisce le liste di leptoni decaduti dalla Z o dalla W (inclusi neutrini)
    '''
    
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
