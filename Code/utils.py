import vector
import numpy as np

def read_file(nome_file):

    with open(nome_file, 'r') as f:
            in_event = False
            event_lines = []
            events = []
            particles = []

            for line in f:
                line = line.strip()
                
                if "<event>" in line:
                    in_event = True
                    event_lines = []
                elif "</event>" in line:
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
                            E = float(parameter[9])
                            particles.append({
                                "pid": pid,
                                "status": status,
                                "mother1": mother1,
                                "mother2": mother2,
                                "px": px,
                                "py": py,
                                "pz": pz,
                                "E" : E
                            })
                    in_event = False
                    events.append(particles)
                    
                    particles = []
                else :
                    if in_event:
                        event_lines.append(line)
               
    return events
 
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

def theta_star(quad_V, quad_H):
    cos_list = []

    for fm_v, fm_h in zip(quad_V, quad_H):
        z_axis = vector.obj(x = 0, y = 0, z = 1)                                            #asse z nel laboratorio
        
        #quadrimomento totale VH
        fm_VH = fm_v + fm_h
        
        #calcolo beta come p/E di VH
        boost_vec = vector.obj(x = -fm_VH.px / fm_VH.E, y = -fm_VH.py / fm_VH.E, z= -fm_VH.pz / fm_VH.E) 
        
        #trovata beta applico -beta come boost
        fm_v_rf = fm_v.boost(boost_vec)             #quadrimomento di V nel sdr VH
        fm_h_rf = fm_h.boost(boost_vec)
        
        p_vec = vector.obj(x = fm_v_rf.px, y = fm_v_rf.py, z = fm_v_rf.pz)      #vettore momento della V nel SDR VH
        
        #angolo tra direzione di volo della V e la direzione del fascio
        cos_theta = compute_angle(p_vec, z_axis)
        cos_list.append(cos_theta)

    return cos_list


