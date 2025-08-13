import vector

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
    
