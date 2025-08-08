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
    
