import pickle
from utils_nonlhe import *


#longitudinal
LHE_file = "H_ττ/unweighted_events_L.lhe"
apply_smearing = True
all_events = read_file(LHE_file)

id_events = apply_selections(all_events)
events_dict = {}

for id in id_events:
    events_dict[id] = all_events[id]

with open("selected_events_L.pkl", "wb") as f:
    pickle.dump(events_dict, f)




#transverse
LHE_file = "H_ττ/unweighted_events_T.lhe"

apply_smearing = True
all_events = read_file(LHE_file)

id_events = apply_selections(all_events)
events_dict = {}

for id in id_events:
    events_dict[id] = all_events[id]

with open("selected_events_T.pkl", "wb") as f:
    pickle.dump(events_dict, f)
