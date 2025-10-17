import pickle
from utils_nonlhe import *



LHE_file = "/home/giorgia/MG/Risultati/H > bb/unweighted_events_L.lhe"
apply_smearing = True

w_longitudinale = build_df_dict_W(LHE_file, apply_smearing = apply_smearing)
z_longitudinale = build_df_dict_Z(LHE_file, apply_smearing = apply_smearing)



LHE_file = "/home/giorgia/MG/Risultati/H > bb/unweighted_events_T.lhe"
z_trasversale = build_df_dict_Z(LHE_file, apply_smearing = apply_smearing)
w_trasversale = build_df_dict_W(LHE_file, apply_smearing = apply_smearing)
   

with open("dataframes_z_long_selected.pkl", "wb") as f:
    pickle.dump(z_longitudinale, f)
    
with open("dataframes_w_long_selected.pkl", "wb") as f:
    pickle.dump(w_longitudinale, f)
    
with open("dataframes_z_tr_selected.pkl", "wb") as f:
    pickle.dump(z_trasversale, f)
    
with open("dataframes_w_tr_selected.pkl", "wb") as f:
    pickle.dump(w_trasversale, f)

