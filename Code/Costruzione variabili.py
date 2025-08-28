import pandas as pd
import pickle
from utils import *


def main():
    
    LHE_file = "/home/giorgia/MG/Risultati/Confronto_L-T/unweighted_events_L.lhe"
    z_longitudinale = build_df_dict_Z(LHE_file)
    w_longitudinale = build_df_dict_W(LHE_file)


    LHE_file = "/home/giorgia/MG/Risultati/Confronto_L-T/unweighted_events_T.lhe"
    z_trasversale = build_df_dict_Z(LHE_file)
    w_trasversale = build_df_dict_W(LHE_file)
       
    
    with open("dataframes_z_long.pkl", "wb") as f:
        pickle.dump(z_longitudinale, f)
        
    with open("dataframes_w_long.pkl", "wb") as f:
        pickle.dump(w_longitudinale, f)
        
    with open("dataframes_z_tr.pkl", "wb") as f:
        pickle.dump(z_trasversale, f)
        
    with open("dataframes_w_tr.pkl", "wb") as f:
        pickle.dump(w_trasversale, f)
    
if __name__ == '__main__':
    main()  
                    
