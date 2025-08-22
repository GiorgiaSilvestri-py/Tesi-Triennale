import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
import pandas as pd
import pickle
from scipy.stats import wasserstein_distance, ks_2samp
from utils import *
import sys

def build_df_Z(LHE_file):
    '''
    Costruisce un dataframe per la Z a partire dal file inserito come argomento
    '''

    pt_l1, pt_l2, eta_l1, eta_l2, phi_l1, phi_l2, pt_b_z = [], [], [], [], [], [], []
    pt_z, pt_h, eta_z, eta_h, phi_z, phi_h, M_zh = [], [], [], [], [], [], []
    cos_star, cos_1, phi_1, delta_phi = [], [], [], []
    
    #costruzione dizionari
    events = read_file(LHE_file)
    
    #isolo gli eventi della Z intanto
    fm_z_dict_all = get_fm_of([23], events)                                                    #costruisce il dizionario di quadrimomenti della Z
    lep_z_dict, antilep_z_dict, _, _ = connect_lep_to_V(events)                                 #dizionario di quadrimomenti di leptone e antileptone dalla Z
    fm_h_all_dict = get_fm_of([25], events)                                                    #dizionario con tutti gli eventi, da filtrare su Z
    
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
    
    for ((e_id1, fm_lep), (e_id2, fm_antilep)) in zip(lep_z_dict.items(), antilep_z_dict.items()):
        pt_l1.append(fm_lep.pt)
        pt_l2.append(fm_antilep.pt)
        
        eta_l1.append(fm_lep.eta)
        eta_l2.append(fm_antilep.eta)
        
        phi_l1.append(fm_lep.phi)
        phi_l2.append(fm_antilep.phi)
        
    pt_b_z, _ = get_lep_ptbalance(events)
    
    pt_z  = [fm_z.pt for fm_z in list(fm_z_dict.values())]
    pt_h  = [fm_h.pt for fm_h in list(fm_h_dict.values())]
    eta_z = [fm_z.eta for fm_z in list(fm_z_dict.values())]
    eta_h = [fm_h.eta for fm_h in list(fm_h_dict.values())]
    phi_z = [fm_z.phi for fm_z in list(fm_z_dict.values())]
    phi_h = [fm_h.phi for fm_h in list(fm_h_dict.values())]
    M_zh  = [fm_tot.M for fm_tot in list(compute_tot_fm(fm_z_dict, fm_h_dict).values())]
    
    cos_star  = get_thetastar_of(fm_z_dict, fm_h_dict)
    cos_1     = get_theta1_of(fm_z_dict, lep_z_dict, antilep_z_dict, fm_h_dict)
    phi_1     = get_phi1_of(fm_z_dict, lep_z_dict, antilep_z_dict, fm_h_dict)
    delta_phi = get_deltaphi_of(fm_z_dict, lep_z_dict, antilep_z_dict, fm_h_dict)

    #creazione dataframe
    log = {
            "Pt leptone"         : pt_l1,
            "Pt antileptone"     : pt_l2,
            "Eta leptone"        : eta_l1,
            "Eta antileptone"    : eta_l2,
            "Phi leptone"        : phi_l1,
            "Phi antileptone"    : phi_l2,
            "Pt Z"               : pt_z,
            "Pt H"               : pt_h,
            "Eta Z"              : eta_z,
            "Eta H"              : eta_h,
            "Phi Z"              : phi_z,
            "Phi H"              : phi_h,
            "Massa invariante ZH": M_zh,
            "Pt balance"         : pt_b_z,
            "cos(θ*)"            : cos_star,
            "cos(θ1)"            : cos_1,
            "φ1"                 : phi_1,
            "Δφ"                 : delta_phi
            }

    df = pd.DataFrame(log)
    
    print(df.shape)
    print(df.head())

    return(df)

def build_df_W(LHE_file):
    '''
    Costruisce un dataframe per la W a partire dal file inserito come argomento
    '''

    pt_l1, pt_l2, eta_l1, eta_l2, phi_l1, phi_l2, pt_b_w = [], [], [], [], [], [], []
    pt_w, pt_h, eta_w, eta_h, phi_w, phi_h, M_wh = [], [], [], [], [], [], []
    cos_star, cos_1, phi_1, delta_phi = [], [], [], []
    
    #costruzione dizionari
    events = read_file(LHE_file)
    
    #isolo gli eventi della w
    fm_w_dict = get_fm_of([24, -24], events)                                                    #costruisce il dizionario di quadrimomenti della Z
    _, _ , lep_w_dict, antilep_w_dict, = connect_alllep_to_V(events)                                 #dizionario di quadrimomenti di leptone e antileptone dalla Z
    fm_h_all_dict = get_fm_of([25], events)                                                    #dizionario con tutti gli eventi, da filtrare su Z
    
    #filtro H
    fm_h_dict = {}                                                                           #dizionario di fm H solo in presenza di z
    
    for event_id in fm_w_dict.keys() & fm_h_all_dict.keys():
        fm_w     = fm_w_dict[event_id]
        fm_higgs = fm_h_all_dict[event_id]        
        
        fm_h_dict[event_id] = fm_higgs
    
    
    #costruzione liste di variabili
    
    for ((e_id1, fm_lep), (e_id2, fm_antilep)) in zip(lep_w_dict.items(), antilep_w_dict.items()):
        pt_l1.append(fm_lep.pt)
        pt_l2.append(fm_antilep.pt)
        
        eta_l1.append(fm_lep.eta)
        eta_l2.append(fm_antilep.eta)
        
        phi_l1.append(fm_lep.phi)
        phi_l2.append(fm_antilep.phi)
        
    _, pt_b_w = get_lep_ptbalance(events)
    
    pt_w  = [fm_w.pt for fm_w in list(fm_w_dict.values())]
    pt_h  = [fm_h.pt for fm_h in list(fm_h_dict.values())]
    eta_w = [fm_w.eta for fm_w in list(fm_w_dict.values())]
    eta_h = [fm_h.eta for fm_h in list(fm_h_dict.values())]
    phi_w = [fm_w.phi for fm_w in list(fm_w_dict.values())]
    phi_h = [fm_h.phi for fm_h in list(fm_h_dict.values())]
    M_wh  = [fm_tot.M for fm_tot in list(compute_tot_fm(fm_w_dict, fm_h_dict).values())]
    
    cos_star  = get_thetastar_of(fm_w_dict, fm_h_dict)
    cos_1     = get_theta1_of(fm_w_dict, lep_w_dict, antilep_w_dict, fm_h_dict)
    phi_1     = get_phi1_of(fm_w_dict, lep_w_dict, antilep_w_dict, fm_h_dict)
    delta_phi = get_deltaphi_of(fm_w_dict, lep_w_dict, antilep_w_dict, fm_h_dict)

    #creazione dataframe
    log = {
            "Pt leptone"         : pt_l1,
            "Pt antileptone"     : pt_l2,
            "Eta leptone"        : eta_l1,
            "Eta antileptone"    : eta_l2,
            "Phi leptone"        : phi_l1,
            "Phi antileptone"    : phi_l2,
            "Pt W"               : pt_w,
            "Pt H"               : pt_h,
            "Eta W"              : eta_w,
            "Eta H"              : eta_h,
            "Phi W"              : phi_w,
            "Phi H"              : phi_h,
            "Massa invariante WH": M_wh,
            "Pt balance"         : pt_b_w,
            "cos(θ*)"            : cos_star,
            "cos(θ1)"            : cos_1,
            "φ1"                 : phi_1,
            "Δφ"                 : delta_phi
            }

    df = pd.DataFrame(log)
    
    print(df.shape)
    print(df.head())

    return(df)



def main():
    
    LHE_file = "/home/giorgia/MG/Risultati/Confronto_L-T/unweighted_events_L.lhe"
    df_long_z = build_df_Z(LHE_file)
    df_long_w = build_df_W(LHE_file)
    
    
    LHE_file = "/home/giorgia/MG/Risultati/Confronto_L-T/unweighted_events_T.lhe"
    df_tr_z = build_df_Z(LHE_file)
    df_tr_w = build_df_W(LHE_file)
    
    df_dict_z = {
              "longitudinale": df_long_z,
              "trasversale"  : df_tr_z
              }
    df_dict_w = {
              "longitudinale": df_long_w,
              "trasversale"  : df_tr_w
              }
    
    
    with open("dataframes_z.pkl", "wb") as f:
        pickle.dump(df_dict_z, f)
        
    with open("dataframes_w.pkl", "wb") as f:
        pickle.dump(df_dict_w, f)
    
if __name__ == '__main__':
    main()  
                    
