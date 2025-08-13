import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils import read_file, contains_particle, build_fm_ZW, boost_to_rf

def main () :
    '''
    Confronta le distribuzioni dell'angolo theta*, definito come la direzione di volo di Z (o W) calcolato nel sdr del centro di massa VH
    '''
    
    #calcolo del quadrimomento VH
    H_fm_L, H_fm_T = [], []
    Z_fm_L, Z_fm_T = [], []
    W_fm_L, W_fm_T = [], []    
    
    #longitudinale
    LHE_L = "unweighted_events_L.lhe"
    events_L = read_file(LHE_L)
    
    quadrivector_ZL, quadrivector_WL = build_fm_ZW(events_L)                                #liste di quadrimomenti per Z e W
    quadrivector_HL = []                                                                    #lista di quadrimomenti per l'higgs
    
    for event in events_L :
        for p in event :
            if p["pid"] == 25 :
                h_vec = vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])
                
                quadrivector_HL.append(h_vec)
                break
    
    #somma dei due
    fm_ZH_L = [x+y for x, y in zip(quadrivector_ZL, quadrivector_HL)]                       #quadrimomento ZH
    fm_WH_L = [x+y for x, y in zip(quadrivector_WL, quadrivector_HL)]                       #quadrimomento WH
    
    
    #trasversale
    LHE_T = "unweighted_events_T.lhe"
    events_T = read_file(LHE_T)

    quadrivector_ZT, quadrivector_WT = build_fm_ZW(events_T)                                #liste di quadrimomenti per Z e W
    quadrivector_HT = []                                                                    #lista di quadrimomenti per l'higgs

    for event in events_T :
        for p in event :
            if p["pid"] == 25 :
                h_vec = vector.obj(px = p["px"], py = p["py"], pz = p["pz"], E = p["E"])
                
                quadrivector_HT.append(h_vec)
                break

    #somma dei due
    fm_ZH_T = [x+y for x, y in zip(quadrivector_ZT, quadrivector_HT)]                       #quadrimomento ZH
    fm_WH_T = [x+y for x, y in zip(quadrivector_WT, quadrivector_HT)]                       #quadrimomento WH
    
        
    #applico boost di lorentz
    fm_ZH_L_rf = boost_to_rf(quadrivector_ZL, quadrivector_HL)
    fm_WH_L_rf = boost_to_rf(quadrivector_WL, quadrivector_HL)
    fm_ZH_T_rf = boost_to_rf(quadrivector_ZT, quadrivector_HT)
    fm_WH_T_rf = boost_to_rf(quadrivector_WT, quadrivector_HT)
    

    
if __name__ == '__main__' :
    main() 
