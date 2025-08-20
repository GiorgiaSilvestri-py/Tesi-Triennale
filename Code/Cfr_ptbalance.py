import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils import *

def get_lep_ptbalance(events_dict):
    '''
    calcola il lepton pt balance
    '''
    
    pt_balance_list_z, pt_balance_list_w = [], []
    
    #dizionari necessari
    z_fm = get_fm_of([23], events_dict)
    zlep, zantilep, wlep, wantilep = connect_lep_to_V(events_dict)
    
    for event_id in z_fm.keys() & zlep.keys() & zantilep.keys() :
        p4z = z_fm[event_id]
        p4lep = zlep[event_id]
        p4antilep = zantilep[event_id]
        
        ptlep = p4lep.pt
        ptantilep = p4antilep.pt
        
        pt_sum = ptlep + ptantilep
        
        if pt_sum == 0:
            continue
            
        pt_balance = abs(ptlep - ptantilep) / pt_sum
    
        pt_balance_list_z.append(pt_balance)
        
    w_fm = get_fm_of([24, -24], events_dict)
    zlep, zantilep, wlep, wantilep = connect_alllep_to_V(events_dict)
    
    for event_id in w_fm.keys() & wlep.keys() & wantilep.keys():
        p4w = w_fm[event_id]
        p4lep = wlep[event_id]
        p4antilep = wantilep[event_id]
        
        ptlep = p4lep.pt
        ptantilep = p4antilep.pt
        
        pt_sum = ptlep + ptantilep
        
        if pt_sum == 0:
            continue
            
        pt_balance = abs(ptlep - ptantilep) / pt_sum
    
        pt_balance_list_w.append(pt_balance)

    return pt_balance_list_z, pt_balance_list_w
    
            
def main():
    '''
    Confronto tramite pt Z e W + grafico scatter
    '''
    
    #analisi polarizzazione longitudinale
    LHE_L = "unweighted_events_L.lhe"
    eventsL_dict = read_file(LHE_L)
    
    #costruzione quadrimomenti
    W_fm_L = get_fm_of([24, -24], eventsL_dict)                                 #W
    Z_fm_L = get_fm_of([23], eventsL_dict)                                      #Z    
    H_fm_L = get_fm_of([25], eventsL_dict)                                      #H
    zlep_l, zantilep_l, wlep_l, wantilep_l = connect_lep_to_V(eventsL_dict)     #dizionari di leptoni e antileptoni da Z/W
    
    
    
    #analisi polarizzazione trasversale
    LHE_T = "unweighted_events_T.lhe"
    eventsT_dict = read_file(LHE_T)
    
    #quadrimomenti
    W_fm_T = get_fm_of([24, -24], eventsT_dict)               
    Z_fm_T = get_fm_of([23], eventsT_dict)               
    zlep_t, zantilep_t, wlep_t, wantilep_t = connect_lep_to_V(eventsT_dict)
          
    
    ptb_l, ptb_lw = get_lep_ptbalance(eventsL_dict)
    ptb_t, ptb_tw = get_lep_ptbalance(eventsT_dict)
       

    data = ptb_l
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)

    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(ptb_l, bins=40, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(ptb_t, bins=40, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversale')
    ax[0].set_xlabel("Pt balance")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione del lepton pt balance (Z)")
    ax[0].legend()
    
    sb.histplot(ptb_lw, bins=40, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(ptb_tw, bins=40, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversale')
    ax[1].set_xlabel("Pt balance")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione del lepton pt balance (W)")
    ax[1].legend()

    plt.savefig("Confronto lepton pt balance.png", dpi=300)
    plt.show()
 
        
if __name__ == '__main__':
    main()  
                    
