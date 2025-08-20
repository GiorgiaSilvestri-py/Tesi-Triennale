import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils import *

def get_theta1(events_dict):

	cos_list_z, cos_list_W = [], []

	#costruisco i dizionari che mi servono
	w_dict = get_fm_of([24, -24], events_dict)                                 #W
	z_dict = get_fm_of([23], events_dict)                                      #Z    
	h_dict = get_fm_of([25], events_dict)                                      #H
	lep_z, antilep_z, lep_w, antilep_w = connect_lep_to_V(events_dict)     #dizionari di leptoni e antileptoni da Z/W
	n_dict = get_fm_of([12, 14, 16], events_dict)                              #neutrini
	antin_dict = get_fm_of([-12, -14, -16], events_dict)                       #antineutrini

	#procedimento per la z
	for event_id in z_dict.keys() & lep_z.keys() & antilep_z.keys() & h_dict.keys():
		fm_z = z_dict[event_id]
		fm_lep = lep_z[event_id]
		fm_antilep = antilep_z[event_id]
		fm_h = h_dict[event_id]
	
		#calcolo del vettore di boost
		boost_vec_z = vector.obj(x = -fm_z.px/fm_z.E, y = -fm_z.py/fm_z.E, z = -fm_z.pz/fm_z.E)
	
		#boost di leptone, H e antileptone nel SDR della Z
		fm_lep_sdrz = fm_lep.boost(boost_vec_z)
		fm_antilep_sdrz = fm_antilep.boost(boost_vec_z)
		fm_h_sdrz = fm_h.boost(boost_vec_z)
		
		#controllo boost
		fm_tot_z = fm_lep_sdrz + fm_antilep_sdrz
		if abs(fm_tot_z.px) > 1e-9 or abs(fm_tot_z.py) > 1e-9 or abs(fm_tot_z.pz) > 1e-9 :
			print("Evento N°:", event_id, ". Trimomento non nullo.")
	
		#identifica le direzioni e calcola l'angolo theta 1
		lep_direction_z = vector.obj(x = fm_lep_sdrz.px, y = fm_lep_sdrz.py, z = fm_lep_sdrz.pz)
		h_direction_z = vector.obj(x = fm_h_sdrz.px, y = fm_h_sdrz.py, z = fm_h_sdrz.pz)
		
		theta_1_z = compute_angle(lep_direction_z, h_direction_z)
		cos_list_z.append(theta_1_z)

	#procedimento per la W
	for event_id in w_dict.keys() & lep_w.keys() & h_dict.keys():
		fm_w = w_dict[event_id]
		fm_lep = lep_w[event_id]
		fm_h = h_dict[event_id]

		#calcolo del vettore di boost
		boost_vec_w = vector.obj(x = -fm_w.px/fm_w.E, y = -fm_w.py/fm_w.E, z = -fm_w.pz/fm_w.E)
	
		#boost di leptone, H nel SDR della w
		fm_lep_sdrw = fm_lep.boost(boost_vec_w)
		fm_h_sdrw = fm_h.boost(boost_vec_w)
		
		#identifica le direzioni e calcola l'angolo theta 1
		lep_direction_w = vector.obj(x = fm_lep_sdrw.px, y = fm_lep_sdrw.py, z = fm_lep_sdrw.pz)
		h_direction_w = vector.obj(x = fm_h_sdrw.px, y = fm_h_sdrw.py, z = fm_h_sdrw.pz)
		
		theta_1_w = compute_angle(lep_direction_w, h_direction_w)
		cos_list_w.append(theta_1_w)
		
	return cos_list_z, cos_list_w


def main () :
    '''
    Confronta le distribuzioni dell'angolo theta 1, definito come l'angolo tra il leptone e la direzione di volo dell'Higgs
    '''
    '''
    Procedimento:
    - SDR comovente con la Z
	    - trovo il quadrimomento della Z da quello dei suoi prodotti di decadimento
	    - calcolo beta come px/E, py/E etc con px, py, pz momento della Z
	    - fatto ciò, applico il boost (con -beta) al leptone
	    - ho ottenuto ora il quadrimomento dei due leptoni e il quadrimomento somma nel SDR della Z
    - applico il boost anche al bosone di Higgs -> ottengo il quadrimomento di h nel SDR della Z
    - a questo punto, calcolo l'angolo, tramite compute_angle, con p1 = p_leptone e p2 = p_Higgs
    '''
    
    #longitudinale
    LHE_L = "unweighted_events_L.lhe"
    eventsL_dict = read_file(LHE_L)
    
    #costruzione quadrimomenti
    W_fm_L = get_fm_of([24, -24], eventsL_dict)                                 #W
    Z_fm_L = get_fm_of([23], eventsL_dict)                                      #Z    
    H_fm_L = get_fm_of([25], eventsL_dict)                                      #H
    lep_ZL, antilep_ZL, _, _ = connect_lep_to_V(eventsL_dict)                    #dizionari di leptoni e antileptoni da Z/W
    l_lep_w, l_al_w = connect_alllep_to_V(eventsL_dict)[2:]
    
    #trasversale
    LHE_T = "unweighted_events_T.lhe"
    eventsT_dict = read_file(LHE_T)
    
    W_fm_T = get_fm_of([24, -24], eventsT_dict)                                 #W
    Z_fm_T = get_fm_of([23], eventsT_dict)                                      #Z
    H_fm_T = get_fm_of([25], eventsT_dict)                                      #H
    lep_ZT,  antilep_ZT, _, _ = connect_lep_to_V(eventsT_dict)                  #lep/antilep
    t_lep_w, t_al_w = connect_alllep_to_V(eventsT_dict)[2:]
    
    #liste di cos(theta)
    Z_cos_list_L = get_theta1_of(Z_fm_L, lep_ZL, antilep_ZL, H_fm_L)
    Z_cos_list_T = get_theta1_of(Z_fm_T, lep_ZT, antilep_ZT, H_fm_T)
    W_cos_list_L = get_theta1_of(W_fm_L, l_lep_w, l_al_w, H_fm_L)
    W_cos_list_T = get_theta1_of(W_fm_T, t_lep_w, t_al_w, H_fm_T)
    
    
    #visualizzazione distribuzioni

    data = Z_cos_list_L
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)
    
    data = W_cos_list_L
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nbins = int((max(data) - min(data)) / bin_width)


    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(Z_cos_list_L, bins=35, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[0],  label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(Z_cos_list_T, bins=35, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='Polarizzazione trasversale')
    ax[0].set_xlabel("cosθ1 (rad)")
    ax[0].set_ylabel("dN/N")
    ax[0].set_title("Distribuzione angolo θ1 (Z)")
    ax[0].legend()
    
    sb.histplot(W_cos_list_L, bins=nbins, color='royalblue', edgecolor = 'steelblue', stat='density', ax=ax[1], label='Polarizzazione longitudinale',  alpha = 0.8)
    sb.histplot(W_cos_list_T, bins=nbins, color='firebrick', stat='density', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='Polarizzazione trasversale')
    ax[1].set_xlabel("cosθ1 (rad)")
    ax[1].set_ylabel("dN/N")
    ax[1].set_title("Distribuzione angolo θ1 (W)")
    ax[1].legend()

    #plt.savefig("Confronto coseno.png", dpi=300)
    plt.show()

  
if __name__ == '__main__' :
    main() 
