import pandas as pd
import pickle
from tabulate import tabulate
import vector
from utils_nonlhe import *

    
def main():

    print("Longitudinal polarisation")

    LHE_file = "unweighted_events_L.lhe"

    events_dict = read_file(LHE_file)


    total_ZL, e_mu_ZL = [], []
    total_WL, e_mu_WL = [], []
    selected_ZL, selected_WL = [], []
    cut1_ZL, cut2_ZL = [], []
    cut1_WL, cut2_WL, cut3_WL, cut12_WL = [], [], [], []
    cut31_WL, cut31_WT = [], []

    for e_id, particle_list in events_dict.items():

        if e_id % 10000 == 0:
            print(f"Processing event {e_id}")

       
        event_type = f_event_type(particle_list)

        #di-lepton
        if event_type == 23:

            total_ZL.append(e_id)

            e_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            mu_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            gamma1_fm = vector.obj(x = 0, y = 0, z = 0, E = 0)
            gamma2_fm = vector.obj(x = 0, y = 0, z = 0, E = 0)
            

            for cinematic in particle_list:
                #leptons
                if abs(cinematic.pid) == 11:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                    e_vec = smeared_momentum

                if  abs(cinematic.pid) == 13:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                    mu_vec = smeared_momentum
                    

            for cinematic in particle_list: 
                    
                if cinematic.pid == 22:
                    non_smeared_momentum1 = cinematic.build_fm()
                    non_smeared_momentum1_bis = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum1, "photon")
                    gamma1_fm = smeared_momentum
                    break

            for cinematic in particle_list: 
                
                if cinematic.pid == 22:
                    non_smeared_momentum2 = cinematic.build_fm()
                    if vector_is_equal(non_smeared_momentum1_bis, non_smeared_momentum2, 0.1):
                        continue
                    smeared_momentum = apply_smearing_to(non_smeared_momentum2, "photon")
                    gamma2_fm = smeared_momentum


            #selections
            #leading / subleading photons pt
            sum_fm = gamma1_fm + gamma2_fm
            m_gg = sum_fm.M
            pt_gamma1 = np.sqrt(gamma1_fm.px**2 + gamma1_fm.py**2)
            if pt_gamma1 == 0 :
                print("Null pt1 in event", e_id)
            
            pt_gamma2 = np.sqrt(gamma2_fm.px**2 + gamma2_fm.py**2)
            if pt_gamma2 == 0 :
                print("Null pt1 in event", e_id)
            
            #print("pt gamma2", pt_gamma2)

            leading_threshold = 3*m_gg / 8
            #print("leading threshold:", leading_threshold)
            subleading_threshold = m_gg / 4
            #print("subleading threshold:", subleading_threshold)

            #lepton
            
            if is_nonzero(e_vec):
          
                e_mu_ZL.append(e_id)

                deltaphi_e_gamma1 = gamma1_fm.deltaphi(e_vec)
                deltaphi_e_gamma2 = gamma2_fm.deltaphi(e_vec)

                deltaeta_e_gamma1 = abs(gamma1_fm.eta - e_vec.eta)
                deltaeta_e_gamma2 = abs(gamma2_fm.eta - e_vec.eta)

                R_gamma1_e = np.sqrt( (deltaeta_e_gamma1)**2 + (deltaphi_e_gamma1)**2 )
                R_gamma2_e = np.sqrt( (deltaeta_e_gamma2)**2 + (deltaphi_e_gamma2)**2 )

                #defining cuts
                cut2 = R_gamma1_e > 1.0 and R_gamma2_e > 1.0
                
                if pt_gamma1 > pt_gamma2:
                    cut1 = pt_gamma1 > leading_threshold and pt_gamma2 > subleading_threshold
                    
                    if cut1:
                        cut1_ZL.append(e_id)

                    if cut2:
                        cut2_ZL.append(e_id)

                    if cut1 and cut2:
                        selected_ZL.append(e_id)
                    
                else:
                    cut1 = pt_gamma2 > leading_threshold and pt_gamma1 > subleading_threshold

                    #if pt_gamma2 > leading_threshold and pt_gamma1 > subleading_threshold:
                    if cut1:
                        cut1_ZL.append(e_id)

                    if cut2:
                        cut2_ZL.append(e_id)
                        
                    if cut1 and cut2:
                        #print(pt_gamma1, pt_gamma2, R_gamma1_e, R_gamma2_e)
                        selected_ZL.append(e_id)
                    


            if is_nonzero(mu_vec):
     

                e_mu_ZL.append(e_id)

                deltaphi_mu_gamma1 = gamma1_fm.deltaphi(mu_vec)
                deltaphi_mu_gamma2 = gamma2_fm.deltaphi(mu_vec)

                deltaeta_mu_gamma1 = abs(gamma1_fm.eta - mu_vec.eta)
                deltaeta_mu_gamma2 = abs(gamma2_fm.eta - mu_vec.eta)

                R_gamma1_mu = np.sqrt( (deltaeta_mu_gamma1)**2 + (deltaphi_mu_gamma1)**2 )
                R_gamma2_mu = np.sqrt( (deltaeta_mu_gamma2)**2 + (deltaphi_mu_gamma2)**2 )

                #defining cuts
                cut2 = R_gamma1_mu > 0.5 and R_gamma2_mu > 0.5

                if pt_gamma1 > pt_gamma2:
                    cut1 = pt_gamma1 > leading_threshold and pt_gamma2 > subleading_threshold

                    if cut1:
                        cut1_ZL.append(e_id)

                    if cut2:
                        cut2_ZL.append(e_id)
                    
                    if cut1 and cut2:
                        selected_ZL.append(e_id)
                    

                else:
                    cut1 = pt_gamma2 > leading_threshold and pt_gamma1 > subleading_threshold

                    if cut1:
                        cut1_ZL.append(e_id)
                    if cut2:
                        cut2_ZL.append(e_id)
                    if cut1 and cut2:
                        selected_ZL.append(e_id)



        #single-lepton
        if event_type == 24 or event_type == -24:

            total_WL.append(e_id)

            lepton_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            neutrino_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            gamma1_fm = vector.obj(x = 0, y = 0, z = 0, E = 0)
            gamma2_fm = vector.obj(x = 0, y = 0, z = 0, E = 0)
            
                        
            #leptons & met
            for cinematic in particle_list:
                if abs(cinematic.pid) in [11, 13]:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                    lepton_vec = smeared_momentum

            for cinematic in particle_list:
                if abs(cinematic.pid) in [12, 14]:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "neutrino")
                    neutrino_vec = compute_neutrino_momentum_from_particles(lepton_vec, smeared_momentum)

            if is_zero(lepton_vec) or is_zero(neutrino_vec):
                continue
                

            e_mu_WL.append(e_id)

            #photons
            for cinematic in particle_list: 
                    
                if cinematic.pid == 22:
                    non_smeared_momentum1 = cinematic.build_fm()
                    non_smeared_momentum1_bis = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum1, "photon")
                    gamma1_fm = smeared_momentum
                    #print("gamma 1 vec:", gamma1_fm)
                    break

            for cinematic in particle_list: 
                
                if cinematic.pid == 22:
                    non_smeared_momentum2 = cinematic.build_fm()
                    
                    if vector_is_equal(non_smeared_momentum1_bis, non_smeared_momentum2, 0.1):
                        
                        continue

                    smeared_momentum = apply_smearing_to(non_smeared_momentum2, "photon")
                    gamma2_fm = smeared_momentum
                    #print("gamma 2 vec:", gamma2_fm)

                
                    
            #selections
            sum_fm = gamma1_fm + gamma2_fm

            m_gg = sum_fm.M

            pt_gamma1 = np.sqrt(gamma1_fm.px**2 + gamma1_fm.py**2)
            #print("pt gamma 1 vec:", pt_gamma1)
            pt_gamma2 = np.sqrt(gamma2_fm.px**2 + gamma2_fm.py**2)
            #print("pt gamma 2 vec:", pt_gamma2)

            leading_threshold = 3*m_gg / 8
            #print("leading threshold:", leading_threshold))
            subleading_threshold = m_gg / 4
            #print("subleading threshold:", subleading_threshold)

            pt_MET = np.sqrt(neutrino_vec.px**2 + neutrino_vec.py**2)

            deltaphi_lep_gamma1 = gamma1_fm.deltaphi(lepton_vec)
            deltaphi_lep_gamma2 = gamma2_fm.deltaphi(lepton_vec)

            deltaeta_lep_gamma1 = abs(gamma1_fm.eta - lepton_vec.eta)
            deltaeta_lep_gamma2 = abs(gamma2_fm.eta - lepton_vec.eta)

            R_gamma1_lep = np.sqrt( (deltaeta_lep_gamma1)**2 + (deltaphi_lep_gamma1)**2 )
            R_gamma2_lep = np.sqrt( (deltaeta_lep_gamma2)**2 + (deltaphi_lep_gamma2)**2 )

            #defining cuts
            cut2 = R_gamma1_lep > 1.0 and R_gamma2_lep > 1.0
            cut3 = pt_MET > 45

            if pt_gamma1 > pt_gamma2:
                cut1 = pt_gamma1 > leading_threshold and pt_gamma2 > subleading_threshold
                
                if cut1:
                    cut1_WL.append(e_id)
                if cut2:
                    cut2_WL.append(e_id)
                if cut3:
                    cut3_WL.append(e_id)

                if cut1 and cut2:
                    cut12_WL.append(e_id)

                if cut3 and cut1:
                    cut31_WL.append(e_id)

                if cut1 and cut2 and cut3:
                    #print(pt_gamma1, pt_gamma2, R_gamma1_lep, R_gamma2_lep, pt_MET)
                    selected_WL.append(e_id)

            else:
                cut1 = pt_gamma2 > leading_threshold and pt_gamma1 > subleading_threshold

                if cut1:
                    cut1_WL.append(e_id)

                if cut2:
                    cut2_WL.append(e_id)
                
                if cut3:
                    cut3_WL.append(e_id)


                if cut1 and cut2:
                    cut12_WL.append(e_id)

                if cut3 and cut1:
                    cut31_WL.append(e_id)
                
                if cut1 and cut2 and cut3:
                    selected_WL.append(e_id)
            

    #----------------------------------------------------------------------------------      


    print("\nTransverse polarisation \n")

    LHE_file = "unweighted_events_T.lhe"
    events_dict = read_file(LHE_file)

    total_ZT, e_mu_ZT = [], []
    total_WT, e_mu_WT = [], []
    selected_ZT, selected_WT = [], []
    cut1_ZT, cut2_ZT = [], []
    cut1_WT, cut2_WT, cut3_WT, cut12_WT = [], [], [], []
    cut31_WT = []

    for e_id, particle_list in events_dict.items():

        if e_id % 10000 == 0:
            print(f"Processing event {e_id}")

        event_type = f_event_type(particle_list)

        #di-lepton
        if event_type == 23:

            total_ZT.append(e_id)

            e_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            mu_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            gamma1_fm = vector.obj(x = 0, y = 0, z = 0, E = 0)
            gamma2_fm = vector.obj(x = 0, y = 0, z = 0, E = 0)
            

            for cinematic in particle_list:
                #leptons
                if abs(cinematic.pid) == 11:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                    e_vec = smeared_momentum

                if  abs(cinematic.pid) == 13:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                    mu_vec = smeared_momentum
                    

            for cinematic in particle_list: 
                    
                if cinematic.pid == 22:
                    non_smeared_momentum1 = cinematic.build_fm()
                    non_smeared_momentum1_bis = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum1, "photon")
                    gamma1_fm = smeared_momentum
                    break

            for cinematic in particle_list: 
                
                if cinematic.pid == 22:
                    non_smeared_momentum2 = cinematic.build_fm()
                    if vector_is_equal(non_smeared_momentum1_bis, non_smeared_momentum2, 0.1):
                        continue
                    smeared_momentum = apply_smearing_to(non_smeared_momentum2, "photon")
                    gamma2_fm = smeared_momentum


            #selections
            #leading / subleading photons pt
            sum_fm = gamma1_fm + gamma2_fm
            m_gg = sum_fm.M
            pt_gamma1 = np.sqrt(gamma1_fm.px**2 + gamma1_fm.py**2)
            if pt_gamma1 == 0 :
                print("Null pt1 in event", e_id)
            #print("pt gamma1", pt_gamma1)
            pt_gamma2 = np.sqrt(gamma2_fm.px**2 + gamma2_fm.py**2)
            if pt_gamma2 == 0 :
                print("Null pt1 in event", e_id)
            
            #print("pt gamma2", pt_gamma2)

            leading_threshold = 3*m_gg / 8
            #print("leading threshold:", leading_threshold)
            subleading_threshold = m_gg / 4
            #print("subleading threshold:", subleading_threshold)

            #lepton
            
            if is_nonzero(e_vec):
                #print(e_id, e_vec)
                
        
                e_mu_ZT.append(e_id)

                deltaphi_e_gamma1 = gamma1_fm.deltaphi(e_vec)
                deltaphi_e_gamma2 = gamma2_fm.deltaphi(e_vec)

                deltaeta_e_gamma1 = abs(gamma1_fm.eta - e_vec.eta)
                deltaeta_e_gamma2 = abs(gamma2_fm.eta - e_vec.eta)

                R_gamma1_e = np.sqrt( (deltaeta_e_gamma1)**2 + (deltaphi_e_gamma1)**2 )
                R_gamma2_e = np.sqrt( (deltaeta_e_gamma2)**2 + (deltaphi_e_gamma2)**2 )

                #defining cuts
                cut2 = R_gamma1_e > 1.0 and R_gamma2_e > 1.0
                
                if pt_gamma1 > pt_gamma2:
                    cut1 = pt_gamma1 > leading_threshold and pt_gamma2 > subleading_threshold
                    
                    #if pt_gamma1 > leading_threshold and pt_gamma2 > subleading_threshold:
                    if cut1:
                        cut1_ZT.append(e_id)


                    #if R_gamma1_e > 1.0 and R_gamma2_e > 1.0:
                    if cut2:
                        cut2_ZT.append(e_id)

                    if cut1 and cut2:
                    
                        selected_ZT.append(e_id)
                    
                else:
                    cut1 = pt_gamma2 > leading_threshold and pt_gamma1 > subleading_threshold

                    #if pt_gamma2 > leading_threshold and pt_gamma1 > subleading_threshold:
                    if cut1:
                        cut1_ZT.append(e_id)

                    if cut2:
                        cut2_ZT.append(e_id)
                        
                    if cut1 and cut2:
                        #print(pt_gamma1, pt_gamma2, R_gamma1_e, R_gamma2_e)
                        selected_ZT.append(e_id)
                    


            if is_nonzero(mu_vec):


                e_mu_ZT.append(e_id)

                deltaphi_mu_gamma1 = gamma1_fm.deltaphi(mu_vec)
                deltaphi_mu_gamma2 = gamma2_fm.deltaphi(mu_vec)

                deltaeta_mu_gamma1 = abs(gamma1_fm.eta - mu_vec.eta)
                deltaeta_mu_gamma2 = abs(gamma2_fm.eta - mu_vec.eta)

                R_gamma1_mu = np.sqrt( (deltaeta_mu_gamma1)**2 + (deltaphi_mu_gamma1)**2 )
                R_gamma2_mu = np.sqrt( (deltaeta_mu_gamma2)**2 + (deltaphi_mu_gamma2)**2 )

                #defining cuts
                cut2 = R_gamma1_mu > 0.5 and R_gamma2_mu > 0.5

                if pt_gamma1 > pt_gamma2:
                    cut1 = pt_gamma1 > leading_threshold and pt_gamma2 > subleading_threshold

                    if cut1:
                        cut1_ZT.append(e_id)

                    if cut2:
                        cut2_ZT.append(e_id)
                    
                    if cut1 and cut2:
                        selected_ZT.append(e_id)
                    

                else:
                    cut1 = pt_gamma2 > leading_threshold and pt_gamma1 > subleading_threshold

                    if cut1:
                        cut1_ZT.append(e_id)
                    if cut2:
                        cut2_ZT.append(e_id)
                    if cut1 and cut2:
                        selected_ZT.append(e_id)


        #single-lepton
        if event_type == 24 or event_type == -24:

            total_WT.append(e_id)

            lepton_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            neutrino_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            gamma1_fm = vector.obj(x = 0, y = 0, z = 0, E = 0)
            gamma2_fm = vector.obj(x = 0, y = 0, z = 0, E = 0)
            
                        
            #leptons & met
            for cinematic in particle_list:
                if abs(cinematic.pid) in [11, 13]:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                    lepton_vec = smeared_momentum

            for cinematic in particle_list:
                if abs(cinematic.pid) in [12, 14]:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "neutrino")
                    neutrino_vec = compute_neutrino_momentum_from_particles(lepton_vec, smeared_momentum)

            if is_zero(lepton_vec) or is_zero(neutrino_vec):
                continue
                

            e_mu_WT.append(e_id)

            #photons
            for cinematic in particle_list: 
                    
                if cinematic.pid == 22:
                    non_smeared_momentum1 = cinematic.build_fm()
                    non_smeared_momentum1_bis = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum1, "photon")
                    gamma1_fm = smeared_momentum
                    #print("gamma 1 vec:", gamma1_fm)
                    break

            for cinematic in particle_list: 
                
                if cinematic.pid == 22:
                    non_smeared_momentum2 = cinematic.build_fm()
                    
                    if vector_is_equal(non_smeared_momentum1_bis, non_smeared_momentum2, 0.1):
                        
                        continue

                    smeared_momentum = apply_smearing_to(non_smeared_momentum2, "photon")
                    gamma2_fm = smeared_momentum
                    #print("gamma 2 vec:", gamma2_fm)

                
                    
            #selections
            sum_fm = gamma1_fm + gamma2_fm

            m_gg = sum_fm.M

            pt_gamma1 = np.sqrt(gamma1_fm.px**2 + gamma1_fm.py**2)
            #print("pt gamma 1 vec:", pt_gamma1)
            pt_gamma2 = np.sqrt(gamma2_fm.px**2 + gamma2_fm.py**2)
            #print("pt gamma 2 vec:", pt_gamma2)

            leading_threshold = 3*m_gg / 8
            #print("leading threshold:", leading_threshold))
            subleading_threshold = m_gg / 4
            #print("subleading threshold:", subleading_threshold)

            pt_MET = np.sqrt(neutrino_vec.px**2 + neutrino_vec.py**2)

            deltaphi_lep_gamma1 = gamma1_fm.deltaphi(lepton_vec)
            deltaphi_lep_gamma2 = gamma2_fm.deltaphi(lepton_vec)

            deltaeta_lep_gamma1 = abs(gamma1_fm.eta - lepton_vec.eta)
            deltaeta_lep_gamma2 = abs(gamma2_fm.eta - lepton_vec.eta)

            R_gamma1_lep = np.sqrt( (deltaeta_lep_gamma1)**2 + (deltaphi_lep_gamma1)**2 )
            R_gamma2_lep = np.sqrt( (deltaeta_lep_gamma2)**2 + (deltaphi_lep_gamma2)**2 )

            #defining cuts
            cut2 = R_gamma1_lep > 1.0 and R_gamma2_lep > 1.0
            cut3 = pt_MET > 45

            if pt_gamma1 > pt_gamma2:
                cut1 = pt_gamma1 > leading_threshold and pt_gamma2 > subleading_threshold
                
                if cut1:
                    cut1_WT.append(e_id)
                if cut2:
                    cut2_WT.append(e_id)
                if cut3:
                    cut3_WT.append(e_id)

                if cut1 and cut2:
                    cut12_WT.append(e_id)

                if cut3 and cut1:
                    cut31_WT.append(e_id)

                if cut1 and cut2 and cut3:
                    #print(pt_gamma1, pt_gamma2, R_gamma1_lep, R_gamma2_lep, pt_MET)
                    selected_WT.append(e_id)

            else:
                cut1 = pt_gamma2 > leading_threshold and pt_gamma1 > subleading_threshold

                if cut1:
                    cut1_WT.append(e_id)
                if cut2:
                    cut2_WT.append(e_id)
                if cut3:
                    cut3_WT.append(e_id)

                if cut1 and cut2:
                    cut12_WT.append(e_id)

                if cut3 and cut1:
                    cut31_WT.append(e_id)
                
                if cut1 and cut2 and cut3:
                    selected_WT.append(e_id)
            
            
    #--------------------------------------------------------------------------------------------

    #ZH dataframe
    N_tot_ZL      = len(total_ZL)
    N_em_ZL       = len(e_mu_ZL)
    N_cut1_ZL     = len(cut1_ZL)
    N_cut2_ZL     = len(cut2_ZL)
    N_cut12_ZL    = len(selected_ZL)

    N_tot_ZT      = len(total_ZT)
    N_em_ZT       = len(e_mu_ZT)
    N_cut1_ZT     = len(cut1_ZT)
    N_cut2_ZT     = len(cut2_ZT)
    N_cut12_ZT    = len(selected_ZT)

    #wh dataframe
    N_tot_WL      = len(total_WL)
    N_em_WL       = len(e_mu_WL)
    N_cut1_WL     = len(cut1_WL)
    N_cut2_WL     = len(cut2_WL)
    N_cut12_WL    = len(cut12_WL)
    N_cut3_WL     = len(cut3_WL)
    N_cut31_WL    = len(cut31_WL)
    N_cut123_WL   = len(selected_WL)

    N_tot_WT      = len(total_WT)
    N_em_WT       = len(e_mu_WT)
    N_cut1_WT     = len(cut1_WT)
    N_cut2_WT     = len(cut2_WT)
    N_cut12_WT    = len(cut12_WT)
    N_cut3_WT     = len(cut3_WT)
    N_cut31_WT    = len(cut31_WT)
    N_cut123_WT   = len(selected_WT)


    #channel dataframe
    data = [
            ["Total events",            N_tot_ZL,   N_tot_ZT,   N_tot_WL,    N_tot_WT],
            ["Events with e/mu",        N_em_ZL,    N_em_ZT,    N_em_WL,     N_em_WT],
            ["Events after cut1 only",  N_cut1_ZL,  N_cut1_ZT,  N_cut1_WL,   N_cut1_WT],
            ["Events after cut2 only",  N_cut2_ZL,  N_cut2_ZT,  N_cut2_WL,   N_cut2_WT],
            ["Events after cut3 only",  "-",        "-",        N_cut3_WL,   N_cut3_WT],
            ["Eventi selezionati",      N_cut12_ZL, N_cut12_ZT, N_cut123_WL, N_cut123_WT]
            ]


    df_channel = pd.DataFrame(data, columns=["", "ZLH", "ZTH", "WLH", "WTH"])

    print("\n" + "="*60)
    print("{:^60}".format("H > γγ"))
    print("="*60 + "\n")
    print(tabulate(df_channel.values.tolist(), headers=df_channel.columns.tolist(), tablefmt="grid"))


    #-------------------------------------------------------------------------------------------------------


    print("Computing efficiencies")
    #ZLH
    #total efficiencies
    epsilon_ZL = N_tot_ZL / 50000
    epsilon_0_ZL = N_em_ZL / 50000

    epsilon_1_ZL = N_cut1_ZL / N_em_ZL
    epsilon_2_ZL = N_cut2_ZL / N_em_ZL
    epsilon_tot_ZL = N_cut12_ZL / N_tot_ZL

    #partial efficiencies
    part_epsilon_0_ZL = N_em_ZL / N_tot_ZL
    part_epsilon_1_ZL = N_cut1_ZL / N_em_ZL
    part_epsilon_2_ZL = N_cut12_ZL / N_cut1_ZL

    #ZTH
    #total efficiencies
    epsilon_ZT = N_tot_ZT / 50000
    epsilon_0_ZT = N_em_ZT / 50000

    epsilon_1_ZT = N_cut1_ZT / N_em_ZT
    epsilon_2_ZT = N_cut2_ZT / N_em_ZT
    epsilon_tot_ZT = N_cut12_ZT / N_tot_ZT

    #partial efficiencies
    part_epsilon_0_ZT = N_em_ZT / N_tot_ZT
    part_epsilon_1_ZT = N_cut1_ZT / N_em_ZT
    part_epsilon_2_ZT = N_cut12_ZT / N_cut1_ZT

    #WLH
    #total efficiencies
    epsilon_WL = N_tot_WL / 50000
    epsilon_0_WL = N_em_WL / 50000

    epsilon_1_WL = N_cut1_WL / N_em_WL
    epsilon_2_WL = N_cut2_WL / N_em_WL
    epsilon_3_WL = N_cut3_WL / N_em_WL
    epsilon_tot_WL = N_cut123_WL / N_tot_WL

    #partial efficiencies
    part_epsilon_0_WL = N_em_WL / N_tot_WL
    part_epsilon_1_WL = N_cut1_WL / N_em_WL
    part_epsilon_2_WL = N_cut12_WL / N_cut1_WL
    part_epsilon_3_WL = N_cut123_WL / N_cut12_WL


    #WTH
    #total efficiencies
    epsilon_WT = N_tot_WT / 50000
    epsilon_0_WT = N_em_WT / 50000

    epsilon_1_WT = N_cut1_WT / N_em_WT
    epsilon_2_WT = N_cut2_WT / N_em_WT
    epsilon_3_WT = N_cut3_WT / N_em_WT
    epsilon_tot_WT = N_cut123_WT / N_tot_WT

    #partial efficiencies
    part_epsilon_0_WT = N_em_WT / N_tot_WT
    part_epsilon_1_WT = N_cut1_WT / N_em_WT
    part_epsilon_2_WT = N_cut12_WT / N_cut1_WT
    part_epsilon_3_WT = N_cut123_WT / N_cut12_WT
    

    efficiencies = np.array([
                            ["channel fraction",        f"{epsilon_ZL:.2f}",     f"{epsilon_ZT:.2f}",     f"{epsilon_WL:.2f}",     f"{epsilon_WT:.2f}"],
                            ["ε₀ = e/μ / total",        f"{epsilon_0_ZL:.2f}",   f"{epsilon_0_ZT:.2f}",   f"{epsilon_0_WL:.2f}",   f"{epsilon_0_WT:.2f}"],
                            ["ε₁ = cut1 / e/μ",         f"{epsilon_1_ZL:.2f}",   f"{epsilon_1_ZT:.2f}",   f"{epsilon_1_WL:.2f}",   f"{epsilon_1_WT:.2f}"],
                            ["ε₂ = cut2 / e/μ",         f"{epsilon_2_ZL:.2f}",   f"{epsilon_2_ZT:.2f}",   f"{epsilon_2_WL:.2f}",   f"{epsilon_2_WT:.2f}"],
                            ["ε₃ = cut3 / e/μ",         "-",                     "-",                     f"{epsilon_3_WL:.2f}",   f"{epsilon_3_WT:.2f}"],
                            ["ε_tot = selected events", f"{epsilon_tot_ZL:.2f}", f"{epsilon_tot_ZT:.2f}", f"{epsilon_tot_WL:.2f}", f"{epsilon_tot_WT:.2f}"]
                        ])
    
    column = ["", "ZLH", "ZTH", "WLH", "WTH"]
    df_efficiencies = pd.DataFrame(efficiencies, columns=column)

    print("\n" + "="*60)
    print("{:^60}".format("Total efficiencies"))
    print("="*60 + "\n")
    print(tabulate(df_efficiencies.values.tolist(), headers=df_efficiencies.columns.tolist(), tablefmt="grid"))


    part_total_eff = np.array([
                                ["ε₀ (e/μ)",      f"{part_epsilon_0_ZL:.2f}", f"{epsilon_0_ZL:.2f}", f"{part_epsilon_0_ZT:.2f}", f"{epsilon_0_ZT:.2f}", f"{part_epsilon_0_WL:.2f}", f"{epsilon_0_WL:.2f}", f"{part_epsilon_0_WT:.2f}", f"{epsilon_0_WT:.2f}"],
                                ["ε₁ (cut1)",     f"{part_epsilon_1_ZL:.2f}", f"{epsilon_1_ZL:.2f}", f"{part_epsilon_1_ZT:.2f}", f"{epsilon_1_ZT:.2f}", f"{part_epsilon_1_WL:.2f}", f"{epsilon_1_WL:.2f}", f"{part_epsilon_1_WT:.2f}", f"{epsilon_1_WT:.2f}"],
                                ["ε₂ (cut2)",     f"{part_epsilon_2_ZL:.2f}", f"{epsilon_2_ZL:.2f}", f"{part_epsilon_2_ZT:.2f}", f"{epsilon_2_ZT:.2f}", f"{part_epsilon_2_WL:.2f}", f"{epsilon_2_WL:.2f}", f"{part_epsilon_2_WT:.2f}", f"{epsilon_2_WT:.2f}"],
                                ["ε₃ (cut3)",     "-",                        "-",                   "-",                        "-",                   f"{part_epsilon_3_WL:.2f}", f"{epsilon_3_WL:.2f}", f"{part_epsilon_3_WT:.2f}", f"{epsilon_3_WT:.2f}"],
                            ])
    
    columns = ["", 
               "ZLH cumulative eff", "ZLH total eff", 
               "ZTH cumulative eff", "ZTH total eff",
               "WLH cumulative eff", "WLH total eff",
               "WTH cumulative eff", "WTH total eff"
               ]

    df_f = pd.DataFrame(part_total_eff, columns = columns)

    print("\n" + "="*120)
    print("{:^120}".format("Cumulative and Total efficiencies"))
    print("="*120 + "\n")

    print(tabulate(df_f.values.tolist(), headers=df_f.columns.tolist(), tablefmt="grid"))


    #-------------------------------------------------------------------------------------------------------


    #computing LHC expected events number

    xsection_L = 282    #fb
    BR         = 0.00270
    xsection_T = 252  #fb
    lum        = 100            #1/fb 

    #ZLH
    N_tot_L_LHC     = xsection_L    * lum                                      #total longitudinal events number
    N_tot_ZL_LHC    = N_tot_L_LHC   * epsilon_ZL                              #ZLH events number
    N_tot_ZL_diph   = N_tot_ZL_LHC  * BR
    N_em_ZL_LHC     = N_tot_ZL_diph * part_epsilon_0_ZL
    N_cut1_ZL_LHC   = N_em_ZL_LHC   * part_epsilon_1_ZL
    N_cut2_ZL_LHC   = N_cut1_ZL_LHC * part_epsilon_2_ZL

    kin_eff_ZL_LHC  = N_cut2_ZL_LHC / N_em_ZL_LHC                           #kinematic cuts efficiency
    chan_eff_ZL_LHC = N_cut2_ZL_LHC / N_tot_ZL_diph                           #channel efficiency
    tot_eff_ZL_LHC  = N_cut2_ZL_LHC / N_tot_ZL_LHC


   # ZTH
    N_tot_T_LHC     = xsection_T    * lum
    N_tot_ZT_LHC    = N_tot_T_LHC   * epsilon_ZT
    N_tot_ZT_diph   = N_tot_ZT_LHC  * BR
    N_em_ZT_LHC     = N_tot_ZT_diph * part_epsilon_0_ZT
    N_cut1_ZT_LHC   = N_em_ZT_LHC   * part_epsilon_1_ZT
    N_cut2_ZT_LHC   = N_cut1_ZT_LHC * part_epsilon_2_ZT

    kin_eff_ZT_LHC  = N_cut2_ZT_LHC / N_em_ZT_LHC
    chan_eff_ZT_LHC = N_cut2_ZT_LHC / N_tot_ZT_diph                           #channel efficiency
    tot_eff_ZT_LHC  = N_cut2_ZT_LHC / N_tot_ZT_LHC

    # WLH
    N_tot_L_LHC     = xsection_L    * lum
    N_tot_WL_LHC    = N_tot_L_LHC   * epsilon_WL
    N_tot_WL_diph   = N_tot_WL_LHC  * BR
    N_em_WL_LHC     = N_tot_WL_diph * part_epsilon_0_WL
    N_cut1_WL_LHC   = N_em_WL_LHC   * part_epsilon_1_WL
    N_cut2_WL_LHC   = N_cut1_WL_LHC * part_epsilon_2_WL
    N_cut3_WL_LHC   = N_cut2_WL_LHC * part_epsilon_3_WL

    kin_eff_WL_LHC  = N_cut2_WL_LHC / N_em_WL_LHC
    chan_eff_WL_LHC = N_cut2_WL_LHC / N_tot_WL_diph                           #channel efficiency
    tot_eff_WL_LHC  = N_cut2_WL_LHC / N_tot_WL_LHC

    # WTH
    N_tot_T_LHC     = xsection_T    * lum
    N_tot_WT_LHC    = N_tot_T_LHC   * epsilon_WT
    N_tot_WT_diph   = N_tot_WT_LHC  * BR
    N_em_WT_LHC     = N_tot_WT_diph * part_epsilon_0_WT
    N_cut1_WT_LHC   = N_em_WT_LHC   * part_epsilon_1_WT
    N_cut2_WT_LHC   = N_cut1_WT_LHC * part_epsilon_2_WT
    N_cut3_WT_LHC   = N_cut2_WT_LHC * part_epsilon_3_WT

    kin_eff_WT_LHC  = N_cut2_WT_LHC / N_em_WT_LHC
    chan_eff_WT_LHC = N_cut2_WT_LHC / N_tot_WT_diph                           #channel efficiency
    tot_eff_WT_LHC  = N_cut2_WT_LHC / N_tot_WT_LHC


    data = [
            ["Total V{i}H events",            f"{N_tot_L_LHC:.0f}",    f"{N_tot_T_LHC:.0f}",    f"{N_tot_L_LHC:.0f}",    f"{N_tot_T_LHC:.0f}"],
            ["Total ch events",               f"{N_tot_ZL_LHC:.0f}",   f"{N_tot_ZT_LHC:.0f}",    f"{N_tot_WL_LHC:.0f}",    f"{N_tot_WT_LHC:.0f}"],
            ["Total γγ events",               f"{N_tot_ZL_diph:.0f}",  f"{N_tot_ZT_diph:.0f}",    f"{N_tot_WL_diph:.0f}",    f"{N_tot_WT_diph:.0f}"],
            ["Events with e/mu",              f"{N_em_ZL_LHC:.0f}",    f"{N_em_ZT_LHC:.0f}",     f"{N_em_WL_LHC:.0f}",     f"{N_em_WT_LHC:.0f}"],
            ["Events after cut1",             f"{N_cut1_ZL_LHC:.0f}",  f"{N_cut1_ZT_LHC:.0f}",   f"{N_cut1_WL_LHC:.0f}",   f"{N_cut1_WT_LHC:.0f}"],
            ["Events after cut1, cut2",       f"{N_cut2_ZL_LHC:.0f}",  f"{N_cut2_ZT_LHC:.0f}",   f"{N_cut2_WL_LHC:.0f}",   f"{N_cut2_WT_LHC:.0f}"],
            ["Events after cut1, cut2, cut3", "-",                      "-",                      f"{N_cut3_WL_LHC:.0f}",   f"{N_cut3_WT_LHC:.0f}"],
            ["Kinematic cuts efficiency",     f"{kin_eff_ZL_LHC:.2f}", f"{kin_eff_ZT_LHC:.2f}",  f"{kin_eff_WL_LHC:.2f}",  f"{kin_eff_WT_LHC:.2f}"],
            ["Channel efficiency",            f"{chan_eff_ZL_LHC:.2f}",f"{chan_eff_ZT_LHC:.2f}",   f"{chan_eff_WL_LHC:.2f}",   f"{chan_eff_WT_LHC:.2f}"],
            ["Total efficiency",              f"{tot_eff_ZL_LHC:.5f}",  f"{tot_eff_ZT_LHC:.5f}",  f"{tot_eff_WL_LHC:.4f}",  f"{tot_eff_WT_LHC:.4f}"]
            ]


    df_channel = pd.DataFrame(data, columns=["", "ZLH", "ZTH", "WLH", "WTH"])

    print("\n" + "="*60)
    print("{:^60}".format("H > γγ (LHC)"))
    print("="*60 + "\n")
    print(tabulate(df_channel.values.tolist(), headers=df_channel.columns.tolist(), tablefmt="grid"))


    #----------------------------------------------------------------------------------------------------------------------------------------------------------------


    print("Ordering cuts")

    cut_list = ["cut1", "cut2", "cut3"]

    for col in column[1:]:
        
        #re-arranging
        values = df_efficiencies.loc[1:4, col].copy()
        num_values = [float(v) for v in values if v != "-"]

        ordered_val = sorted(num_values)
        print("Ordered total efficiencies:", ordered_val)
        sorted_indices = sorted(range(len(num_values)), key=lambda i: num_values[i])

        ordered_cut_list = [cut_list[i] for i in sorted_indices]
        print(f"Ordered cuts for {col}:", ordered_cut_list)

        j = 0
        for m in range(1, 5): 

            if df_efficiencies.at[m, col] == "-":
                continue

            df_efficiencies.at[m, col] = f"{ordered_val[j]:.2f}"
            j += 1

    #ordering cumulative efficiencies
    #ZLH

    #ZTH

    N_cut31_WL = len(cut31_WL)
    N_cut31_WT = len(cut31_WT)

    #WLH
    #partial efficiencies
    part_epsilon_1_WL = N_cut3_WL / N_em_WL
    part_epsilon_2_WL = N_cut31_WL / N_cut3_WL
    part_epsilon_3_WL = N_cut123_WL / N_cut31_WL

    #WTH
    #partial efficiencies
    part_epsilon_1_WT = N_cut3_WT / N_em_WT
    part_epsilon_2_WT = N_cut31_WT / N_cut3_WT
    print(f"part_epsilon_2_WT = {N_cut31_WT} / {N_cut3_WT}")
    part_epsilon_3_WT = N_cut123_WT / N_cut31_WT
    

    part_total_eff_ordered = np.array([
                                ["ε₀ = e/μ / total",    f"{epsilon_0_ZL:.2f}",      f"{epsilon_0_ZL:.2f}", f"{epsilon_0_ZT:.2f}",      f"{epsilon_0_ZT:.2f}", f"{epsilon_0_WL:.2f}",      f"{epsilon_0_WL:.2f}", f"{epsilon_0_WT:.2f}",      f"{epsilon_0_WT:.2f}"],
                                ["ε₁ (cut1)",           f"{part_epsilon_1_ZL:.2f}", f"{epsilon_1_ZL:.2f}", f"{part_epsilon_1_ZT:.2f}", f"{epsilon_1_ZT:.2f}", f"{part_epsilon_1_WL:.2f}", f"{epsilon_3_WL:.2f}", f"{part_epsilon_1_WT:.2f}", f"{epsilon_3_WT:.2f}"],
                                ["ε₂ (cut2)",           f"{part_epsilon_2_ZL:.2f}", f"{epsilon_2_ZL:.2f}", f"{part_epsilon_2_ZT:.2f}", f"{epsilon_2_ZT:.2f}", f"{part_epsilon_2_WL:.2f}", f"{epsilon_1_WL:.2f}", f"{part_epsilon_2_WT:.2f}", f"{epsilon_1_WT:.2f}"],
                                ["ε₃ (cut3)",           "-",                        "-",                   "-",                        "-",                   f"{part_epsilon_3_WL:.2f}", f"{epsilon_2_WL:.2f}", f"{part_epsilon_3_WT:.2f}", f"{epsilon_2_WT:.2f}"],
                            ])
    
    columns = ["", 
               "ZLH cumulative eff.", "ZLH total eff", 
               "ZTH cumulative eff.", "ZTH total eff",
               "WLH cumulative eff", "WLH total eff",
               "WTH cumulative eff", "WTH total eff"
               ]

    df_ordered = pd.DataFrame(part_total_eff_ordered, columns = columns)

    print("\n" + "="*120)
    print("{:^120}".format("Cumulative and Total efficiencies - [ordered by increasing total efficiency]"))
    print("="*120 + "\n")
    print(tabulate(df_ordered.values.tolist(), headers=df_ordered.columns.tolist(), tablefmt="grid"))


    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    
    

if __name__ == '__main__':
    main()