from array import array
import os
import pandas as pd
import pickle
from tabulate import tabulate
from utils_nonlhe import *
import vector

def compute_deltaR(vec1, vec2):

    deltaphi = vec1.deltaphi(vec2)
    deltaeta = abs(vec1.eta - vec2.eta)

    return np.sqrt( (deltaphi)**2 + (deltaeta)**2 )

def main():

    print("Longitudinal polarisation")

    LHE_file = "unweighted_events_L.lhe"
    events_dict = read_file(LHE_file)

    tree_ZLH = ROOT.TTree("Events", "ZLH channel")
    pt_leading_ZLH    = array('f', [0.])
    pt_subleading_ZLH = array('f', [0.])
    leading_th_ZLH    = array('f', [0.])
    subleading_th_ZLH = array('f', [0.])
    R1_ZLH            = array('f', [0.])
    R2_ZLH            = array('f', [0.])
    event_type_ZLH    = ROOT.std.string()

    tree_ZLH.Branch("pt_leading", pt_leading_ZLH, "pt_leading/F")
    tree_ZLH.Branch("pt_subleading", pt_subleading_ZLH, "pt_subleading/F")
    tree_ZLH.Branch("leading_th", leading_th_ZLH, "leading_th/F")
    tree_ZLH.Branch("subleading_th", subleading_th_ZLH, "subleading_th/F")
    tree_ZLH.Branch("R1", R1_ZLH, "R1/F")
    tree_ZLH.Branch("R2", R2_ZLH, "R2/F")
    tree_ZLH.Branch("event_type", event_type_ZLH)

    # TTree WLH
    tree_WLH = ROOT.TTree("Events", "WLH channel")

    pt_leading_WLH    = array('f', [0.])
    pt_subleading_WLH = array('f', [0.])
    leading_th_WLH    = array('f', [0.])
    subleading_th_WLH = array('f', [0.])
    R1_WLH            = array('f', [0.])
    R2_WLH            = array('f', [0.])
    pt_MET_WLH        = array('f', [0.])
    event_type_WLH    = ROOT.std.string()

    tree_WLH.Branch("pt_leading", pt_leading_WLH, "pt_leading/F")
    tree_WLH.Branch("pt_subleading", pt_subleading_WLH, "pt_subleading/F")
    tree_WLH.Branch("leading_th", leading_th_WLH, "leading_th/F")
    tree_WLH.Branch("subleading_th", subleading_th_WLH, "subleading_th/F")
    tree_WLH.Branch("R1", R1_WLH, "R1/F")
    tree_WLH.Branch("R2", R2_WLH, "R2/F")
    tree_WLH.Branch("pt_MET", pt_MET_WLH, "pt_MET/F")
    tree_WLH.Branch("event_type", event_type_WLH)

    #--------------------------------------------------------------------------

    #events loop

    for e_id, particle_list in events_dict.items():

        if e_id % 10000 == 0:
            print(f"Processing event {e_id}")

        #if e_id > 1000:    break

        event_type = f_event_type(particle_list)

        #di-lepton
        if event_type == 23:

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
            m_gg   = sum_fm.M
            
            if gamma1_fm.pt == 0 :                      #check
                print("Null pt1 in event", e_id)
            
            if gamma2_fm.pt == 0 :
                print("Null pt1 in event", e_id)
            
            if is_nonzero(e_vec):
                lep_vec = e_vec
                event_type_ZLH.replace(0, len(event_type_ZLH), "e")

            elif is_nonzero(mu_vec):
                lep_vec = mu_vec
                event_type_ZLH.replace(0, len(event_type_ZLH), "mu")
            else:
                event_type_ZLH.replace(0, len(event_type_ZLH), "neutrino/tau")
        
            R_gamma1_lep = compute_deltaR(gamma1_fm, lep_vec)
            R_gamma2_lep = compute_deltaR(gamma2_fm, lep_vec)
        
            pt_leading_ZLH[0]    = max(gamma1_fm.pt, gamma2_fm.pt)
            pt_subleading_ZLH[0] = min(gamma1_fm.pt, gamma2_fm.pt)
            leading_th_ZLH[0]    = 3*m_gg / 8
            subleading_th_ZLH[0] = m_gg / 4
            R1_ZLH[0]            = R_gamma1_lep
            R2_ZLH[0]            = R_gamma2_lep

            tree_ZLH.Fill()
             

        #single-lepton
        if event_type == 24 or event_type == -24:

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

            #photons
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
            sum_fm = gamma1_fm + gamma2_fm
            m_gg   = sum_fm.M
            pt_MET = np.sqrt(neutrino_vec.px**2 + neutrino_vec.py**2)

            if gamma1_fm.pt == 0 :                      #check
                print("Null pt1 in event", e_id)
            
            if gamma2_fm.pt == 0 :
                print("Null pt1 in event", e_id)

            #lepton            
            if is_nonzero(lepton_vec):
                lep_vec = lepton_vec
                event_type_WLH.replace(0, len(event_type_WLH), "e/mu")
            else:
                event_type_WLH.replace(0, len(event_type_WLH), "tau")
          
            R_gamma1_lep = compute_deltaR(gamma1_fm, lep_vec)
            R_gamma2_lep = compute_deltaR(gamma2_fm, lep_vec)
        
            pt_leading_WLH[0]    = max(gamma1_fm.pt, gamma2_fm.pt)
            pt_subleading_WLH[0] = min(gamma1_fm.pt, gamma2_fm.pt)
            leading_th_WLH[0]    = 3*m_gg / 8
            subleading_th_WLH[0] = m_gg / 4
            R1_WLH[0]            = R_gamma1_lep
            R2_WLH[0]            = R_gamma2_lep
            pt_MET_WLH[0]        = pt_MET

            tree_WLH.Fill()


    #----------------------------------------------------------------------

    print(f"Displaying VLH output example dataframe\n")

    df_ZLH = ROOT.RDataFrame(tree_ZLH)
    df_WLH = ROOT.RDataFrame(tree_WLH)
    
    print("Example ZLH dataframe")
    df_ZLH.Display().Print()
    
    print("Example WLH dataframe")
    df_WLH.Display().Print()

    #----------------------------------------------------------------------------------      


    print("\nTransverse polarisation \n")

    LHE_file = "unweighted_events_T.lhe"
    events_dict = read_file(LHE_file)

    #TTree initialisation - ZTH
    tree_ZTH = ROOT.TTree("Events", "ZTH channel")

    pt_leading_ZTH    = array('f', [0.])
    pt_subleading_ZTH = array('f', [0.])
    leading_th_ZTH    = array('f', [0.])
    subleading_th_ZTH = array('f', [0.])
    R1_ZTH            = array('f', [0.])
    R2_ZTH            = array('f', [0.])
    event_type_ZTH    = ROOT.std.string()

    tree_ZTH.Branch("pt_leading",    pt_leading_ZTH, "pt_leading/F")
    tree_ZTH.Branch("pt_subleading", pt_subleading_ZTH, "pt_subleading/F")
    tree_ZTH.Branch("leading_th",    leading_th_ZTH, "leading_th/F")
    tree_ZTH.Branch("subleading_th", subleading_th_ZTH, "subleading_th/F")
    tree_ZTH.Branch("R1",            R1_ZTH, "R1/F")
    tree_ZTH.Branch("R2",            R2_ZTH, "R2/F")
    tree_ZTH.Branch("event_type",    event_type_ZTH)

    #TTree initialisation - WTH
    tree_WTH = ROOT.TTree("Events", "WTH channel")

    pt_leading_WTH    = array('f', [0.])
    pt_subleading_WTH = array('f', [0.])
    leading_th_WTH    = array('f', [0.])
    subleading_th_WTH = array('f', [0.])
    R1_WTH            = array('f', [0.])
    R2_WTH            = array('f', [0.])
    pt_MET_WTH        = array('f', [0.])
    event_type_WTH    = ROOT.std.string()

    tree_WTH.Branch("pt_leading",    pt_leading_WTH, "pt_leading/F")
    tree_WTH.Branch("pt_subleading", pt_subleading_WTH, "pt_subleading/F")
    tree_WTH.Branch("leading_th",    leading_th_WTH, "leading_th/F")
    tree_WTH.Branch("subleading_th", subleading_th_WTH, "subleading_th/F")
    tree_WTH.Branch("R1",            R1_WTH, "R1/F")
    tree_WTH.Branch("R2",            R2_WTH, "R2/F")
    tree_WTH.Branch("pt_MET",        pt_MET_WTH, "pt_MET/F")
    tree_WTH.Branch("event_type",    event_type_WTH)

    #--------------------------------------------------------------------------------

    #events loop
    for e_id, particle_list in events_dict.items():

        if e_id % 10000 == 0:
            print(f"Processing event {e_id}")

        event_type = f_event_type(particle_list)

        #di-lepton
        if event_type == 23:

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
                    
            #photons
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
            sum_fm = gamma1_fm + gamma2_fm
            m_gg   = sum_fm.M

            if gamma1_fm.pt == 0 :                      #check
                print("Null pt1 in event", e_id)
            
            if gamma2_fm.pt == 0 :
                print("Null pt1 in event", e_id)
            
            if is_nonzero(e_vec):
                lep_vec = e_vec
                event_type_ZTH.replace(0, len(event_type_ZTH), "e")
            elif is_nonzero(mu_vec):
                lep_vec = mu_vec
                event_type_ZTH.replace(0, len(event_type_ZTH), "mu")
            else:
                event_type_ZTH.replace(0, len(event_type_ZTH), "neutrino/tau")
        
            R_gamma1_lep = compute_deltaR(gamma1_fm, lep_vec)
            R_gamma2_lep = compute_deltaR(gamma2_fm, lep_vec)
        
            pt_leading_ZTH[0]    = max(gamma1_fm.pt, gamma2_fm.pt)
            pt_subleading_ZTH[0] = min(gamma1_fm.pt, gamma2_fm.pt)
            leading_th_ZTH[0]    = 3*m_gg / 8
            subleading_th_ZTH[0] = m_gg / 4
            R1_ZTH[0]            = R_gamma1_lep
            R2_ZTH[0]            = R_gamma2_lep

            tree_ZTH.Fill()

            
        #single-lepton
        if event_type == 24 or event_type == -24:

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

            #photons
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
            sum_fm = gamma1_fm + gamma2_fm
            m_gg   = sum_fm.M
            pt_MET = np.sqrt(neutrino_vec.px**2 + neutrino_vec.py**2)

            if gamma1_fm.pt == 0 :                      #check
                print("Null pt1 in event", e_id)
            
            if gamma2_fm.pt == 0 :
                print("Null pt1 in event", e_id)
            
            #lepton            
            if is_nonzero(lepton_vec):
                lep_vec = lepton_vec
                event_type_WTH.replace(0, len(event_type_WTH), "e/mu")
            else:
                event_type_WTH.replace(0, len(event_type_WTH), "tau")
          
            R_gamma1_lep = compute_deltaR(gamma1_fm, lep_vec)
            R_gamma2_lep = compute_deltaR(gamma2_fm, lep_vec)
        
            pt_leading_WTH[0]    = max(gamma1_fm.pt, gamma2_fm.pt)
            pt_subleading_WTH[0] = min(gamma1_fm.pt, gamma2_fm.pt)
            leading_th_WTH[0]    = 3*m_gg / 8
            subleading_th_WTH[0] = m_gg / 4
            R1_WTH[0]            = R_gamma1_lep
            R2_WTH[0]            = R_gamma2_lep
            pt_MET_WTH[0]        = pt_MET

            tree_WTH.Fill()

    #-------------------------------------------------------------------------------------
             
    #saving TTree in ROOT file
    df_ZLH.Snapshot("tree_ZLH", "complete_ZLH.root")
    df_WLH.Snapshot("tree_WLH", "complete_WLH.root")

    df_ZTH = ROOT.RDataFrame(tree_ZTH)
    df_WTH = ROOT.RDataFrame(tree_WTH)
  
    print("Example ZTH dataframe")
    df_ZTH.Display().Print()

    print("Example WTH dataframe")    
    df_WTH.Display().Print()

    #creating TTree ROOT file
    df_ZTH.Snapshot("tree_ZTH", "complete_ZTH.root")
    df_WTH.Snapshot("tree_WTH", "complete_WTH.root")


if __name__ == '__main__':
    main()