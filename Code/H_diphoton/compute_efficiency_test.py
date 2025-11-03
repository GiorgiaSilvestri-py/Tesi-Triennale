import pandas as pd
import pickle
from array import array
from tabulate import tabulate
import vector
from utils_nonlhe import *

def compute_deltaR(vec1, vec2):

    deltaphi = vec1.deltaphi(vec2)
    deltaeta = abs(vec1.eta - vec2.eta)

    return np.sqrt( (deltaphi)**2 + (deltaeta)**2 )

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


    for e_id, particle_list in events_dict.items():

        if e_id % 10000 == 0:
            print(f"Processing event {e_id}")

        if e_id > 1000:    break

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

            
            e_mu_WL.append(e_id)

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
            
            pt_MET = np.sqrt(neutrino_vec.px**2 + neutrino_vec.py**2)

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

         
    df_ZLH = ROOT.RDataFrame(tree_ZLH)
  
    df_WLH = ROOT.RDataFrame(tree_WLH)
  
    df_ZLH.Display().Print()
    df_WLH.Display().Print()
    
    N_tot_ZLH = df_ZLH.Count().GetValue()
    N_tot_WLH = df_WLH.Count().GetValue()
    
    df_ZLH = df_ZLH\
            .Define("cut0_ZLH", '(event_type == "e") || (event_type == "mu")')\
            .Define("cut1_ZLH", 'pt_leading > leading_th && pt_subleading > subleading_th')\
            .Define("cut2_ZLH", '(event_type == "e" && R1 > 1.0 && R2 > 1.0) || (event_type == "mu" && R1 > 0.5 && R2 > 0.5)')
    
    df_WLH = df_WLH\
            .Define("cut0_WLH", 'event_type == "e/mu"')\
            .Define("cut1_WLH", 'pt_leading > leading_th && pt_subleading > subleading_th')\
            .Define("cut2_WLH", 'event_type == "e/mu" && R1 > 1.0 && R2 > 1.0')\
            .Define("cut3_WLH", 'pt_MET > 45')
    
    
    N_lep_ZLH    = df_ZLH.Filter("cut0_ZLH").Count().GetValue()
    N_cut1_ZLH   = df_ZLH.Filter("cut0_ZLH && cut1_ZLH").Count().GetValue()
    N_cut2_ZLH   = df_ZLH.Filter("cut0_ZLH && cut2_ZLH").Count().GetValue()
    N_cut12_ZLH  = df_ZLH.Filter("cut0_ZLH && cut1_ZLH && cut2_ZLH").Count().GetValue()

    N_lep_WLH    = df_WLH.Filter("cut0_WLH").Count().GetValue()
    N_cut1_WLH   = df_WLH.Filter("cut0_WLH && cut1_WLH").Count().GetValue()
    N_cut2_WLH   = df_WLH.Filter("cut0_WLH && cut2_WLH").Count().GetValue()
    N_cut3_WLH   = df_WLH.Filter("cut0_WLH && cut3_WLH").Count().GetValue()
    N_cut123_WLH = df_WLH.Filter("cut0_WLH && cut1_WLH && cut2_WLH && cut3_WLH").Count().GetValue()

    

    table = ROOT.TTree("Summary", "Event counts per channel")

    label = ROOT.std.string()
    ZLH = array('i', [0])
    WLH = array('i', [0])

    table.Branch("label", label)
    table.Branch("ZLH", ZLH, "ZLH/I")
    table.Branch("WLH", WLH, "WLH/I")

    def fill_row(name, zlh_val, wlh_val):
        label.assign(name)
        ZLH[0] = zlh_val
        WLH[0] = wlh_val
        table.Fill()

    fill_row("Total events",            N_tot_ZLH,   N_tot_WLH)
    fill_row("Events with e/mu",        N_lep_ZLH,   N_lep_WLH)
    fill_row("Events after cut1 only",  N_cut1_ZLH,  N_cut1_WLH)
    fill_row("Events after cut2 only",  N_cut2_ZLH,  N_cut2_WLH)
    fill_row("Events after cut3 only",  -1,          N_cut3_WLH)
    fill_row("Eventi selezionati",      N_cut12_ZLH, N_cut123_WLH)

    df_table = ROOT.RDataFrame(table)
    df_np = df_table.AsNumpy(["label", "ZLH", "WLH"])
    df_pd = pd.DataFrame(df_np)
    print(df_pd)

    print("Entries in table:", table.GetEntries())

    for entry in table:
        print(f"{entry.label} → ZLH: {entry.ZLH}, WLH: {entry.WLH}")

    
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


    for e_id, particle_list in events_dict.items():

        if e_id % 10000 == 0:
            print(f"Processing event {e_id}")

        if e_id > 1000: break

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
            m_gg   = sum_fm.M

            if gamma1_fm.pt == 0 :                      #check
                print("Null pt1 in event", e_id)
            
            if gamma2_fm.pt == 0 :
                print("Null pt1 in event", e_id)
            
            pt_MET = np.sqrt(neutrino_vec.px**2 + neutrino_vec.py**2)

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


             
    df_ZTH = ROOT.RDataFrame(tree_ZTH)
  
    df_WTH = ROOT.RDataFrame(tree_WTH)
  
    df_ZTH.Display().Print()
    df_WTH.Display().Print()
    
    N_tot_ZTH = df_ZTH.Count().GetValue()
    N_tot_WTH = df_WTH.Count().GetValue()

    df_ZTH = df_ZTH\
            .Define("cut0_ZTH", '(event_type == "e") || (event_type == "mu")')\
            .Define("cut1_ZTH", 'pt_leading > leading_th && pt_subleading > subleading_th')\
            .Define("cut2_ZTH", '(event_type == "e" && R1 > 1.0 && R2 > 1.0) || (event_type == "mu" && R1 > 0.5 && R2 > 0.5)')
    
    df_WTH = df_WTH\
            .Define("cut0_WTH", 'event_type == "e/mu"')\
            .Define("cut1_WTH", 'pt_leading > leading_th && pt_subleading > subleading_th')\
            .Define("cut2_WTH", 'event_type == "e/mu" && R1 > 1.0 && R2 > 1.0')\
            .Define("cut3_WTH", 'pt_MET > 45')
    
    
    N_lep_ZTH    = df_ZTH.Filter("cut0_ZTH").Count().GetValue()
    N_cut1_ZTH   = df_ZTH.Filter("cut0_ZTH && cut1_ZTH").Count().GetValue()
    N_cut2_ZTH   = df_ZTH.Filter("cut0_ZTH && cut2_ZTH").Count().GetValue()
    N_cut12_ZTH  = df_ZTH.Filter("cut0_ZTH && cut1_ZTH && cut2_ZTH").Count().GetValue()

    N_lep_WTH    = df_WTH.Filter("cut0_WTH").Count().GetValue()
    N_cut1_WTH   = df_WTH.Filter("cut0_WTH && cut1_WTH").Count().GetValue()
    N_cut2_WTH   = df_WTH.Filter("cut0_WTH && cut2_WTH").Count().GetValue()
    N_cut3_WTH   = df_WTH.Filter("cut0_WTH && cut3_WTH").Count().GetValue()
    N_cut123_WTH = df_WTH.Filter("cut0_WTH && cut1_WTH && cut2_WTH && cut3_WTH").Count().GetValue()

    
    table = ROOT.TTree("Summary", "Event counts per channel")

    label = ROOT.std.string()
    ZTH = array('i', [0])
    WTH = array('i', [0])

    table.Branch("label", label)
    table.Branch("ZTH", ZTH, "ZTH/I")
    table.Branch("WTH", WTH, "WTH/I")

    def fill_row(name, zth_val, wth_val):
        label.assign(name)
        ZTH[0] = zth_val
        WTH[0] = wth_val
        table.Fill()

    fill_row("Total events",            N_tot_ZTH,   N_tot_WTH)
    fill_row("Events with e/mu",        N_lep_ZTH,   N_lep_WTH)
    fill_row("Events after cut1 only",  N_cut1_ZTH,  N_cut1_WTH)
    fill_row("Events after cut2 only",  N_cut2_ZTH,  N_cut2_WTH)
    fill_row("Events after cut3 only",  -1,          N_cut3_WTH)
    fill_row("Selected events",         N_cut12_ZTH, N_cut123_WTH)

    df_table = ROOT.RDataFrame(table)
    df_np = df_table.AsNumpy(["label", "ZTH", "WTH"])
    df_pd = pd.DataFrame(df_np)
    print(df_pd)

    print("Entries in table:", table.GetEntries())

    for entry in table:
        print(f"{entry.label} → ZTH: {entry.ZTH}, WTH: {entry.WTH}")
    
    #--------------------------------------------------------------------------------------------

    '''
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
    '''

    #channel dataframe
    data = [
            ["Total events",            N_tot_ZLH,   N_tot_ZTH,   N_tot_WLH,    N_tot_WTH],
            ["Events with e/mu",        N_lep_ZLH,   N_lep_ZTH,   N_lep_WLH,    N_lep_WTH],
            ["Events after cut1 only",  N_cut1_ZLH,  N_cut1_ZTH,  N_cut1_WLH,   N_cut1_WTH],
            ["Events after cut2 only",  N_cut2_ZLH,  N_cut2_ZTH,  N_cut2_WLH,   N_cut2_WTH],
            ["Events after cut3 only",  "-",         "-",         N_cut3_WLH,   N_cut3_WTH],
            ["Eventi selezionati",      N_cut12_ZLH, N_cut12_ZTH, N_cut123_WLH, N_cut123_WTH]
            ]


    df_channel = pd.DataFrame(data, columns=["", "ZLH", "ZTH", "WLH", "WTH"])

    print("\n" + "="*60)
    print("{:^60}".format("H > γγ"))
    print("="*60 + "\n")
    print(tabulate(df_channel.values.tolist(), headers=df_channel.columns.tolist(), tablefmt="grid"))

    N_lep_ZLH    = df_ZLH.Filter("cut0_ZLH").Count().GetValue()
    N_cut1_ZLH   = df_ZLH.Filter("cut0_ZLH && cut1_ZLH").Count().GetValue()
    N_cut12_ZLH  = df_ZLH.Filter("cut0_ZLH && cut1_ZLH && cut2_ZLH").Count().GetValue()

    N_lep_WLH    = df_WLH.Filter("cut0_WLH").Count().GetValue()
    N_cut1_WLH   = df_WLH.Filter("cut0_WLH && cut1_WLH").Count().GetValue()
    N_cut12_WLH   = df_WLH.Filter("cut0_WLH && cut1_WLH && cut2_WLH").Count().GetValue()
    N_cut123_WLH = df_WLH.Filter("cut0_WLH && cut1_WLH && cut2_WLH && cut3_WLH").Count().GetValue()

    N_lep_ZTH    = df_ZTH.Filter("cut0_ZTH").Count().GetValue()
    N_cut1_ZTH   = df_ZTH.Filter("cut0_ZTH && cut1_ZTH").Count().GetValue()
    N_cut12_ZTH  = df_ZTH.Filter("cut0_ZTH && cut1_ZTH && cut2_ZTH").Count().GetValue()

    N_lep_WTH    = df_WTH.Filter("cut0_WTH").Count().GetValue()
    N_cut1_WTH   = df_WTH.Filter("cut0_WTH && cut1_WTH").Count().GetValue()
    N_cut12_WTH   = df_WTH.Filter("cut0_WTH && cut1_WTH && cut2_WTH").Count().GetValue()
    N_cut123_WTH = df_WTH.Filter("cut0_WTH && cut1_WTH && cut2_WTH && cut3_WTH").Count().GetValue()

    data_cumulative = [
                    ["Total events",            N_tot_ZLH,   N_tot_ZTH,   N_tot_WLH,    N_tot_WTH],
                    ["Events with e/mu",        N_lep_ZLH,   N_lep_ZTH,   N_lep_WLH,    N_lep_WTH],
                    ["Events after cut1",       N_cut1_ZLH,  N_cut1_ZTH,  N_cut1_WLH,   N_cut1_WTH],
                    ["Events after cut2",       N_cut12_ZLH, N_cut12_ZTH, N_cut12_WLH,  N_cut12_WTH],
                    ["Events after cut3",       "-",         "-",         N_cut123_WLH, N_cut123_WTH],
                    ]
    
    df_cumulative = pd.DataFrame(data_cumulative, columns=["", "ZLH", "ZTH", "WLH", "WTH"])

    print("\n" + "="*60)
    print("{:^60}".format("H > γγ"))
    print("="*60 + "\n")
    print(tabulate(df_cumulative.values.tolist(), headers=df_cumulative.columns.tolist(), tablefmt="grid"))
    #-------------------------------------------------------------------------------------------------------


    print("Computing efficiencies")
    #ZLH
    #total efficiencies
    tot_eff_list_channel  = []
    part_eff_list_channel = []

    for ch in ["ZLH", "ZTH", "WLH", "WTH"]:
        epsilon     = df_channel.iloc[0][ch] / 5000
        epsilon_0   = df_channel.iloc[1][ch] / 5000
        epsilon_1   = df_channel.iloc[2][ch] / df_channel.iloc[1][ch]
        epsilon_2   = df_channel.iloc[3][ch] / df_channel.iloc[1][ch]

        part_epsilon_0  = df_cumulative.iloc[1][ch] / df_cumulative.iloc[0][ch]
        part_epsilon_1  = df_cumulative.iloc[2][ch] / df_cumulative.iloc[1][ch]
        part_epsilon_2  = df_cumulative.iloc[3][ch] / df_cumulative.iloc[2][ch]

        if df_channel.iloc[4][ch] == "-":
            epsilon_3      = -1
            part_epsilon_3 = -1
        else:
            epsilon_3      = df_channel.iloc[4][ch]    / df_channel.iloc[2][ch]
            part_epsilon_3 = df_cumulative.iloc[4][ch] / df_cumulative.iloc[3][ch]

        epsilon_tot = df_channel.iloc[5][ch] / df_channel.iloc[0][ch]

        tot_eff_list_channel.append(
                                [epsilon, epsilon_0, epsilon_1, epsilon_2, epsilon_3, epsilon_tot]
                                )
        part_eff_list_channel.append(
                                [epsilon, part_epsilon_0, part_epsilon_1, part_epsilon_2, part_epsilon_3, epsilon_tot]
                                )
        
    tot_eff_array_channel = np.array(tot_eff_list_channel)
    part_eff_array_channel = np.array(part_eff_list_channel)

    columns = ["",
               "ZLH cumulative eff", "ZLH total eff", 
               "ZTH cumulative eff", "ZTH total eff",
               "WLH cumulative eff", "WLH total eff",
               "WTH cumulative eff", "WTH total eff"
               ]

    rows = []

    for eps in range(6): #efficiencies

        epsilon_row = []
        if eps == 0:
                epsilon_row.append("epsilon channel")
        elif eps == 5:
            epsilon_row.append("epsilon tot")
        else:
            epsilon_row.append(f"epsilon_{eps - 1}")

        for ch in range(4): #channels

            partial = part_eff_array_channel[ch][eps]
            total = tot_eff_array_channel[ch][eps]

            

            epsilon_row.append(partial)
            epsilon_row.append(total)

        rows.append(epsilon_row)
        

    df_f = pd.DataFrame(rows, columns = columns)

    print("\n" + "="*120)
    print("{:^120}".format("Cumulative and Total efficiencies"))
    print("="*120 + "\n")

    print(tabulate(df_f.values.tolist(), headers=df_f.columns.tolist(), tablefmt="grid"))

    return 0

    print("Ordering cuts")

    cut_list = ["cut1", "cut2", "cut3"]
    channels = ["ZLH total eff", "ZTH total eff", "WLH total eff", "WTH total eff"]
    columns = ["",
                "ZLH cumulative eff", "ZLH total eff", 
                "ZTH cumulative eff", "ZTH total eff",
                "WLH cumulative eff", "WLH total eff",
                "WTH cumulative eff", "WTH total eff"
                ]
                
    for col in columns[2::2]:
        
        #re-arranging
        values = df_f.loc[2:4, col].copy()
        num_values = [float(v) for v in values if v != -1]

        ordered_val = sorted(num_values)
        print("Ordered total efficiencies:", ordered_val)
        sorted_indices = sorted(range(len(num_values)), key=lambda i: num_values[i])

        ordered_cut_list = [cut_list[i] for i in sorted_indices]
        print(f"Ordered cuts for {col}:", ordered_cut_list)

        j = 0
        for m in range(2, 5): 

            if df_f.at[m, col] == -1:
                continue

            df_f.at[m, col] = f"{ordered_val[j]:.2f}"
            j += 1

    #for col in columns[1::2]:


    '''
    epsilon_ZL = df_channel[1, "ZLH"] / 50000
    #N_tot_ZL / 50000
    epsilon_0_ZL = df_channel[2, "ZLH"] / 50000

    epsilon_1_ZL = df_channel[3, "ZLH"] / df_channel[2, "ZLH"]
    epsilon_2_ZL = df_channel[4, "ZLH"] / df_channel[2, "ZLH"]
    epsilon_tot_ZL = df_channel[5, "ZLH"] / df_channel[1, "ZLH"]

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
    '''

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