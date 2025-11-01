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
    cut1_WL, cut2_WL, cut3_WL, cut4_WL, cut5_WL = [], [], [], [], []
    cut12_WL, cut123_WL, cut1234_WL = [], [], []
    cut13_WL, cut135_WL, cut1352_WL = [], [], []


    for e_id, part_dict in events_dict.items():

        if e_id % 10000 == 0:
            print(f"Processing event {e_id}")

        
        #if e_id > 5000:
            #break
        
        event_type = f_event_type(part_dict)
             
        #di-lepton
        if event_type == 23:

            total_ZL.append(e_id)

            fm_vec_1 = vector.obj(x = 0, y = 0, z = 0, E = 0)
            fm_vec_2 = vector.obj(x = 0, y = 0, z = 0, E = 0)
                
            for pid, cinematic in part_dict.items():
                if (pid) in [11, 13]:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                    fm_vec_1 = smeared_momentum
                if (pid) in [-11, -13]:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                    fm_vec_2 = smeared_momentum
                    
            if not is_nonzero(fm_vec_1) or not is_nonzero(fm_vec_2):
                continue
                
            fm_Z = fm_vec_1 + fm_vec_2

            e_mu_ZL.append(e_id)
            
            for pid, cinematic in part_dict.items():
                if (pid) == 5:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    jet1_fm = smeared_momentum
                    
                if (pid) == -5:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    jet2_fm = smeared_momentum
                    
            fm_jj = jet1_fm + jet2_fm
            
            #selections
            pt_V = np.sqrt(fm_Z.px**2 + fm_Z.py**2)
            delta_phi = fm_Z.deltaphi(fm_jj)

            cut1 = pt_V > 75
            cut2 = delta_phi > 2.5

            if cut1:
                cut1_ZL.append(e_id)

            if cut2:
                cut2_ZL.append(e_id)
            
            if cut1 and cut2:
                selected_ZL.append(e_id)

                
    
        #single-lepton
        if event_type == 24 or event_type == -24:

            total_WL.append(e_id)
              
            fm_vec_1 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
            fm_vec_2 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
            
            for pid, cinematic in part_dict.items():
                if abs(pid) in [11, 13]:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                    fm_vec_1 = smeared_momentum
                    
            for pid, cinematic in part_dict.items():
                if abs(pid) in [12, 14]:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "neutrino")
                    fm_vec_2 = compute_neutrino_momentum_from_particles(fm_vec_1, smeared_momentum)
                    
            if is_zero(fm_vec_1) or is_zero(fm_vec_2):
                continue

            fm_W = fm_vec_1 + fm_vec_2

            #if is_nonzero(fm_vec_1) or is_nonzero(fm_vec_2):
            e_mu_WL.append(e_id)

            #jets
            for pid, cinematic in part_dict.items():
                if (pid) == 5:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    jet1_fm = smeared_momentum
                    
                if (pid) == -5:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    jet2_fm = smeared_momentum
                    
            if is_zero(jet1_fm) or is_zero(jet2_fm):
                print(f"Null jet vector in event {e_id}")
                    
            fm_jj = jet1_fm + jet2_fm
            
            
            #selections
            MET = vector.obj(x = fm_vec_2.px, y =  fm_vec_2.py, z = 0)
            pt_W = np.sqrt(fm_W.px**2 + fm_W.py**2)
            pt_j1 = np.sqrt(jet1_fm.px**2 + jet1_fm.py**2)
            pt_j2 = np.sqrt(jet2_fm.px**2 + jet2_fm.py**2)
            pt_jj = np.sqrt(fm_jj.px**2 + fm_jj.py**2)            
            delta_phi_lep_met = fm_vec_1.deltaphi(MET) 
            delta_phi_w_jj = fm_W.deltaphi(fm_jj)

            cut1 = pt_W > 150
            cut2 = pt_j1 > 25 and pt_j2 > 25
            cut3 = pt_jj > 100
            cut4 = delta_phi_lep_met < 2.0
            cut5 = delta_phi_w_jj > 2.5

            if cut1:
                cut1_WL.append(e_id)

            if cut2:
                cut2_WL.append(e_id)

            if cut3:
                cut3_WL.append(e_id)

            if cut4:
                cut4_WL.append(e_id)

            if cut5:
                cut5_WL.append(e_id)

            if cut1 and cut2:
                cut12_WL.append(e_id)
            
            if cut1 and cut2 and cut3:
                cut123_WL.append(e_id)
            
            if cut1 and cut2 and cut3 and cut4:
                cut1234_WL.append(e_id)          

                        
            if cut1 and cut2 and cut3 and cut4 and cut5:
                selected_WL.append(e_id)

            #ordered partial -> 1, 3, 5, 2, 4
            if cut1 and cut3:
                cut13_WL.append(e_id)

            if cut1 and cut3 and cut5:
                cut135_WL.append(e_id)
                 
            if cut1 and cut3 and cut5 and cut2:
                cut1352_WL.append(e_id)



    #----------------------------------------------------------------------------------      

    print("\n Transverse polarisation \n")

    LHE_file = "unweighted_events_T.lhe"
    events_dict = read_file(LHE_file)

    total_ZT, e_mu_ZT = [], []
    total_WT, e_mu_WT = [], []
    selected_ZT, selected_WT = [], []
    cut1_ZT, cut2_ZT = [], []
    cut1_WT, cut2_WT, cut3_WT, cut4_WT, cut5_WT = [], [], [], [], []
    cut12_WT, cut123_WT, cut1234_WT = [], [], [] 
    cut13_WT, cut135_WT, cut1354_WT = [], [], []

    for e_id, part_dict in events_dict.items():

        
        if e_id % 10000 == 0:
            print(f"Processing event {e_id}")

        #if e_id > 5000:
            #break
        
        event_type = f_event_type(part_dict)
             
        #di-lepton
        if event_type == 23:

            total_ZT.append(e_id)

            fm_vec_1 = vector.obj(x = 0, y = 0, z = 0, E = 0)
            fm_vec_2 = vector.obj(x = 0, y = 0, z = 0, E = 0)
                
            for pid, cinematic in part_dict.items():
                if (pid) in [11, 13]:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                    fm_vec_1 = smeared_momentum
                if (pid) in [-11, -13]:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                    fm_vec_2 = smeared_momentum
                    
            if is_zero(fm_vec_1) or is_zero(fm_vec_2):
                #print(f"Null lepton vector in event {e_id}")
                continue
                
            fm_Z = fm_vec_1 + fm_vec_2
            e_mu_ZT.append(e_id)
            
            for pid, cinematic in part_dict.items():
                if (pid) == 5:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    jet1_fm = smeared_momentum
                    
                if (pid) == -5:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    jet2_fm = smeared_momentum
                    
            fm_jj = jet1_fm + jet2_fm
            
            #selections
            pt_V = np.sqrt(fm_Z.px**2 + fm_Z.py**2)
            delta_phi = fm_Z.deltaphi(fm_jj)

            cut1 = pt_V > 75
            cut2 = delta_phi > 2.5
            
            if cut1:
                cut1_ZT.append(e_id)

            if cut2:
                cut2_ZT.append(e_id)

            if cut1 and cut2:
                selected_ZT.append(e_id)
                
    
        #single-lepton
        if event_type == 24 or event_type == -24:

            total_WT.append(e_id)
            
            fm_vec_1 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
            fm_vec_2 = vector.obj(px = 0, py = 0, pz = 0, E = 0)
            
            for pid, cinematic in part_dict.items():
                if abs(pid) in [11, 13]:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "lepton")
                    fm_vec_1 = smeared_momentum
                    
            for pid, cinematic in part_dict.items():
                if abs(pid) in [12, 14]:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "neutrino")
                    fm_vec_2 = compute_neutrino_momentum_from_particles(fm_vec_1, smeared_momentum)
                    
            if is_zero(fm_vec_1) or is_zero(fm_vec_2):
                continue

            fm_W = fm_vec_1 + fm_vec_2
            e_mu_WT.append(e_id)

            #jets
            for pid, cinematic in part_dict.items():
                if (pid) == 5:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    jet1_fm = smeared_momentum
                    
                if (pid) == -5:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    jet2_fm = smeared_momentum
                    
            if not is_nonzero(jet1_fm) or not is_nonzero(jet2_fm):
                print(f"Null jet vector in event {e_id}")
                    
            fm_jj = jet1_fm + jet2_fm
            
            #selections
            MET = vector.obj(x = fm_vec_2.px, y =  fm_vec_2.py, z = 0)
            pt_W = np.sqrt(fm_W.px**2 + fm_W.py**2)
            pt_j1 = np.sqrt(jet1_fm.px**2 + jet1_fm.py**2)
            pt_j2 = np.sqrt(jet2_fm.px**2 + jet2_fm.py**2)
            pt_jj = np.sqrt(fm_jj.px**2 + fm_jj.py**2)            
            delta_phi_lep_met = fm_vec_1.deltaphi(MET) 
            delta_phi_w_jj = fm_W.deltaphi(fm_jj)

            cut1 = pt_W > 150
            cut2 = pt_j1 > 25 and pt_j2 > 25
            cut3 = pt_jj > 100
            cut4 = delta_phi_lep_met < 2.0
            cut5 = delta_phi_w_jj > 2.5

            if cut1:
                cut1_WT.append(e_id)
            
            if cut2:
                cut2_WT.append(e_id)
            
            if cut3:
                cut3_WT.append(e_id)
            
            if cut4:
                cut4_WT.append(e_id)

            if cut5:
                cut5_WT.append(e_id)

            if cut1 and cut2:
                cut12_WT.append(e_id)
            
            if cut1 and cut2 and cut3:
                cut123_WT.append(e_id)
            
            if cut1 and cut2 and cut3 and cut4:
                cut1234_WT.append(e_id)
            
            if cut1 and cut2 and cut3 and cut4 and cut5:
                selected_WT.append(e_id)            

            #partial eff -> 1, 3, 5, 4, 2
            if cut1 and cut3:
                cut13_WT.append(e_id)

            if cut1 and cut3 and cut5:
                cut135_WT.append(e_id)

            if cut1 and cut3 and cut5 and cut4:
                cut1354_WT.append(e_id)

            
    #----------------------------------------------------------------------------------      

            
    #ZH dataframe
    N_tot_ZL = len(total_ZL)
    N_em_ZL = len(e_mu_ZL)
    N_cut1_ZL = len(cut1_ZL)
    N_cut2_ZL = len(cut2_ZL)
    N_cut12_ZL = len(selected_ZL)

    N_tot_ZT = len(total_ZT)
    N_em_ZT = len(e_mu_ZT)
    N_cut1_ZT = len(cut1_ZT)
    N_cut2_ZT = len(cut2_ZT)
    N_cut12_ZT = len(selected_ZT)


    #wh dataframe
    N_tot_WL = len(total_WL)
    N_em_WL = len(e_mu_WL)
    N_cut1_WL = len(cut1_WL)
    N_cut2_WL = len(cut2_WL)
    N_cut3_WL = len(cut3_WL)
    N_cut4_WL = len(cut4_WL)
    N_cut5_WL = len(cut5_WL)
    N_cut12345_WL = len(selected_WL)

    N_cut12_WL = len(cut12_WL)
    N_cut123_WL = len(cut123_WL)
    N_cut1234_WL = len(cut1234_WL)

    #WTH
    N_tot_WT = len(total_WT)
    N_em_WT = len(e_mu_WT)
    N_cut1_WT = len(cut1_WT)
    N_cut2_WT = len(cut2_WT)
    N_cut3_WT = len(cut3_WT)
    N_cut4_WT = len(cut4_WT)
    N_cut5_WT = len(cut5_WT)
    N_cut12345_WT = len(selected_WT)

    N_cut12_WT = len(cut12_WT)
    N_cut123_WT = len(cut123_WT)
    N_cut1234_WT = len(cut1234_WT)

    #channel dataframe
    data = [
            ["Total events",           N_tot_ZL,   N_tot_ZT,   N_tot_WL,      N_tot_WT],
            ["Events with e/μ",        N_em_ZL,    N_em_ZT,    N_em_WL,       N_em_WT],
            ["Events after cut1 only", N_cut1_ZL,  N_cut1_ZT,  N_cut1_WL,     N_cut1_WT],
            ["Events after cut2 only", N_cut2_ZL,  N_cut2_ZT,  N_cut2_WL,     N_cut2_WT],
            ["Events after cut3 only", "-",        "-",        N_cut3_WL,     N_cut3_WT],
            ["Events after cut4 only", "-",        "-",        N_cut4_WL,     N_cut4_WT],
            ["Events after cut5 only", "-",        "-",        N_cut5_WL,     N_cut5_WT],
            ["Selected events",        N_cut12_ZL, N_cut12_ZT, N_cut12345_WL, N_cut12345_WT]
            ]


    df_channel = pd.DataFrame(data, columns=["", "ZLH", "ZTH", "WLH", "WTH"])

    print("\n" + "="*60)
    print("{:^60}".format("H > b b~"))
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
    epsilon_WL  = N_tot_WL / 50000
    epsilon_0_WL = N_em_WL / 50000
    epsilon_1_WL = N_cut1_WL / N_em_WL
    epsilon_2_WL = N_cut2_WL / N_em_WL
    epsilon_3_WL = N_cut3_WL / N_em_WL
    epsilon_4_WL = N_cut4_WL / N_em_WL
    epsilon_5_WL = N_cut5_WL / N_em_WL
    epsilon_tot_WL = N_cut12345_WL / N_tot_WL

    #partial efficiencies
    part_epsilon_0_WL = N_em_WL / N_tot_WL
    part_epsilon_1_WL = N_cut1_WL / N_em_WL
    part_epsilon_2_WL = N_cut12_WL / N_cut1_WL
    part_epsilon_3_WL = N_cut123_WL / N_cut12_WL
    part_epsilon_4_WL = N_cut1234_WL / N_cut123_WL
    part_epsilon_5_WL = N_cut12345_WL / N_cut1234_WL


    #WTH
    #total efficiencies
    epsilon_WT = N_tot_WT / 50000
    epsilon_0_WT = N_em_WT / 50000
    epsilon_1_WT = N_cut1_WT / N_em_WT
    epsilon_2_WT = N_cut2_WT / N_em_WT
    epsilon_3_WT = N_cut3_WT / N_em_WT
    epsilon_4_WT = N_cut4_WT / N_em_WT
    epsilon_5_WT = N_cut5_WT / N_em_WT
    epsilon_tot_WT = N_cut12345_WT / N_tot_WT

    #partial efficiencies
    part_epsilon_0_WT = N_em_WT / N_tot_WT
    part_epsilon_1_WT = N_cut1_WT / N_em_WT
    part_epsilon_2_WT = N_cut12_WT / N_cut1_WT
    part_epsilon_3_WT = N_cut123_WT / N_cut12_WT
    part_epsilon_4_WT = N_cut1234_WT / N_cut123_WT
    part_epsilon_5_WT = N_cut12345_WT / N_cut1234_WT

    

    efficiencies = np.array([
                            ["channel fraction",        f"{epsilon_ZL:.2f}",   f"{epsilon_ZT:.2f}",   f"{epsilon_WL:.2f}",   f"{epsilon_WT:.2f}"],
                            ["ε₀ = e/μ / total",        f"{epsilon_0_ZL:.2f}", f"{epsilon_0_ZT:.2f}", f"{epsilon_0_WL:.2f}", f"{epsilon_0_WT:.2f}"],
                            ["ε₁ = cut1 / e/μ",         f"{epsilon_1_ZL:.2f}", f"{epsilon_1_ZT:.2f}", f"{epsilon_1_WL:.2f}", f"{epsilon_1_WT:.2f}"],
                            ["ε₂ = cut2 / e/μ",         f"{epsilon_2_ZL:.2f}", f"{epsilon_2_ZT:.2f}", f"{epsilon_2_WL:.2f}", f"{epsilon_2_WT:.2f}"],
                            ["ε₃ = cut3 / e/μ",         "-",                   "-",                   f"{epsilon_3_WL:.2f}", f"{epsilon_3_WT:.2f}"],
                            ["ε₄ = cut4 / e/μ",         "-",                   "-",                   f"{epsilon_4_WL:.2f}", f"{epsilon_4_WT:.2f}"],
                            ["ε₅ = cut5 / e/μ",         "-",                   "-",                   f"{epsilon_5_WL:.2f}", f"{epsilon_5_WT:.2f}"],
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
                                ["ε₄ (cut4)",     "-",                        "-",                   "-",                        "-",                   f"{part_epsilon_4_WL:.2f}", f"{epsilon_4_WL:.2f}", f"{part_epsilon_4_WT:.2f}", f"{epsilon_4_WT:.2f}"],
                                ["ε₅ (cut5)",     "-",                        "-",                   "-",                        "-",                   f"{part_epsilon_5_WL:.2f}", f"{epsilon_5_WL:.2f}", f"{part_epsilon_5_WT:.2f}", f"{epsilon_5_WT:.2f}"],
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


    #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    #computing LHC expected events number

    #xsection_L = 0.42    #fb
    xsection_L = 282    #fb
    BR = 0.5824
    #xsection_T = 0.3739  #fb
    xsection_T = 252  #fb

    lum = 100            #1/fb 

    # ZLH
    N_tot_L_LHC     = xsection_L * lum                                      #total longitudinal events number
    N_tot_ZL_LHC    = N_tot_L_LHC * epsilon_ZL                              #ZLH events number
    N_tot_ZL_bb     = N_tot_ZL_LHC * BR                                     #h > bb events number
    N_em_ZL_LHC     = N_tot_ZL_bb * part_epsilon_0_ZL
    N_cut1_ZL_LHC   = N_em_ZL_LHC * part_epsilon_1_ZL
    N_cut2_ZL_LHC   = N_cut1_ZL_LHC * part_epsilon_2_ZL
    kin_eff_ZL_LHC  = N_cut2_ZL_LHC / N_em_ZL_LHC                           #kinematic cuts efficiency
    chan_eff_ZL_LHC = N_cut2_ZL_LHC / N_tot_ZL_bb                           #channel efficiency
    tot_eff_ZL_LHC  = N_cut2_ZL_LHC / N_tot_ZL_LHC


    # ZTH
    N_tot_T_LHC     = xsection_T * lum
    N_tot_ZT_LHC    = N_tot_T_LHC * epsilon_ZT
    N_tot_ZT_bb     = N_tot_ZT_LHC * BR
    N_em_ZT_LHC     = N_tot_ZT_bb * part_epsilon_0_ZT
    N_cut1_ZT_LHC   = N_em_ZT_LHC * part_epsilon_1_ZT
    N_cut2_ZT_LHC   = N_cut1_ZT_LHC * part_epsilon_2_ZT
    kin_eff_ZT_LHC  = N_cut2_ZT_LHC / N_em_ZT_LHC
    chan_eff_ZT_LHC = N_cut2_ZT_LHC / N_tot_ZT_bb                           #channel efficiency
    tot_eff_ZT_LHC  = N_cut2_ZT_LHC / N_tot_ZT_LHC

    # WLH
    N_tot_L_LHC    = xsection_L * lum
    N_tot_WL_LHC    = N_tot_L_LHC * epsilon_WL
    N_tot_WL_bb     = N_tot_WL_LHC * BR
    N_em_WL_LHC     = N_tot_WL_bb * part_epsilon_0_WL
    N_cut1_WL_LHC   = N_em_WL_LHC * part_epsilon_1_WL
    N_cut2_WL_LHC   = N_cut1_WL_LHC * part_epsilon_2_WL
    N_cut3_WL_LHC   = N_cut2_WL_LHC * part_epsilon_3_WL
    N_cut4_WL_LHC   = N_cut3_WL_LHC * part_epsilon_4_WL
    N_cut5_WL_LHC   = N_cut4_WL_LHC * part_epsilon_5_WL
    kin_eff_WL_LHC  = N_cut5_WL_LHC / N_em_WL_LHC
    chan_eff_WL_LHC = N_cut5_WL_LHC / N_tot_WL_bb                           #channel efficiency
    tot_eff_WL_LHC  = N_cut5_WL_LHC / N_tot_WL_LHC

    # WTH
    N_tot_T_LHC     = xsection_T * lum
    N_tot_WT_LHC    = N_tot_T_LHC * epsilon_WT
    N_tot_WT_bb     = N_tot_WT_LHC * BR
    N_em_WT_LHC     = N_tot_WT_bb * part_epsilon_0_WT
    N_cut1_WT_LHC   = N_em_WT_LHC * part_epsilon_1_WT
    N_cut2_WT_LHC   = N_cut1_WT_LHC * part_epsilon_2_WT
    N_cut3_WT_LHC   = N_cut2_WT_LHC * part_epsilon_3_WT
    N_cut4_WT_LHC   = N_cut3_WT_LHC * part_epsilon_4_WT
    N_cut5_WT_LHC   = N_cut4_WT_LHC * part_epsilon_5_WT

    kin_eff_WT_LHC  = N_cut5_WT_LHC / N_em_WT_LHC
    chan_eff_WT_LHC = N_cut5_WT_LHC / N_tot_WT_bb                           #channel efficiency
    tot_eff_WT_LHC  = N_cut5_WT_LHC / N_tot_WT_LHC

    data = [
            ["Total L/T events",                  f"{N_tot_L_LHC:.0f}",     f"{N_tot_T_LHC:.0f}",     f"{N_tot_L_LHC:.0f}",     f"{N_tot_T_LHC:.0f}"],
            ["Total V{i}H events",                  f"{N_tot_ZL_LHC:.0f}",    f"{N_tot_ZT_LHC:.0f}",    f"{N_tot_WL_LHC:.0f}",    f"{N_tot_WT_LHC:.0f}"],
            ["Total bb events",                     f"{N_tot_ZL_bb:.0f}",     f"{N_tot_ZT_bb:.0f}",     f"{N_tot_WL_bb:.0f}",     f"{N_tot_WT_bb:.0f}"],
            ["Events with e/mu",                    f"{N_em_ZL_LHC:.0f}",     f"{N_em_ZT_LHC:.0f}",     f"{N_em_WL_LHC:.0f}",     f"{N_em_WT_LHC:.0f}"],
            ["Events after cut1",                   f"{N_cut1_ZL_LHC:.0f}",   f"{N_cut1_ZT_LHC:.0f}",   f"{N_cut1_WL_LHC:.0f}",   f"{N_cut1_WT_LHC:.0f}"],
            ["Events after cut1, cut2",             f"{N_cut2_ZL_LHC:.0f}",   f"{N_cut2_ZT_LHC:.0f}",   f"{N_cut2_WL_LHC:.0f}",   f"{N_cut2_WT_LHC:.0f}"],
            ["Events after cut1, cut2, cut3",       "-",                      "-",                      f"{N_cut3_WL_LHC:.0f}",   f"{N_cut3_WT_LHC:.0f}"],
            ["Events after cut1, cut2, cut3, cut4", "-",                      "-",                      f"{N_cut4_WL_LHC:.0f}",   f"{N_cut4_WT_LHC:.0f}"],
            ["Events after cut1, cut2, cut3, cut4, cut5", "-",                "-",                      f"{N_cut5_WL_LHC:.0f}",   f"{N_cut5_WT_LHC:.0f}"],
            ["Kinematic cuts efficiency",           f"{kin_eff_ZL_LHC:.2f}",  f"{kin_eff_ZT_LHC:.2f}",  f"{kin_eff_WL_LHC:.2f}",  f"{kin_eff_WT_LHC:.2f}"],
            ["Channel efficiency",                  f"{chan_eff_ZL_LHC:.2f}", f"{chan_eff_ZT_LHC:.2f}", f"{chan_eff_WL_LHC:.2f}", f"{chan_eff_WT_LHC:.2f}"],
            ["Total efficiency",                    f"{tot_eff_ZL_LHC:.2f}",  f"{tot_eff_ZT_LHC:.2f}",  f"{tot_eff_WL_LHC:.2f}",  f"{tot_eff_WT_LHC:.4f}"]
        ]
    
    df_channel = pd.DataFrame(data, columns=["", "ZLH", "ZTH", "WLH", "WTH"])

    print("\n" + "="*60)
    print("{:^60}".format("H > bb~ (LHC)"))
    print("="*60 + "\n")
    print(tabulate(df_channel.values.tolist(), headers=df_channel.columns.tolist(), tablefmt="grid"))


    #----------------------------------------------------------------------------------------------------------------------------------------------------------------


    print("Ordering cuts")

    cut_list = ["cut1", "cut2", "cut3", "cut4", "cut5"]

    for col in column[1:]:
        
        #re-arranging
        values = df_efficiencies.loc[2:6, col].copy()
        num_values = [float(v) for v in values if v != "-"]

        ordered_val = sorted(num_values)
        print("Ordered total efficiencies:", ordered_val)
        sorted_indices = sorted(range(len(num_values)), key=lambda i: num_values[i])

        ordered_cut_list = [cut_list[i] for i in sorted_indices]
        print(f"Ordered cuts for {col}:", ordered_cut_list)

        j = 0
        for m in range(2, 7): 

            if df_efficiencies.at[m, col] == "-":
                continue

            df_efficiencies.at[m, col] = f"{ordered_val[j]:.2f}"
            j += 1

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    #ordering cumulative efficiencies

    N_cut13_WL = len(cut13_WL)
    N_cut135_WL = len(cut135_WL)
    N_cut1352_WL = len(cut1352_WL)

    N_cut13_WT = len(cut13_WT)
    N_cut135_WT = len(cut135_WT)
    N_cut1354_WT = len(cut1354_WT)

    #ZLH
    #partial efficiencies
    part_epsilon_1_ZL = N_cut1_ZL / N_em_ZL
    part_epsilon_2_ZL = N_cut12_ZL / N_cut1_ZL

    #ZTH
    part_epsilon_1_ZT = N_cut1_ZT / N_em_ZT
    part_epsilon_2_ZT = N_cut12_ZT / N_cut1_ZT

    #WLH
    #partial efficiencies
    part_epsilon_1_WL = N_cut1_WL / N_em_WL
    part_epsilon_2_WL = N_cut13_WL / N_cut1_WL
    part_epsilon_3_WL = N_cut135_WL / N_cut13_WL
    part_epsilon_4_WL = N_cut1352_WL / N_cut135_WL
    part_epsilon_5_WL = N_cut12345_WL / N_cut1352_WL

    #WTH
    #partial efficiencies
    part_epsilon_1_WT = N_cut1_WT / N_em_WT
    part_epsilon_2_WT = N_cut13_WT / N_cut1_WT
    part_epsilon_3_WT = N_cut135_WT / N_cut13_WT
    part_epsilon_4_WT = N_cut1354_WT / N_cut135_WT
    part_epsilon_5_WT = N_cut12345_WT / N_cut1354_WT
    

    part_total_eff_ordered = np.array([
                                ["ε₀ (e/μ)",            f"{part_epsilon_0_ZL:.2f}", f"{epsilon_0_ZL:.2f}", f"{part_epsilon_0_ZT:.2f}",      f"{epsilon_0_ZT:.2f}", f"{part_epsilon_0_WL:.2f}",      f"{epsilon_0_WL:.2f}", f"{part_epsilon_0_WT:.2f}",      f"{epsilon_0_WT:.2f}"],
                                ["ε₁ (cut1)",           f"{part_epsilon_1_ZL:.2f}", f"{epsilon_1_ZL:.2f}", f"{part_epsilon_1_ZT:.2f}", f"{epsilon_1_ZT:.2f}", f"{part_epsilon_1_WL:.2f}", f"{epsilon_1_WL:.2f}", f"{part_epsilon_1_WT:.2f}", f"{epsilon_1_WT:.2f}"],
                                ["ε₂ (cut2)",           f"{part_epsilon_2_ZL:.2f}", f"{epsilon_2_ZL:.2f}", f"{part_epsilon_2_ZT:.2f}", f"{epsilon_2_ZT:.2f}", f"{part_epsilon_2_WL:.2f}", f"{epsilon_3_WL:.2f}", f"{part_epsilon_2_WT:.2f}", f"{epsilon_3_WT:.2f}"],
                                ["ε₃ (cut3)",           "-",                        "-",                   "-",                        "-",                   f"{part_epsilon_3_WL:.2f}", f"{epsilon_5_WL:.2f}", f"{part_epsilon_3_WT:.2f}", f"{epsilon_5_WT:.2f}"],
                                ["ε₄ (cut4)",           "-",                        "-",                   "-",                        "-",                   f"{part_epsilon_4_WL:.2f}", f"{epsilon_2_WL:.2f}", f"{part_epsilon_4_WT:.2f}", f"{epsilon_4_WT:.2f}"],
                                ["ε₅ (cut5)",           "-",                        "-",                   "-",                        "-",                   f"{part_epsilon_5_WL:.2f}", f"{epsilon_4_WL:.2f}", f"{part_epsilon_5_WT:.2f}", f"{epsilon_2_WT:.2f}"],
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

   

if __name__ == '__main__':
    main()