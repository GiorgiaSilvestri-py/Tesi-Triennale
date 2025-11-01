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
    cut1_WL, cut2_WL = [], []


    for e_id, particle_list in events_dict.items():

        if e_id % 10000 == 0:
            print(f"Processing event {e_id}")

        #if e_id > 5000: break

        event_type = f_event_type(particle_list)

        #di-lepton
        if event_type == 23:

            total_ZL.append(e_id)

            e_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            mu_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
           
            for cinematic in particle_list:
                #leptons
                if abs(cinematic.pid) == 11:
                    e_vec = cinematic.build_fm()

                if  abs(cinematic.pid) == 13:
                    mu_vec = cinematic.build_fm()

            if is_zero(e_vec) and is_zero(mu_vec):
                continue

            e_mu_ZL.append(e_id)

            tau_list, antitau_list = [], []
            for cinematic in particle_list:
                if cinematic.pid == 15:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    fm_vec_1 = smeared_momentum
                    tau_list.append(fm_vec_1)

                if cinematic.pid == -15:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    fm_vec_2 = smeared_momentum
                    antitau_list.append(fm_vec_2)
            
            for fm_tau in tau_list:

                for fm_anti_tau in antitau_list:
                    fm_sum_0 = fm_tau + fm_anti_tau
                    tollerance = 10

                    if (fm_sum_0.M - 125) < tollerance:
                        fm_sum = fm_sum_0
                    

            #selections
            pt_ditau = np.sqrt(fm_sum.px**2 + fm_sum.py**2)
            eta_ditau = fm_sum.eta

            cut1 = pt_ditau > 40
            cut2 = abs(eta_ditau) < 2.1
            
            if cut1:
                cut1_ZL.append(e_id)
            if cut2:
                cut2_ZL.append(e_id)

            if cut1 and cut2:
                selected_ZL.append(e_id)


        #di-lepton
        if abs(event_type) == 24:

            total_WL.append(e_id)

            e_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            mu_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            
            for cinematic in particle_list:
                #leptons
                if abs(cinematic.pid) == 11:
                    e_vec = cinematic.build_fm()

                if  abs(cinematic.pid) == 13:
                    mu_vec = cinematic.build_fm()
                    
            if is_zero(e_vec) and is_zero(mu_vec):
                continue

            e_mu_WL.append(e_id)

            tau_list, antitau_list = [], []
            for cinematic in particle_list:
                if cinematic.pid == 15:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    fm_vec_1 = smeared_momentum
                    tau_list.append(fm_vec_1)

                if cinematic.pid == -15:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    fm_vec_2 = smeared_momentum
                    antitau_list.append(fm_vec_2)
            
            for fm_tau in tau_list:

                for fm_anti_tau in antitau_list:
                    fm_sum_0 = fm_tau + fm_anti_tau
                    tollerance = 10

                    if (fm_sum_0.M - 125) < tollerance:
                        fm_sum = fm_sum_0
                    

            #selections
            pt_ditau = np.sqrt(fm_sum.px**2 + fm_sum.py**2)
            eta_ditau = fm_sum.eta

            cut1 = pt_ditau > 40
            cut2 = abs(eta_ditau) < 2.1
            
            if cut1:
                cut1_WL.append(e_id)
            if cut2:
                cut2_WL.append(e_id)

            if cut1 and cut2:
                selected_WL.append(e_id)
            

    #----------------------------------------------------------------------------------      


    print("\nTransverse polarisation\n")

    LHE_file = "unweighted_events_T.lhe"
    events_dict = read_file(LHE_file)

    total_ZT, e_mu_ZT = [], []
    total_WT, e_mu_WT = [], []
    selected_ZT, selected_WT = [], []

    cut1_ZT, cut2_ZT = [], []
    cut1_WT, cut2_WT = [], []


    for e_id, particles_list in events_dict.items():

        if e_id % 10000 == 0:
            print(f"Processing event {e_id}")

        #if e_id > 5000: break

        event_type = f_event_type(particles_list)

        #di-lepton
        if event_type == 23:

            total_ZT.append(e_id)

            e_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            mu_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
           
            for cinematic in particles_list:
                #leptons
                if abs(cinematic.pid) == 11:
                    e_vec = cinematic.build_fm()
                    
                if  abs(cinematic.pid) == 13:
                    mu_vec = cinematic.build_fm()
                  
            if is_zero(e_vec) and is_zero(mu_vec):
                continue

            e_mu_ZT.append(e_id)

            tau_list, antitau_list = [], []
            for cinematic in particles_list:
                if cinematic.pid == 15:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    fm_vec_1 = smeared_momentum
                    tau_list.append(fm_vec_1)

                if cinematic.pid == -15:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    fm_vec_2 = smeared_momentum
                    antitau_list.append(fm_vec_2)
            
            for fm_tau in tau_list:

                for fm_anti_tau in antitau_list:
                    fm_sum_0 = fm_tau + fm_anti_tau
                    tollerance = 10

                    if (fm_sum_0.M - 125) < tollerance:
                        fm_sum = fm_sum_0
                    

            #selections
            pt_ditau = np.sqrt(fm_sum.px**2 + fm_sum.py**2)
            eta_ditau = fm_sum.eta

            cut1 = pt_ditau > 40
            cut2 = abs(eta_ditau) < 2.1
            
            if cut1:
                cut1_ZT.append(e_id)
            if cut2:
                cut2_ZT.append(e_id)

            if cut1 and cut2:
                selected_ZT.append(e_id)


        #di-lepton
        if abs(event_type) == 24:

            total_WT.append(e_id)

            e_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            mu_vec = vector.obj(x = 0, y = 0, z = 0, E = 0)
            
            for cinematic in particles_list:
                #leptons
                if abs(cinematic.pid) == 11:
                    e_vec = cinematic.build_fm()
                    
                if  abs(cinematic.pid) == 13:
                    mu_vec = cinematic.build_fm()
                   
            if is_zero(e_vec) and is_zero(mu_vec):
                continue

            e_mu_WT.append(e_id)

            tau_list, antitau_list = [], []
            for cinematic in particles_list:
                if cinematic.pid == 15:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    fm_vec_1 = smeared_momentum
                    tau_list.append(fm_vec_1)

                if cinematic.pid == -15:
                    non_smeared_momentum = cinematic.build_fm()
                    smeared_momentum = apply_smearing_to(non_smeared_momentum, "jet")
                    fm_vec_2 = smeared_momentum
                    antitau_list.append(fm_vec_2)
            
            for fm_tau in tau_list:

                for fm_anti_tau in antitau_list:
                    fm_sum_0 = fm_tau + fm_anti_tau
                    tollerance = 10

                    if (fm_sum_0.M - 125) < tollerance:
                        fm_sum = fm_sum_0
                    

            #selections
            pt_ditau = np.sqrt(fm_sum.px**2 + fm_sum.py**2)
            eta_ditau = fm_sum.eta

            cut1 = pt_ditau > 40
            cut2 = abs(eta_ditau) < 2.1
            
            if cut1:
                cut1_WT.append(e_id)
            if cut2:
                cut2_WT.append(e_id)

            if cut1 and cut2:
                selected_WT.append(e_id)
            

    #----------------------------------------------------------------------------------      

            
            
    #ZH dataframe
    N_tot_ZL   = len(total_ZL)
    N_em_ZL    = len(e_mu_ZL)
    N_cut1_ZL  = len(cut1_ZL)
    N_cut2_ZL  = len(cut2_ZL)
    N_cut12_ZL = len(selected_ZL)

    N_tot_ZT   = len(total_ZT)
    N_em_ZT    = len(e_mu_ZT)
    N_cut1_ZT  = len(cut1_ZT)
    N_cut2_ZT  = len(cut2_ZT)
    N_cut12_ZT = len(selected_ZT)


    #wh dataframe
    N_tot_WL   = len(total_WL)
    N_em_WL    = len(e_mu_WL)
    N_cut1_WL  = len(cut1_WL)
    N_cut2_WL  = len(cut2_WL)
    N_cut12_WL = len(selected_WL)

    N_tot_WT   = len(total_WT)
    N_em_WT    = len(e_mu_WT)
    N_cut1_WT  = len(cut1_WT)
    N_cut2_WT  = len(cut2_WT)
    N_cut12_WT = len(selected_WT)


    #channel dataframe
    data = [
            ["Total events",            N_tot_ZL,   N_tot_ZT,   N_tot_WL,   N_tot_WT],
            ["Events with e/mu",        N_em_ZL,    N_em_ZT,    N_em_WL,    N_em_WT],
            ["Events after cut1 only",  N_cut1_ZL,  N_cut1_ZT,  N_cut1_WL,  N_cut1_WT],
            ["Events after cut2 only",  N_cut2_ZL,  N_cut2_ZT,  N_cut2_WL,  N_cut2_WT],
            ["Selected events",         N_cut12_ZL, N_cut12_ZT, N_cut12_WL, N_cut12_WT]
            ]


    df_channel = pd.DataFrame(data, columns=["", "ZLH", "ZTH", "WLH", "WTH"])

    print("\n" + "="*60)
    print("{:^60}".format("H > ττ"))
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
    epsilon_tot_WL = N_cut12_WL / N_tot_WL

    #partial efficiencies
    part_epsilon_0_WL = N_em_WL / N_tot_WL
    part_epsilon_1_WL = N_cut1_WL / N_em_WL
    part_epsilon_2_WL = N_cut12_WL / N_cut1_WL


    #WTH
    #total efficiencies
    epsilon_WT = N_tot_WT / 50000
    epsilon_0_WT = N_em_WT / 50000
    epsilon_1_WT = N_cut1_WT / N_em_WT
    epsilon_2_WT = N_cut2_WT / N_em_WT
    epsilon_tot_WT = N_cut12_WT / N_tot_WT

    #partial efficiencies
    part_epsilon_0_WT = N_em_WT / N_tot_WT
    part_epsilon_1_WT = N_cut1_WT / N_em_WT
    part_epsilon_2_WT = N_cut12_WT / N_cut1_WT
    

    efficiencies = np.array([
                            ["channel fraction",        f"{epsilon_ZL:.2f}",    f"{epsilon_ZT:.2f}",      f"{epsilon_WL:.2f}",     f"{epsilon_WT:.2f}"],
                            ["ε₀ = e/μ / total",        f"{epsilon_0_ZL:.2f}",   f"{epsilon_0_ZT:.2f}",   f"{epsilon_0_WL:.2f}",   f"{epsilon_0_WT:.2f}"],
                            ["ε₁ = cut1 / e/μ",         f"{epsilon_1_ZL:.2f}",   f"{epsilon_1_ZT:.2f}",   f"{epsilon_1_WL:.2f}",   f"{epsilon_1_WT:.2f}"],
                            ["ε₂ = cut2 / e/μ",         f"{epsilon_2_ZL:.2f}",   f"{epsilon_2_ZT:.2f}",   f"{epsilon_2_WL:.2f}",   f"{epsilon_2_WT:.2f}"],
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
                                ["ε₁ (cut1)",           f"{part_epsilon_1_ZL:.2f}", f"{epsilon_1_ZL:.2f}", f"{part_epsilon_1_ZT:.2f}", f"{epsilon_1_ZT:.2f}", f"{part_epsilon_1_WL:.2f}", f"{epsilon_1_WL:.2f}", f"{part_epsilon_1_WT:.2f}", f"{epsilon_1_WT:.2f}"],
                                ["ε₂ (cut2)",           f"{part_epsilon_2_ZL:.2f}", f"{epsilon_2_ZL:.2f}", f"{part_epsilon_2_ZT:.2f}", f"{epsilon_2_ZT:.2f}", f"{part_epsilon_2_WL:.2f}", f"{epsilon_2_WL:.2f}", f"{part_epsilon_2_WT:.2f}", f"{epsilon_2_WT:.2f}"],
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
    BR         = 0.06272
    xsection_T = 252  #fb
    lum        = 100            #1/fb 

    #ZLH
    N_tot_L_LHC     = xsection_L     * lum                                      #total longitudinal events number
    N_tot_ZL_LHC    = N_tot_L_LHC    * epsilon_ZL                              #ZLH events number
    N_tot_ZL_ditau  = N_tot_ZL_LHC   * BR                                     #h > bb events number
    N_em_ZL_LHC     = N_tot_ZL_ditau * part_epsilon_0_ZL
    N_cut1_ZL_LHC   = N_em_ZL_LHC    * part_epsilon_1_ZL
    N_cut2_ZL_LHC   = N_cut1_ZL_LHC  * part_epsilon_2_ZL

    kin_eff_ZL_LHC  = N_cut2_ZL_LHC / N_em_ZL_LHC                           #kinematic cuts efficiency
    chan_eff_ZL_LHC = N_cut2_ZL_LHC / N_tot_ZL_ditau                           #channel efficiency
    tot_eff_ZL_LHC  = N_cut2_ZL_LHC / N_tot_ZL_LHC



    # ZTH
    N_tot_T_LHC     = xsection_T     * lum
    N_tot_ZT_LHC    = N_tot_T_LHC    * epsilon_ZT
    N_tot_ZT_ditau  = N_tot_ZT_LHC   * BR
    N_em_ZT_LHC     = N_tot_ZT_ditau * part_epsilon_0_ZT
    N_cut1_ZT_LHC   = N_em_ZT_LHC    * part_epsilon_1_ZT
    N_cut2_ZT_LHC   = N_cut1_ZT_LHC  * part_epsilon_2_ZT

    kin_eff_ZT_LHC  = N_cut2_ZT_LHC / N_em_ZT_LHC
    chan_eff_ZT_LHC = N_cut2_ZT_LHC / N_tot_ZT_ditau                           #channel efficiency
    tot_eff_ZT_LHC  = N_cut2_ZT_LHC / N_tot_ZT_LHC



    # WLH
    N_tot_L_LHC     = xsection_L     * lum
    N_tot_WL_LHC    = N_tot_L_LHC    * epsilon_WL
    N_tot_WL_ditau  = N_tot_WL_LHC   * BR
    N_em_WL_LHC     = N_tot_WL_ditau * part_epsilon_0_WL
    N_cut1_WL_LHC   = N_em_WL_LHC    * part_epsilon_1_WL
    N_cut2_WL_LHC   = N_cut1_WL_LHC  * part_epsilon_2_WL
    
    kin_eff_WL_LHC  = N_cut2_WL_LHC / N_em_WL_LHC
    chan_eff_WL_LHC = N_cut2_WL_LHC / N_tot_WL_ditau                           #channel efficiency
    tot_eff_WL_LHC  = N_cut2_WL_LHC / N_tot_WL_LHC

    # WTH
    N_tot_T_LHC     = xsection_T     * lum
    N_tot_WT_LHC    = N_tot_T_LHC    * epsilon_WT
    N_tot_WT_ditau  = N_tot_WT_LHC   * BR
    N_em_WT_LHC     = N_tot_WT_ditau * part_epsilon_0_WT
    N_cut1_WT_LHC   = N_em_WT_LHC    * part_epsilon_1_WT
    N_cut2_WT_LHC   = N_cut1_WT_LHC  * part_epsilon_2_WT

    kin_eff_WT_LHC  = N_cut2_WT_LHC / N_em_WT_LHC
    chan_eff_WT_LHC = N_cut2_WT_LHC / N_tot_WT_ditau                           #channel efficiency
    tot_eff_WT_LHC  = N_cut2_WT_LHC / N_tot_WT_LHC


    data = [
            ["Total L/T events",             f"{N_tot_L_LHC:.0f}",      f"{N_tot_T_LHC:.0f}",       f"{N_tot_L_LHC:.0f}",       f"{N_tot_T_LHC:.0f}"],
            ["Total V{i}H events",           f"{N_tot_ZL_LHC:.0f}",     f"{N_tot_ZT_LHC:.0f}",      f"{N_tot_WL_LHC:.0f}",      f"{N_tot_WT_LHC:.0f}"],
            ["Total bb events",              f"{N_tot_ZL_ditau:.0f}",   f"{N_tot_ZT_ditau:.0f}",    f"{N_tot_WL_ditau:.0f}",    f"{N_tot_WT_ditau:.0f}"],
            ["Events with e/mu",             f"{N_em_ZL_LHC:.0f}",      f"{N_em_ZT_LHC:.0f}",       f"{N_em_WL_LHC:.0f}",       f"{N_em_WT_LHC:.0f}"],
            ["Events after cut1",            f"{N_cut1_ZL_LHC:.0f}",    f"{N_cut1_ZT_LHC:.0f}",     f"{N_cut1_WL_LHC:.0f}",     f"{N_cut1_WT_LHC:.0f}"],
            ["Events after cut1, cut2",      f"{N_cut2_ZL_LHC:.0f}",    f"{N_cut2_ZT_LHC:.0f}",     f"{N_cut2_WL_LHC:.0f}",     f"{N_cut2_WT_LHC:.0f}"],
            ["Kinematic cuts efficiency",    f"{kin_eff_ZL_LHC:.2f}",   f"{kin_eff_ZT_LHC:.2f}",    f"{kin_eff_WL_LHC:.2f}",    f"{kin_eff_WT_LHC:.2f}"],
            ["Channel efficiency",           f"{chan_eff_ZL_LHC:.2f}",  f"{chan_eff_ZT_LHC:.2f}",   f"{chan_eff_WL_LHC:.2f}",   f"{chan_eff_WT_LHC:.2f}"],
            ["Total efficiency",             f"{tot_eff_ZL_LHC:.4f}",   f"{tot_eff_ZT_LHC:.4f}",    f"{tot_eff_WL_LHC:.4f}",    f"{tot_eff_WT_LHC:.4f}"]
        ]


    df_channel = pd.DataFrame(data, columns=["", "ZLH", "ZTH", "WLH", "WTH"])

    print("\n" + "="*60)
    print("{:^60}".format("H > ττ (LHC)"))
    print("="*60 + "\n")
    print(tabulate(df_channel.values.tolist(), headers=df_channel.columns.tolist(), tablefmt="grid"))


    #----------------------------------------------------------------------------------------------------------------------------------------------------------------


    print("Ordering cuts")

    cut_list = ["cut1", "cut2"]

    for col in column[1:]:
        
        #re-arranging
        values = df_efficiencies.loc[1:2, col].copy()
        num_values = [float(v) for v in values if v != "-"]

        ordered_val = sorted(num_values)
        print("Ordered total efficiencies:", ordered_val)
        sorted_indices = sorted(range(len(num_values)), key=lambda i: num_values[i])

        ordered_cut_list = [cut_list[i] for i in sorted_indices]
        print(f"Ordered cuts for {col}:", ordered_cut_list)

        j = 0
        for m in range(1, 3): 

            if df_efficiencies.at[m, col] == "-":
                continue

            df_efficiencies.at[m, col] = f"{ordered_val[j]:.2f}"
            j += 1

    
    #ordering cumulative efficiencies
    #ZLH
    part_epsilon_1_ZL = N_cut2_ZL / N_em_ZL
    part_epsilon_2_ZL = N_cut12_ZL / N_cut2_ZL

    #ZTH
    part_epsilon_1_ZT = N_cut2_ZT / N_em_ZT
    part_epsilon_2_ZT = N_cut12_ZT / N_cut2_ZT

    #WLH
    part_epsilon_1_WL = N_cut2_WL / N_em_WL
    part_epsilon_2_WL = N_cut12_WL / N_cut2_WL

    #WTH
    part_epsilon_1_WT = N_cut2_WT / N_em_WT
    part_epsilon_2_WT = N_cut12_WT / N_cut2_WT

    

    part_total_eff_ordered = np.array([
                                ["ε₀ = e/μ / total",    f"{epsilon_0_ZL:.2f}",      f"{epsilon_0_ZL:.2f}", f"{epsilon_0_ZT:.2f}",      f"{epsilon_0_ZT:.2f}", f"{epsilon_0_WL:.2f}",      f"{epsilon_0_WL:.2f}", f"{epsilon_0_WT:.2f}",      f"{epsilon_0_WT:.2f}"],
                                ["ε₁ (cut1)",           f"{part_epsilon_1_ZL:.2f}", f"{epsilon_2_ZL:.2f}", f"{part_epsilon_1_ZT:.2f}", f"{epsilon_2_ZT:.2f}", f"{part_epsilon_1_WL:.2f}", f"{epsilon_2_WL:.2f}", f"{part_epsilon_1_WT:.2f}", f"{epsilon_2_WT:.2f}"],
                                ["ε₂ (cut2)",           f"{part_epsilon_2_ZL:.2f}", f"{epsilon_1_ZL:.2f}", f"{part_epsilon_2_ZT:.2f}", f"{epsilon_1_ZT:.2f}", f"{part_epsilon_2_WL:.2f}", f"{epsilon_1_WL:.2f}", f"{part_epsilon_2_WT:.2f}", f"{epsilon_1_WT:.2f}"],
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