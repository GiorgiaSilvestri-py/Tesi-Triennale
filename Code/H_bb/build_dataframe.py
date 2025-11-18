from array import array
from utils_nonlhe import *
import vector

def main():

    print("Longitudinal polarisation")

    LHE_file = "unweighted_events_L.lhe"
    events_dict = read_file(LHE_file)

    #TTree ZLH
    tree_ZLH = ROOT.TTree("Events", "ZLH channel")

    pt_Z_ZLH          = array('f', [0.])
    deltaphi_zjj_ZLH  = array('f', [0.])
    event_type_ZLH    = ROOT.std.string()

    tree_ZLH.Branch("pt_z", pt_Z_ZLH, "pt_z/F")
    tree_ZLH.Branch("deltaphi_zjj", deltaphi_zjj_ZLH, "deltaphi_zjj/F")
    tree_ZLH.Branch("event_type", event_type_ZLH)

    # TTree WLH
    tree_WLH = ROOT.TTree("Events", "WLH channel")

    pt_W_WLH          = array('f', [0.])
    pt_j1_WLH         = array('f', [0.])
    pt_j2_WLH         = array('f', [0.])
    pt_jj_WLH         = array('f', [0.])
    deltaphi_Wjj_WLH  = array('f', [0.])
    deltaphi_lmet_WLH = array('f', [0.])
    event_type_WLH    = ROOT.std.string()

    tree_WLH.Branch("pt_w", pt_W_WLH, "pt_w/F")
    tree_WLH.Branch("pt_j1", pt_j1_WLH, "pt_j1/F")
    tree_WLH.Branch("pt_j2", pt_j2_WLH, "pt_j2/F")
    tree_WLH.Branch("pt_jj", pt_jj_WLH, "pt_jj/F")
    tree_WLH.Branch("deltaphi_Wjj", deltaphi_Wjj_WLH, "deltaphi_Wjj/F")
    tree_WLH.Branch("deltaphi_l_met", deltaphi_lmet_WLH, "deltaphi_l_met/F")
    tree_WLH.Branch("event_type", event_type_WLH)

    #--------------------------------------------------------------------------

    #events loop

    for e_id, part_dict in events_dict.items():

        if e_id % 10000 == 0:
            print(f"Processing event {e_id}")

        if e_id > 1000:    break

        event_type = f_event_type(part_dict)

        #di-lepton
        if event_type == 23:

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
                    
            if is_nonzero(fm_vec_1) or is_nonzero(fm_vec_2):
                event_type_ZLH.replace(0, len(event_type_ZLH), "e/mu")
            else:
                event_type_ZLH.replace(0, len(event_type_ZLH), "tau/neutrino")

                
            fm_Z = fm_vec_1 + fm_vec_2

            
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
            pt_V      = np.sqrt(fm_Z.px**2 + fm_Z.py**2)
            delta_phi = fm_Z.deltaphi(fm_jj)

            pt_Z_ZLH[0]         = pt_V
            deltaphi_zjj_ZLH[0] = delta_phi

            tree_ZLH.Fill()
             

        #single-lepton
        if event_type == 24 or event_type == -24:

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
                    
            if is_nonzero(fm_vec_1):
                event_type_WLH.replace(0, len(event_type_WLH), "e/mu")
            else:
                event_type_WLH.replace(0, len(event_type_WLH), "tau")

            fm_W = fm_vec_1 + fm_vec_2

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
            MET               = vector.obj(x = fm_vec_2.px, y =  fm_vec_2.py, z = 0)
            pt_W              = np.sqrt(fm_W.px**2 + fm_W.py**2)
            pt_j1             = np.sqrt(jet1_fm.px**2 + jet1_fm.py**2)
            pt_j2             = np.sqrt(jet2_fm.px**2 + jet2_fm.py**2)
            pt_jj             = np.sqrt(fm_jj.px**2 + fm_jj.py**2)            
            delta_phi_lep_met = fm_vec_1.deltaphi(MET) 
            delta_phi_w_jj    = fm_W.deltaphi(fm_jj)

        
            pt_W_WLH[0]           = pt_W
            pt_j1_WLH[0]          = pt_j1
            pt_j2_WLH[0]          = pt_j2
            pt_jj_WLH[0]          = pt_jj
            deltaphi_Wjj_WLH[0]   = delta_phi_lep_met
            deltaphi_lmet_WLH[0]  = delta_phi_w_jj

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

    pt_Z_ZTH           = array('f', [0.])
    deltaphi_zjj_ZTH   = array('f', [0.])
    event_type_ZTH     = ROOT.std.string()

    tree_ZTH.Branch("pt_z",         pt_Z_ZTH, "pt_z/F")
    tree_ZTH.Branch("deltaphi_zjj", deltaphi_zjj_ZTH, "deltaphi_zjj/F")
    tree_ZTH.Branch("event_type",   event_type_ZTH)

    #TTree initialisation - WTH
    tree_WTH = ROOT.TTree("Events", "WTH channel")

    pt_W_WTH          = array('f', [0.])
    pt_j1_WTH         = array('f', [0.])
    pt_j2_WTH         = array('f', [0.])
    pt_jj_WTH         = array('f', [0.])
    deltaphi_Wjj_WTH  = array('f', [0.])
    deltaphi_lmet_WTH = array('f', [0.])
    event_type_WTH    = ROOT.std.string()

    tree_WTH.Branch("pt_w",           pt_W_WTH,          "pt_w/F")
    tree_WTH.Branch("pt_j1",          pt_j1_WTH,         "pt_j1/F")
    tree_WTH.Branch("pt_j2",          pt_j2_WTH,         "pt_j2/F")
    tree_WTH.Branch("pt_jj",          pt_jj_WTH,         "pt_jj/F")
    tree_WTH.Branch("deltaphi_Wjj",   deltaphi_Wjj_WTH,  "deltaphi_Wjj/F")
    tree_WTH.Branch("deltaphi_l_met", deltaphi_lmet_WTH, "deltaphi_l_met/F")
    tree_WTH.Branch("event_type",     event_type_WTH)

    
    #--------------------------------------------------------------------------------

    #events loop
    for e_id, part_dict in events_dict.items():

        if e_id % 10000 == 0:
            print(f"Processing event {e_id}")

        if e_id > 1000:    break

        event_type = f_event_type(part_dict)

        #di-lepton
        if event_type == 23:

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
                    
            if not is_zero(fm_vec_1) or not is_zero(fm_vec_2):
                event_type_ZTH.replace(0, len(event_type_ZTH), "e/mu")
            else:
                event_type_ZTH.replace(0, len(event_type_ZTH), "tau/neutrino")
                
            fm_Z = fm_vec_1 + fm_vec_2

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
                    
            fm_jj = jet1_fm + jet2_fm
            

            #selections
            pt_V = np.sqrt(fm_Z.px**2 + fm_Z.py**2)
            delta_phi = fm_Z.deltaphi(fm_jj)

            pt_Z_ZTH[0]         = pt_V
            deltaphi_zjj_ZTH[0] = delta_phi

            tree_ZTH.Fill()

            
        #single-lepton
        if event_type == 24 or event_type == -24:

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
                    
            if not is_zero(fm_vec_1) or not is_zero(fm_vec_2):          #non-null vec = 11/13
                event_type_WTH.replace(0, len(event_type_WTH), "e/mu")
            else:                                                       #null vec = 15
                event_type_WTH.replace(0, len(event_type_WTH), "tau")

            fm_W = fm_vec_1 + fm_vec_2

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
            MET               = vector.obj(x = fm_vec_2.px, y =  fm_vec_2.py, z = 0)
            pt_W              = np.sqrt(fm_W.px**2 + fm_W.py**2)
            pt_j1             = np.sqrt(jet1_fm.px**2 + jet1_fm.py**2)
            pt_j2             = np.sqrt(jet2_fm.px**2 + jet2_fm.py**2)
            pt_jj             = np.sqrt(fm_jj.px**2 + fm_jj.py**2)            
            delta_phi_lep_met = fm_vec_1.deltaphi(MET) 
            delta_phi_w_jj    = fm_W.deltaphi(fm_jj)

            pt_W_WTH[0]           = pt_W
            pt_j1_WTH[0]          = pt_j1
            pt_j2_WTH[0]          = pt_j2
            pt_jj_WTH[0]          = pt_jj
            deltaphi_Wjj_WTH[0]   = delta_phi_lep_met
            deltaphi_lmet_WTH[0]  = delta_phi_w_jj

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