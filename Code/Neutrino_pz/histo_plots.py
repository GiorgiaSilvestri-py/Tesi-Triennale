import matplotlib.pyplot as plt
import matplotlib
import math
from matplotlib.colors import LogNorm
import numpy as np
import vector
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import os
import ROOT
import seaborn as sb
from scipy.stats import binned_statistic_2d
from utils_nonlhe import solve_equation, read_file, f_event_type, get_fm_of, get_thetastar_of


def compute_neutrino_momentum(events_dict, e_id):
   
    vec_lepton            = vector.obj(px = 0, py = 0, pz = 0, E = 0)
    vec_neutrino_expected = vector.obj(px = 0, py = 0, pz = 0, E = 0)
    
    p_dict = events_dict[e_id]
    for pid, cinematic in p_dict.items():
        if abs(pid) in [11, 13, 15]:
            vec_lepton = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
            
        if abs(pid) in [12, 14, 16]:
            vec_neutrino_expected = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)

    #skip event if either one is None
    if not is_nonzero(vec_lepton) or not is_nonzero(vec_neutrino_expected):
        return
        
  
    return compute_neutrino_momentum_from_particles(vec_lepton, vec_neutrino_expected)
    
def compute_neutrino_momentum_diff(events_dict, e_id):
   
    vec_lepton            = vector.obj(px = 0, py = 0, pz = 0, E = 0)
    vec_neutrino_expected = vector.obj(px = 0, py = 0, pz = 0, E = 0)
    
    p_dict = events_dict[e_id]
    for pid, cinematic in p_dict.items():
        if abs(pid) in [11, 13, 15]:
            vec_lepton = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
            
        if abs(pid) in [12, 14, 16]:
            vec_neutrino_expected = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)

    #skip event if either one is None
    if not is_nonzero(vec_lepton) or not is_nonzero(vec_neutrino_expected):
        return
        
  
    return compute_neutrino_momentum_from_particles_diff(vec_lepton, vec_neutrino_expected)
    

def compute_neutrino_momentum_from_particles(vec_lepton, vec_neutrino_expected):
    
    #compute pz with maximum angle 3d
    
    expected_pz_neutrino = vec_neutrino_expected.pz
    
    #solvin second grade equation
    M = 80
    pt_lepton_dot_pt_neutrino = vec_lepton.px * vec_neutrino_expected.px + vec_lepton.py * vec_neutrino_expected.py
    pt_neutrino = vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2
    
    alpha = M**2 + 2 * pt_lepton_dot_pt_neutrino
    a = 4 * (vec_lepton.pz ** 2) - 4 * vec_lepton.E ** 2
    b = 4 * alpha * vec_lepton.pz
    c = alpha**2 - 4 * vec_lepton.E**2 * pt_neutrino
    
    #second grade equation solutions
    pz_neutrino_delta_plus, pz_neutrino_delta_minus = solve_equation(a, b, c)
    

    #computing angles between lepton and neutrino
    neutrino_energy_deltaplus = np.sqrt( vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2 + pz_neutrino_delta_plus**2 )
    built_neutrino_vector_plus = vector.obj(
                                            px = vec_neutrino_expected.px,
                                            py = vec_neutrino_expected.py, 
                                            pz = pz_neutrino_delta_plus,
                                            E  = neutrino_energy_deltaplus
                                            )
  
    neutrino_direction_plus = built_neutrino_vector_plus.to_3D()
    lepton_direction        = vec_lepton.to_3D()
    
    
    #first solution
    cos_theta_plus = compute_angle_between(lepton_direction, neutrino_direction_plus)
    
    neutrino_energy_deltaminus = np.sqrt( vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2 + pz_neutrino_delta_minus**2 )
    built_neutrino_vector_minus = vector.obj(
                                            px = vec_neutrino_expected.px,
                                            py = vec_neutrino_expected.py, 
                                            pz = pz_neutrino_delta_minus,
                                            E  = neutrino_energy_deltaminus
                                            )
            
    neutrino_direction_minus = built_neutrino_vector_minus.to_3D()
    #second solution
    cos_theta_minus = compute_angle_between(lepton_direction, neutrino_direction_minus)         #coseno
    
    
    #selecting correct solution
    if (cos_theta_plus) > (cos_theta_minus):
        built_neutrino_vector = built_neutrino_vector_plus
        #pz_neutrino_built = built_neutrino_vector_plus.pz
        cos = cos_theta_plus
        
        
    else:
        built_neutrino_vector = built_neutrino_vector_minus
        #pz_neutrino_built = built_neutrino_vector_minus.pz
        cos = cos_theta_minus
   
        
    return built_neutrino_vector, cos


def compute_neutrino_momentum_from_particles_diff(vec_lepton, vec_neutrino_expected):
   
    #compute pz cosidering lepton pz   
   
    expected_pz_neutrino = vec_neutrino_expected.pz
    
    #solvin second grade equation
    M = 80
    pt_lepton_dot_pt_neutrino = vec_lepton.px * vec_neutrino_expected.px + vec_lepton.py * vec_neutrino_expected.py
    pt_neutrino = vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2
    
    alpha = M**2 + 2 * pt_lepton_dot_pt_neutrino
    a = 4 * (vec_lepton.pz ** 2) - 4 * vec_lepton.E ** 2
    b = 4 * alpha * vec_lepton.pz
    c = alpha**2 - 4 * vec_lepton.E**2 * pt_neutrino
    
    #second grade equation solutions
    pz_neutrino_delta_plus, pz_neutrino_delta_minus = solve_equation(a, b, c)
    
    neutrino_energy_deltaplus = np.sqrt( vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2 + pz_neutrino_delta_plus**2 )
    built_neutrino_vector_plus = vector.obj(
                                            px = vec_neutrino_expected.px,
                                            py = vec_neutrino_expected.py, 
                                            pz = pz_neutrino_delta_plus,
                                            E  = neutrino_energy_deltaplus
                                            )
  
    neutrino_direction_plus = built_neutrino_vector_plus.to_3D()
    lepton_direction        = vec_lepton.to_3D()
    
    
  
    neutrino_energy_deltaminus = np.sqrt( vec_neutrino_expected.px**2 + vec_neutrino_expected.py**2 + pz_neutrino_delta_minus**2 )
    built_neutrino_vector_minus = vector.obj(
                                            px = vec_neutrino_expected.px,
                                            py = vec_neutrino_expected.py, 
                                            pz = pz_neutrino_delta_minus,
                                            E  = neutrino_energy_deltaminus
                                            )
            
    
    #computing differences with lepton pz
    diff_plus = abs(pz_neutrino_delta_plus - vec_lepton.pz)
    diff_minus = abs(pz_neutrino_delta_minus - vec_lepton.pz)
   
    #selecting correct solution

    if diff_plus < diff_minus:
        built_neutrino_vector = built_neutrino_vector_plus
        
    else:
        built_neutrino_vector = built_neutrino_vector_minus
        
    return built_neutrino_vector

def compute_angle_between(v1, v2):
    """
    Calcola il coseno dell'angolo tra due vettori tridimensionali
    """
    
    p1 = np.array([v1.x, v1.y, v1.z])
    p2 = np.array([v2.x, v2.y, v2.z])

    cos_theta = np.dot(p1, p2) / (np.linalg.norm(p1) * np.linalg.norm(p2))
    
    return cos_theta
    
def main():


    LHE_L = "unweighted_events_L.lhe"
    eventsL_dict = read_file(LHE_L)
    
    #first method - minimum angle
    vec_lepton = {}
    vec_neutrino = {}
    w_exp = {}
    fm_h_dict = {}
    fm_w_dict = {}
    vec_built = {}
    cos_list_l = []
    
    for evt_id, p_dict in eventsL_dict.items():
        e_type = f_event_type(p_dict)
        
        if e_type == 23:
            continue
        for pid, cinematic in p_dict.items():
            if abs(pid) in [12, 14, 16]:
                vec_neutrino[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
        for pid, cinematic in p_dict.items():     
            if abs(pid) in [11, 13, 15]:
                vec_lepton[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
                vec_built[evt_id], cos = compute_neutrino_momentum_from_particles(vec_lepton[evt_id], vec_neutrino[evt_id])
                cos_list_l.append(cos)
                
            if abs(pid) in [24]:
                w_exp[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
        for pid, cinematic in p_dict.items():
            if abs(pid) in [25]:
                fm_h_dict[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
            
        fm_w_dict[evt_id] = vec_lepton[evt_id] + vec_built[evt_id]
        
    #other variables
    fm_pz_list_l, fm_pt_list_l = [], []
    decay_angle_list, decay_angle_list_built = [], []
    decay_cos_list, decay_cos_list_built = [], []
    pz_diff, pz_diff_built_angle_l = [], []
    pt_w_l = []
    
    R_l_angle = []
    
    for e_id in w_exp.keys():
        fm_w = w_exp[e_id]
        fm_n = vec_neutrino[e_id]
        fm_n_built = vec_built[e_id]
        fm_l = vec_lepton[e_id]
        
        fm_pz_list_l.append(abs(fm_l.pz))
        fm_pt_list_l.append(abs(np.sqrt(fm_l.px**2 + fm_l.py**2)))
        pz_diff.append(abs(fm_n.pz - fm_l.pz))
        pz_diff_built_angle_l.append(abs(fm_n_built.pz - fm_l.pz))
        pt_w_l.append(abs(np.sqrt(fm_w.px**2 + fm_w.py**2)))
        
        R_l_angle.append(fm_n_built.pz / fm_n.pz)
        
    
    #ii method - pz difference
    
    vec_lepton = {}
    vec_neutrino = {}
    w_exp = {}
    fm_h_dict = {}
    vec_built = {}
    cos_list_l_diff = []
    
    for evt_id, p_dict in eventsL_dict.items():
        e_type = f_event_type(p_dict)
        
        if e_type == 23:
            continue
        for pid, cinematic in p_dict.items():
            if abs(pid) in [12, 14, 16]:
                vec_neutrino[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
        for pid, cinematic in p_dict.items():     
            if abs(pid) in [11, 13, 15]:
                vec_lepton[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
                vec_built[evt_id] = compute_neutrino_momentum_from_particles_diff(vec_lepton[evt_id], vec_neutrino[evt_id])
                
            if abs(pid) in [24]:
                w_exp[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
                
        for pid, cinematic in p_dict.items():
            if abs(pid) in [25]:
                fm_h_dict[evt_id] = vector.obj(px = cinematic.px, py = cinematic.py, pz = cinematic.pz, E = cinematic.E)
        
        fm_w_dict[evt_id] = vec_lepton[evt_id] + vec_built[evt_id]
                
                       
   
    #other variables
    pz_diff, pz_diff_built = [], []
    fm_pt_list_l = []
    R_l_diff = []
    pt_h, pz_h = [], []
    eta_lep, eta_H = [], []
    
    for e_id in w_exp.keys():
        fm_w = w_exp[e_id]
        fm_h = fm_h_dict[e_id]
        fm_n = vec_neutrino[e_id]
        fm_n_built = vec_built[e_id]
        fm_l = vec_lepton[e_id]
                
        pz_diff.append(abs(fm_n.pz - fm_l.pz))
        pz_diff_built.append(abs(fm_n_built.pz - fm_l.pz))
        
        fm_pt_list_l.append(abs(np.sqrt(fm_l.px**2 + fm_l.py**2)))
        
        eta_lep.append(fm_l.eta)
        eta_H.append(fm_h.eta)
        
        pz_h.append(fm_h.pz)
        pt_h.append(abs(np.sqrt(fm_h.px**2 + fm_h.py**2)))
        
        
        neutrino_direction = fm_n.to_3D()
        n_built_direction = fm_n_built.to_3D()
        lepton_direction   = fm_l.to_3D()
        
        cos_list_l_diff.append(compute_angle_between(n_built_direction, lepton_direction))
        
        R_l_diff.append(fm_n_built.pz / fm_n.pz)
        
        

    #Pz
    fm_pz_array = np.array(fm_pz_list_l)
    fm_pz_h_array = np.array(pz_h)
    
    #Pt
    fm_pt_array = np.array(fm_pt_list_l)
    fm_pt_w_array = np.array(pt_w_l)
    fm_pt_h_array = np.array(pt_h)
    
    #eta
    eta_lep_array = np.array(eta_lep)
    eta_h_array = np.array(eta_H)
    
    #cos_array = np.array(cos_list_l)
    R_array_angle = np.array(R_l_angle)
    R_array_diff  = np.array(R_l_diff)
    
    
    #------------------------------------------------------------------------------------------------------------------------------------
    
    #HISTO1D
    
    h1 = ROOT.TProfile("Profile 1", "Minimum Angle",
                           40, 0, 1400,
                           -2, 10)
    for xi, yi in zip(fm_pz_array, R_array_angle):
        h1.Fill(xi, yi)
        
    h2 = ROOT.TProfile("Profile 2", "Pz difference",
                           40, 0, 1400,
                           -2, 10)
    for xi, yi in zip(fm_pz_array, R_array_diff):
        h2.Fill(xi, yi)
        
    N = 50
    stops_raw = np.logspace(-2, 0, N, base=10)
    stops = (stops_raw - stops_raw.min()) / (stops_raw.max() - stops_raw.min())

    cmap = matplotlib.colormaps.get_cmap('inferno')
    red   = np.array([cmap(x)[0] for x in stops], dtype=np.float64)
    green = np.array([cmap(x)[1] for x in stops], dtype=np.float64)
    blue  = np.array([cmap(x)[2] for x in stops], dtype=np.float64)

    ROOT.TColor.CreateGradientColorTable(N, stops, red, green, blue, 255)
    ROOT.gStyle.SetNumberContours(255)


    ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetPalette(ROOT.kRainbow)

    canvas = ROOT.TCanvas("canvas", "Confronto Profile1D", 1600, 1000)
    canvas.Divide(2, 1)

    canvas.cd(1)
    h1.Draw("COLZ")
    
    h1.GetXaxis().SetTitle("pz lepton (GeV)")
    h1.GetYaxis().SetTitle("Ratio")
    
    h2.GetXaxis().SetTitle("pz lepton (GeV)")
    h2.GetYaxis().SetTitle("Ratio")

    canvas.cd(2)
    h2.Draw("COLZ")

    canvas.SaveAs("pzlepton_profile1d.png")
    
    
    
    #------------------------------------------------------------------------------------------------------------------------------------
    
    #HISTO1D
    
    h1 = ROOT.TProfile("Profile 1", "Minimum Angle",
                           40, 0, 400,
                           -2, 10)
    for xi, yi in zip(fm_pt_array, R_array_angle):
        h1.Fill(xi, yi)
        
    h2 = ROOT.TProfile("Profile 2", "Pz difference",
                           40, 0, 400,
                           -2, 10)
    for xi, yi in zip(fm_pt_array, R_array_diff):
        h2.Fill(xi, yi)
        
    N = 50
    stops_raw = np.logspace(-2, 0, N, base=10)
    stops = (stops_raw - stops_raw.min()) / (stops_raw.max() - stops_raw.min())

    cmap = matplotlib.colormaps.get_cmap('inferno')
    red   = np.array([cmap(x)[0] for x in stops], dtype=np.float64)
    green = np.array([cmap(x)[1] for x in stops], dtype=np.float64)
    blue  = np.array([cmap(x)[2] for x in stops], dtype=np.float64)

    ROOT.TColor.CreateGradientColorTable(N, stops, red, green, blue, 255)
    ROOT.gStyle.SetNumberContours(255)


    ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetPalette(ROOT.kRainbow)

    canvas = ROOT.TCanvas("canvas", "Confronto Profile1D", 1600, 1000)
    canvas.Divide(2, 1)

    canvas.cd(1)
    h1.Draw("COLZ")
    
    h1.GetXaxis().SetTitle("pt lepton (GeV)")
    h1.GetYaxis().SetTitle("Ratio")
    
    h2.GetXaxis().SetTitle("pt lepton (GeV)")
    h2.GetYaxis().SetTitle("Ratio")

    canvas.cd(2)
    h2.Draw("COLZ")

    canvas.SaveAs("ptlepton_profile1d.png")


    #HISTO1D
    
    h1 = ROOT.TProfile("Profile 1", "Minimum Angle",
                           40, -6, 6,
                           -2, 10)
    for xi, yi in zip(eta_lep_array, R_array_angle):
        h1.Fill(xi, yi)
        
    h2 = ROOT.TProfile("Profile 2", "Pz difference",
                           40, -6, 6,
                           -2, 10)
    for xi, yi in zip(eta_lep_array, R_array_diff):
        h2.Fill(xi, yi)
        
    N = 50
    stops_raw = np.logspace(-2, 0, N, base=10)
    stops = (stops_raw - stops_raw.min()) / (stops_raw.max() - stops_raw.min())

    cmap = matplotlib.colormaps.get_cmap('inferno')
    red   = np.array([cmap(x)[0] for x in stops], dtype=np.float64)
    green = np.array([cmap(x)[1] for x in stops], dtype=np.float64)
    blue  = np.array([cmap(x)[2] for x in stops], dtype=np.float64)

    ROOT.TColor.CreateGradientColorTable(N, stops, red, green, blue, 255)
    ROOT.gStyle.SetNumberContours(255)


    ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetPalette(ROOT.kRainbow)

    canvas = ROOT.TCanvas("canvas", "Confronto Profile1D", 1600, 1000)
    canvas.Divide(2, 1)

    canvas.cd(1)
    h1.Draw("COLZ")
    
    h1.GetXaxis().SetTitle("eta lepton (rad)")
    h1.GetYaxis().SetTitle("Ratio")
    
    h2.GetXaxis().SetTitle("eta lepton (rad)")
    h2.GetYaxis().SetTitle("Ratio")

    canvas.cd(2)
    h2.Draw("COLZ")

    canvas.SaveAs("etalepton_profile1d.png")


    
    #------------------------------------------------------------------------------------------------------------
    
    #PROFILE3D
        
    
    h1 = ROOT.TProfile3D("Profile 1", "Minimum Angle",
                           20, 0, 1400,
                           20, -6, 6,
                           25, 0, 400,)
                           
    h1.SetMinimum(-2)
    h1.SetMaximum(10)

    for xi, yi, zi, ti in zip(fm_pz_array, eta_lep, fm_pt_array,  R_array_angle):
        h1.Fill(xi, yi, zi, ti)
        
    h2 = ROOT.TProfile3D("Profile 2", "Pz difference",
                           20, 0, 1400,
                           20, -6, 6,
                           25, 0, 400,)
                           
    h2.SetMinimum(-2)
    h2.SetMaximum(10)

                           
                       
    for xi, yi, zi, ti in zip(fm_pz_array, eta_lep, fm_pt_array, R_array_diff):
        h2.Fill(xi, yi, zi, ti)
        
    N = 50
    stops_raw = np.logspace(-2, 0, N, base=10)
    stops = (stops_raw - stops_raw.min()) / (stops_raw.max() - stops_raw.min())

    cmap = matplotlib.colormaps.get_cmap('inferno')
    red   = np.array([cmap(x)[0] for x in stops], dtype=np.float64)
    green = np.array([cmap(x)[1] for x in stops], dtype=np.float64)
    blue  = np.array([cmap(x)[2] for x in stops], dtype=np.float64)

    ROOT.TColor.CreateGradientColorTable(N, stops, red, green, blue, 255)
    ROOT.gStyle.SetNumberContours(255)
        
    ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetPalette(ROOT.kRainbow)

    canvas = ROOT.TCanvas("canvas", "Confronto Profile2D", 1900, 1000)
    canvas.Divide(2, 1)

    canvas.cd(1)
    h1.Draw("COLZ")
    
    h1.GetXaxis().SetTitle("pz lepton (GeV)")
    h1.GetZaxis().SetTitle("pt lepton (GeV)")
    h1.GetYaxis().SetTitle("eta lepton (rad)")
    
    h2.GetXaxis().SetTitle("pz lepton (GeV)")
    h2.GetZaxis().SetTitle("pt lepton (GeV)")
    h2.GetYaxis().SetTitle("eta lepton (rad)")
    
    canvas.cd(2)
    h2.Draw("COLZ")

    canvas.SaveAs("profile3d.png")
    
    output_file = ROOT.TFile("profile3d_output.root", "RECREATE")
    h1.Write()
    h2.Write()
    canvas.Write()
    output_file.Close()
    

    #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    #HISTOGRAMS - root
    
    #Pz
    fm_pz_array = np.array(fm_pz_list_l)
    fm_pz_h_array = np.array(pz_h)
    
    #Pt
    fm_pt_array = np.array(fm_pt_list_l)
    fm_pt_w_array = np.array(pt_w_l)
    fm_pt_h_array = np.array(pt_h)
    
    #eta
    eta_lep_array = np.array(eta_lep)
    eta_h_array = np.array(eta_H)
    
    #cos_array = np.array(cos_list_l)
    R_array_angle = np.array(R_l_angle)
    R_array_diff  = np.array(R_l_diff)
    
    array_list = [fm_pz_array, fm_pz_h_array, fm_pt_array, fm_pt_w_array, fm_pt_h_array, eta_lep_array, eta_h_array]
    var_name = ["pz lepton (GeV)", "pz higgs (GeV)", "pt lepton (GeV)", "pt W (GeV)", "pt higgs (GeV)", "eta lepton", "eta higgs"]
    
    for (a_list1, var1) in zip(array_list, var_name):
        x_bins = 20
        y_bins = 25
        xmin = np.min(a_list1)
        xmax = np.max(a_list1)
        
        
        
        for (a_list2, var2) in zip(array_list, var_name):
            if var1 == var2:
                continue
                
            output_dir = "/home/giorgia/MG/Risultati/Confronto_L-T/Pz_neutrino/histo_root"
            os.makedirs(output_dir, exist_ok=True)
            filename_inverted = os.path.join(output_dir, f"{var2}_{var1}_histo2d_median_.png")

            if os.path.exists(filename_inverted):
                print(f"{filename_inverted} already created")
                continue
                
            ymin = np.min(a_list2)
            ymax = np.max(a_list2)
                
            h1 = ROOT.TProfile2D("Profile 1", "Minimum Angle",
                               x_bins, xmin, xmax,
                               y_bins, ymin, ymax,
                               0.01, 10)
                               
            for xi, yi, zi in zip(a_list1, a_list2, R_array_angle):
                h1.Fill(xi, yi, zi)
                
            h2 = ROOT.TProfile2D("Profile 1", "Pz difference",
                                x_bins, xmin, xmax,
                                y_bins, ymin, ymax,
                                0.01, 10)
                                
            for xi, yi, zi in zip(a_list1, a_list2, R_array_diff):
                h2.Fill(xi, yi, zi)
                
            N = 50
            stops_raw = np.logspace(-2, 0, N, base=10)
            stops = (stops_raw - stops_raw.min()) / (stops_raw.max() - stops_raw.min())

            cmap = matplotlib.colormaps.get_cmap('inferno')
            red   = np.array([cmap(x)[0] for x in stops], dtype=np.float64)
            green = np.array([cmap(x)[1] for x in stops], dtype=np.float64)
            blue  = np.array([cmap(x)[2] for x in stops], dtype=np.float64)

            ROOT.TColor.CreateGradientColorTable(N, stops, red, green, blue, 255)
            ROOT.gStyle.SetNumberContours(255)
            
            ROOT.gStyle.SetOptStat(0)
            #ROOT.gStyle.SetPalette(ROOT.kDarkBodyRadiator)

            canvas = ROOT.TCanvas("canvas", "Confronto Profile2D", 1800, 1000)
            canvas.Divide(2, 1)

            canvas.cd(1)
            h1.Draw("COLZ")
            
            h1.GetXaxis().SetTitle(f"{var1}")
            h1.GetYaxis().SetTitle(f"{var2}")
            
            h2.GetXaxis().SetTitle(f"{var1}")
            h2.GetYaxis().SetTitle(f"{var2}")

            canvas.cd(2)
            h2.Draw("COLZ")
            
            filename = os.path.join(output_dir, f"{var1}_{var2}_histo2d_median_.png")

            canvas.SaveAs(filename)
    
    #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    #SINGLE HISTOGRAM
    
    #I method - minimum angle
    #Pz
    fm_pz_array = np.array(fm_pz_list_l)
    fm_pz_h_array = np.array(pz_h)
    
    #Pt
    fm_pt_array = np.array(fm_pt_list_l)
    fm_pt_w_array = np.array(pt_w_l)
    fm_pt_h_array = np.array(pt_h)
    
    #eta
    eta_lep_array = np.array(eta_lep)
    eta_h_array = np.array(eta_H)
    
    #cos_array = np.array(cos_list_l)
    R_array_angle = np.array(R_l_angle)
    R_array_diff  = np.array(R_l_diff)
    
    print([attr for attr in dir(ROOT) if attr.startswith("k")])
    
    h1 = ROOT.TProfile2D("Profile 1", "Minimum Angle",
                           20, 0, 1400,
                           25, 0, 400,
                           -2, 10)
    for xi, yi, zi in zip(fm_pz_array, fm_pt_w_array, R_array_angle):
        h1.Fill(xi, yi, zi)
        
    h2 = ROOT.TProfile2D("Profile 2", "Pz difference",
                           20, 0, 1400,
                           25, 0, 400,
                           -2, 10)
    for xi, yi, zi in zip(fm_pz_array, fm_pt_w_array, R_array_diff):
        h2.Fill(xi, yi, zi)
        
    N = 50
    stops_raw = np.logspace(-2, 0, N, base=10)
    stops = (stops_raw - stops_raw.min()) / (stops_raw.max() - stops_raw.min())

    cmap = matplotlib.colormaps.get_cmap('inferno')
    red   = np.array([cmap(x)[0] for x in stops], dtype=np.float64)
    green = np.array([cmap(x)[1] for x in stops], dtype=np.float64)
    blue  = np.array([cmap(x)[2] for x in stops], dtype=np.float64)

    ROOT.TColor.CreateGradientColorTable(N, stops, red, green, blue, 255)
    ROOT.gStyle.SetNumberContours(255)
    ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetPalette(ROOT.kRainbow)

    canvas = ROOT.TCanvas("canvas", "Confronto Profile2D", 1600, 1000)
    canvas.Divide(2, 1)

    canvas.cd(1)
    h1.Draw("COLZ")
    
    h1.GetXaxis().SetTitle("pz lepton (GeV)")
    h1.GetYaxis().SetTitle("pt W (GeV)")
    
    h2.GetXaxis().SetTitle("pz lepton (GeV)")
    h2.GetYaxis().SetTitle("pt W (GeV)")

    canvas.cd(2)
    h2.Draw("COLZ")

    canvas.SaveAs("ptW_pzlepton_profile2d.png")

    #--------------------------------------------------------------------------------------------------------------------

    #pt lepton, pz lepton - root    
    
    h1 = ROOT.TProfile2D("Profile 1", "Minimum Angle",
                           25, 0, 1400,
                           25, 0, 400,
                           -2, 10)
    for xi, yi, zi in zip(fm_pz_array, fm_pt_array, R_array_angle):
        h1.Fill(xi, yi, zi)
        
    h2 = ROOT.TProfile2D("Profile 1", "Pz difference",
                           25, 0, 1400,
                           25, 0, 400,
                           -2, 10)
    for xi, yi, zi in zip(fm_pz_array, fm_pt_array, R_array_diff):
        h2.Fill(xi, yi, zi)
        
    ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetPalette(ROOT.kDarkBodyRadiator)

    canvas = ROOT.TCanvas("canvas", "Confronto Profile2D", 1600, 1000)
    canvas.Divide(2, 1)

    canvas.cd(1)
    h1.Draw("COLZ")
    
    h1.GetXaxis().SetTitle("pz lepton (GeV)")
    h1.GetYaxis().SetTitle("pt lepton (GeV)")
    
    h2.GetXaxis().SetTitle("pz lepton (GeV)")
    h2.GetYaxis().SetTitle("pt lepton (GeV)")

    canvas.cd(2)
    h2.Draw("COLZ")

    canvas.SaveAs("ptlepton_pzlepton_profile2d.png")

    #--------------------------------------------------------------------------------------------------------------------

    #eta lepton & pz lepton

    h1 = ROOT.TProfile2D("Profile 1", "Minimum Angle",
                           25, 0, 1400,
                           25, -6, 6,
                           -2, 10)
    for xi, yi, zi in zip(fm_pz_array, eta_lep_array, R_array_angle):
        h1.Fill(xi, yi, zi)
        
    h2 = ROOT.TProfile2D("Profile 1", "Pz difference",
                           25, 0, 1400,
                           25, -6, 6,
                           -2, 10)
    for xi, yi, zi in zip(fm_pz_array, eta_lep_array, R_array_diff):
        h2.Fill(xi, yi, zi)
        
    ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetPalette(ROOT.kDarkBodyRadiator)

    canvas = ROOT.TCanvas("canvas", "Confronto Profile2D", 1600, 1000)
    canvas.Divide(2, 1)

    canvas.cd(1)
    h1.Draw("COLZ")
    
    h1.GetXaxis().SetTitle("pz lepton (GeV)")
    h1.GetYaxis().SetTitle("eta lepton (rad)")
    
    h2.GetXaxis().SetTitle("pz lepton (GeV)")
    h2.GetYaxis().SetTitle("eta lepton (rad)")

    canvas.cd(2)
    h2.Draw("COLZ TEXT")

    canvas.SaveAs("etalepton_pzlepton_profile2d.png")
    
    
    #-------------------------------------------------------------------------------------------------
    
    #pz lepton pt lepton - median
    
    mask = (fm_pz_array >= 0) & (fm_pz_array <= 1500)

    x_filtered = fm_pz_array[mask]
    y_filtered = fm_pt_array[mask]
    #y_filtered = fm_pt_w_array[mask]
    z_filtered_1 = R_array_angle[mask]
    z_filtered_2 = R_array_diff[mask]

    x_bins = 20
    y_bins = 25
    x_edges = np.linspace(0, 1500, x_bins + 1)
    y_edges = np.linspace(0, 400, y_bins + 1)
    
    x_idx = np.searchsorted(x_edges, x_filtered, side='right') - 1
    y_idx = np.searchsorted(y_edges, y_filtered, side='right') - 1

    valid = (x_idx >= 0) & (x_idx < len(x_edges)-1) & (y_idx >= 0) & (y_idx < len(y_edges)-1)

    x_idx = x_idx[valid]
    y_idx = y_idx[valid]
    z_valid_1 = z_filtered_1[valid]
    z_valid_2 = z_filtered_2[valid]

    bin_values_1 = [[[] for _ in range(len(y_edges)-1)] for _ in range(len(x_edges)-1)]

    for xi, yi, zi in zip(x_idx, y_idx, z_valid_1):
        bin_values_1[xi][yi].append(zi)

    median_z_1 = np.full((len(x_edges)-1, len(y_edges)-1), np.nan)
    for i in range(len(x_edges)-1):
        for j in range(len(y_edges)-1):
            if bin_values_1[i][j]:
                median_z_1[i, j] = np.median(bin_values_1[i][j])
                
    bin_values_2 = [[[] for _ in range(len(y_edges)-1)] for _ in range(len(x_edges)-1)]

    for xi, yi, zi in zip(x_idx, y_idx, z_valid_2):
        bin_values_2[xi][yi].append(zi)

    median_z_2 = np.full((len(x_edges)-1, len(y_edges)-1), np.nan)
    for i in range(len(x_edges)-1):
        for j in range(len(y_edges)-1):
            if bin_values_2[i][j]:
                median_z_2[i, j] = np.median(bin_values_2[i][j])

    
    vmin = 0.1
    vmax = 10

    #title
    plt.figure(figsize=(16, 9))
    
    plt.suptitle('Pz ratio - minimum angle & Pz difference', fontsize=16)

    plt.subplot(1, 2, 1)
    plt.pcolormesh(x_edges, y_edges, median_z_1.T, norm=LogNorm(vmin=vmin, vmax=vmax), cmap='inferno')
    plt.colorbar()
    plt.xlabel('Pz lepton (GeV)', fontsize = 15)
    plt.ylabel('Pt lepton (GeV)', fontsize = 15)
    #plt.ylabel('Pt W (GeV)', fontsize = 15)
    plt.title('Minimum Angle', fontsize = 15)


    plt.subplot(1, 2, 2)
    plt.pcolormesh(x_edges, y_edges, median_z_2.T, norm=LogNorm(vmin=vmin, vmax=vmax), cmap='inferno')
    plt.colorbar()
    plt.xlabel('Pz lepton (GeV)', fontsize = 15)
    plt.ylabel('Pt lepton (GeV)', fontsize = 15)
    #plt.ylabel('Pt W (GeV)', fontsize = 15)
    plt.title('Pz difference', fontsize = 15)

    plt.tight_layout()
    plt.savefig("/home/giorgia/MG/Risultati/Confronto_L-T/Pz_neutrino/pz_lepton_pt_lepton_histogram", dpi=300) #plt.savefig("Mean_median_angle", dpi=300)
    plt.show()
    
  
    #--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
    
    #HISTOGRAMS - matplotlib
    
    #Pz
    fm_pz_array = np.array(fm_pz_list_l)
    fm_pz_h_array = np.array(pz_h)
    
    #Pt
    fm_pt_array = np.array(fm_pt_list_l)
    fm_pt_w_array = np.array(pt_w_l)
    fm_pt_h_array = np.array(pt_h)
    
    #eta
    eta_lep_array = np.array(eta_lep)
    eta_h_array = np.array(eta_H)
    
    #cos_array = np.array(cos_list_l)
    R_array_angle = np.array(R_l_angle)
    R_array_diff  = np.array(R_l_diff)
    
    array_list = [fm_pz_array, fm_pz_h_array, fm_pt_array, fm_pt_w_array, fm_pt_h_array, eta_lep_array, eta_h_array]
    var_name = ["pz lepton (GeV)", "pz higgs (GeV)", "pt lepton (GeV)", "pt W (GeV)", "pt higgs (GeV)", "eta lepton", "eta higgs"]
    
    for (a_list1, var1) in zip(array_list, var_name):
        x_bins = 40
        y_bins = 40
        x_edges = np.linspace(np.min(a_list1), np.max(a_list1), x_bins + 1)
        
        
        for (a_list2, var2) in zip(array_list, var_name):
            if var1 == var2:
                continue
                
            output_dir = "/home/giorgia/MG/Risultati/Confronto_L-T/Pz_neutrino/histo_2D"
            os.makedirs(output_dir, exist_ok=True)
            filename_inverted = os.path.join(output_dir, f"{var2}_{var1}_histo2d_median_.png")

            if os.path.exists(filename_inverted):
                print(f"{filename_inverted} already created")
                continue

            
            y_edges = np.linspace(np.min(a_list2), np.max(a_list2), y_bins + 1)
            
            x = a_list1
            y = a_list2
            z1 = R_array_angle
            z2 = R_array_diff
            
            x_idx = np.searchsorted(x_edges, x, side='right') - 1
            y_idx = np.searchsorted(y_edges, y, side='right') - 1

            valid = (x_idx >= 0) & (x_idx < len(x_edges)-1) & (y_idx >= 0) & (y_idx < len(y_edges)-1)

            x_idx = x_idx[valid]
            y_idx = y_idx[valid]
            z_valid1 = z1[valid]
            z_valid2 = z2[valid]

            bin_values1 = [[[] for _ in range(len(y_edges)-1)] for _ in range(len(x_edges)-1)]

            for xi, yi, zi in zip(x_idx, y_idx, z_valid1):
                bin_values1[xi][yi].append(zi)
                
            bin_values2 = [[[] for _ in range(len(y_edges)-1)] for _ in range(len(x_edges)-1)]

            for xi, yi, zi in zip(x_idx, y_idx, z_valid2):
                bin_values2[xi][yi].append(zi)

            
            median_z1 = np.full((len(x_edges)-1, len(y_edges)-1), np.nan)
            for i in range(len(x_edges)-1):
                for j in range(len(y_edges)-1):
                    if bin_values1[i][j]:
                        median_z1[i, j] = np.median(bin_values1[i][j])
                        
            median_z2 = np.full((len(x_edges)-1, len(y_edges)-1), np.nan)
            for i in range(len(x_edges)-1):
                for j in range(len(y_edges)-1):
                    if bin_values2[i][j]:
                        median_z2[i, j] = np.median(bin_values2[i][j])
                        
            vmin = 0.01
            vmax = 10
            
            plt.figure(figsize=(16, 9))
            plt.suptitle('Pz ratio - minimum angle & pz difference', fontsize=16)
            plt.subplot(1, 2, 1)
            plt.pcolormesh(x_edges, y_edges, median_z1.T, norm=LogNorm(vmin=vmin, vmax=vmax), cmap='inferno')
            plt.colorbar()
            plt.xlabel(f"{var1}", fontsize=14)
            plt.ylabel(f"{var2}", fontsize=14)
            plt.title('Minimum angle', fontsize=15)
        
            plt.subplot(1, 2, 2)
            plt.pcolormesh(x_edges, y_edges, median_z2.T, norm=LogNorm(vmin=max(vmin, 1e-3), vmax=vmax), cmap='inferno')
            plt.colorbar()
            plt.xlabel(f"{var1}", fontsize=14)
            plt.ylabel(f"{var2}", fontsize=14)
            plt.title('Pz difference', fontsize=15)
              

            filename = os.path.join(output_dir, f"{var1}_{var2}_histo2d_median_.png")

            plt.tight_layout()
            plt.savefig(filename, dpi=300)
            plt.show()
    
    
if __name__ == '__main__':
    main()  
                    
