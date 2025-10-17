import pandas as pd
import pickle
from utils_nonlhe import *
import matplotlib.pyplot as plt
import os
import ROOT


def main():
    
 
    #read pkl
    
    with open('/home/giorgia/MG/Risultati/H > bb/dataframes_z_long_selected.pkl', 'rb') as f:
        z_longitudinale = pickle.load(f)
        
    with open('/home/giorgia/MG/Risultati/H > bb/dataframes_w_long_selected.pkl', 'rb') as f:
        w_longitudinale = pickle.load(f)

    
    print("READING PICKLE")
    events = z_longitudinale["events"]    
    names  = z_longitudinale["names"]     
    events_w = w_longitudinale["events"]  
    names_w  = w_longitudinale["names"]   


    #z
    idx_pt_leptone_z = names.index("Pt lepton")
    idx_pt_aleptone_z = names.index("Pt antilepton")
    
    idx_eta_leptone_z = names.index("Eta lepton")
    idx_eta_aleptone_z = names.index("Eta antilepton")
    
    idx_phi_leptone_z = names.index("Phi lepton")
    idx_phi_aleptone_z = names.index("Phi antilepton")
    
    idx_pt_z = names.index("Pt Z")
    #idx_pt_h = names.index("Pt H")
    
    idx_eta_z = names.index("Eta Z")
    #idx_eta_h = names.index("Eta H")
    
    idx_phi_z = names.index("Phi Z")
    #idx_phi_h = names.index("Phi H")
    
    idx_m_zh = names.index("ZH invariant mass")
    
    idx_pt_balance_z = names.index("Pt balance")
    idx_theta_star_z = names.index("cos(θ*)")
    idx_theta_1_z = names.index("cos(θ1)")
    idx_phi_1_z = names.index("φ1")
    idx_deltaphi_z = names.index("Δφ")
    
   
    #w
    idx_pt_leptone_w     = names_w.index("Pt lepton")
    idx_pt_aleptone_w    = names_w.index("Pt antilepton")

    idx_eta_leptone_w    = names_w.index("Eta lepton")
    idx_eta_aleptone_w   = names_w.index("Eta antilepton")

    idx_phi_leptone_w    = names_w.index("Phi lepton")
    idx_phi_aleptone_w   = names_w.index("Phi antilepton")

    idx_pt_w             = names_w.index("Pt W")
    idx_eta_w            = names_w.index("Eta W")
    idx_phi_w            = names_w.index("Phi W")

    idx_m_wh             = names_w.index("WH invariant mass")

    idx_pt_balance_w     = names_w.index("Pt balance")
    idx_theta_star_w     = names_w.index("cos(θ*)")
    idx_theta_1_w        = names_w.index("cos(θ1)")
    idx_phi_1_w          = names_w.index("φ1")
    idx_deltaphi_w       = names_w.index("Δφ")
    
    
    
    # Pt
    pt_leptone_l_z     = events[:, idx_pt_leptone_z]
    pt_aleptone_l_z    = events[:, idx_pt_aleptone_z]
    pt_z_l             = events[:, idx_pt_z]
    #pt_h_l             = events[:, idx_pt_h]
    pt_balance_l_z     = events[:, idx_pt_balance_z]

    # Eta
    eta_leptone_l_z    = events[:, idx_eta_leptone_z]
    eta_aleptone_l_z   = events[:, idx_eta_aleptone_z]
    eta_z_l            = events[:, idx_eta_z]
    #eta_h_l            = events[:, idx_eta_h]

    # Phi
    phi_leptone_l_z    = events[:, idx_phi_leptone_z]
    phi_aleptone_l_z   = events[:, idx_phi_aleptone_z]
    phi_z_l            = events[:, idx_phi_z]
    #phi_h_l            = events[:, idx_phi_h]
    phi_1_l_z          = events[:, idx_phi_1_z]
    deltaphi_l_z       = events[:, idx_deltaphi_z]

    # Angoli e massa
    theta_star_l_z     = events[:, idx_theta_star_z]
    theta_1_l_z        = events[:, idx_theta_1_z]
    m_zh_l             = events[:, idx_m_zh]
        
    # Pt
    pt_leptone_l_w     = events_w[:, idx_pt_leptone_w]
    pt_aleptone_l_w    = events_w[:, idx_pt_aleptone_w]
    pt_w_l             = events_w[:, idx_pt_w]
    pt_balance_l_w     = events_w[:, idx_pt_balance_w]

    # Eta
    eta_leptone_l_w    = events_w[:, idx_eta_leptone_w]
    eta_aleptone_l_w   = events_w[:, idx_eta_aleptone_w]
    eta_w_l            = events_w[:, idx_eta_w]

    # Phi
    phi_leptone_l_w    = events_w[:, idx_phi_leptone_w]
    phi_aleptone_l_w   = events_w[:, idx_phi_aleptone_w]
    phi_w_l            = events_w[:, idx_phi_w]
    phi_1_l_w          = events_w[:, idx_phi_1_w]
    deltaphi_l_w       = events_w[:, idx_deltaphi_w]

    # Angoli e massa
    theta_star_l_w     = events_w[:, idx_theta_star_w]
    theta_1_l_w        = events_w[:, idx_theta_1_w]
    m_wh_l             = events_w[:, idx_m_wh]
        
        
    
    #transverse
    with open('/home/giorgia/MG/Risultati/H > bb/dataframes_z_tr_selected.pkl', 'rb') as f:
        z_trasversale = pickle.load(f)
        
    with open('/home/giorgia/MG/Risultati/H > bb/dataframes_w_tr_selected.pkl', 'rb') as f:
        w_trasversale = pickle.load(f)
    
    events = z_trasversale["events"]     
    names  = z_trasversale["names"]      
    events_w = w_trasversale["events"]   
    names_w  = w_trasversale["names"]  
    
    
    # Pt
    pt_leptone_t_z     = events[:, idx_pt_leptone_z]
    pt_aleptone_t_z    = events[:, idx_pt_aleptone_z]
    pt_z_t             = events[:, idx_pt_z]
    #pt_h_t             = events[:, idx_pt_h]
    pt_balance_t_z     = events[:, idx_pt_balance_z]

    # Eta
    eta_leptone_t_z    = events[:, idx_eta_leptone_z]
    eta_aleptone_t_z   = events[:, idx_eta_aleptone_z]
    eta_z_t            = events[:, idx_eta_z]
    #eta_h_t            = events[:, idx_eta_h]

    # Phi
    phi_leptone_t_z    = events[:, idx_phi_leptone_z]
    phi_aleptone_t_z   = events[:, idx_phi_aleptone_z]
    phi_z_t            = events[:, idx_phi_z]
    #phi_h_t            = events[:, idx_phi_h]
    phi_1_t_z          = events[:, idx_phi_1_z]
    deltaphi_t_z       = events[:, idx_deltaphi_z]

    # Angoli e massa
    theta_star_t_z     = events[:, idx_theta_star_z]
    theta_1_t_z        = events[:, idx_theta_1_z]
    m_zh_t             = events[:, idx_m_zh]
    
    
    # Pt
    pt_leptone_t_w     = events_w[:, idx_pt_leptone_w]
    pt_aleptone_t_w    = events_w[:, idx_pt_aleptone_w]
    pt_w_t             = events_w[:, idx_pt_w]
    pt_balance_t_w     = events_w[:, idx_pt_balance_w]

    # Eta
    eta_leptone_t_w    = events_w[:, idx_eta_leptone_w]
    eta_aleptone_t_w   = events_w[:, idx_eta_aleptone_w]
    eta_w_t            = events_w[:, idx_eta_w]

    # Phi
    phi_leptone_t_w    = events_w[:, idx_phi_leptone_w]
    phi_aleptone_t_w   = events_w[:, idx_phi_aleptone_w]
    phi_w_t            = events_w[:, idx_phi_w]
    phi_1_t_w          = events_w[:, idx_phi_1_w]
    deltaphi_t_w       = events_w[:, idx_deltaphi_w]

    # Angoli e massa
    theta_star_t_w     = events_w[:, idx_theta_star_w]
    theta_1_t_w        = events_w[:, idx_theta_1_w]
    m_wh_t             = events_w[:, idx_m_wh]
        
    
    
    print("BUILDING DICTIONARIES")
    variables_z = {
        "Pt lepton": (pt_leptone_l_z, pt_leptone_t_z, "#it{p}_{T} (GeV)"),
        "Pt antilepton": (pt_aleptone_l_z, pt_aleptone_t_z, "#it{p}_{T} (GeV)"),
        "Eta lepton": (eta_leptone_l_z, eta_leptone_t_z, "#eta"),
        "Eta antilepton": (eta_aleptone_l_z, eta_aleptone_t_z, "#eta"),
        "Phi lepton": (phi_leptone_l_z, phi_leptone_t_z, "#phi (rad)"),
        "Phi antilepton": (phi_aleptone_l_z, phi_aleptone_t_z, "#phi (rad)"),
        "Pt Z": (pt_z_l, pt_z_t, "#it{p}_{T}^{Z} (GeV)"),
        # "Pt H": (pt_h_l, pt_h_t, "#it{p}_{T}^{H} (GeV)"),
        "Eta Z": (eta_z_l, eta_z_t, "#eta^{Z}"),
        # "Eta H": (eta_h_l, eta_h_t, "#eta^{H}"),
        "Phi Z": (phi_z_l, phi_z_t, "#phi^{Z} (rad)"),
        # "Phi H": (phi_h_l, phi_h_t, "#phi^{H} (rad)"),
        "Pt balance": (pt_balance_l_z, pt_balance_t_z, "lepton #it{p}_{T} balance"),
        "cos(θ*)": (theta_star_l_z, theta_star_t_z, "cos(#theta^{*})"),
        "cos(θ1)": (theta_1_l_z, theta_1_t_z, "cos(#theta_{1})"),
        "φ1": (phi_1_l_z, phi_1_t_z, "#phi_{1} (rad)"),
        "Δφ": (deltaphi_l_z, deltaphi_t_z, "|#Delta#phi| (rad)"),
        "ZH invariant mass": (m_zh_l, m_zh_t, "M_{ZH} (GeV)")
    }

    variables_w = {
        "Pt lepton": (pt_leptone_l_w, pt_leptone_t_w, "#it{p}_{T} (GeV)"),
        "Pt antilepton": (pt_aleptone_l_w, pt_aleptone_t_w, "#it{p}_{T} (GeV)"),
        "Eta lepton": (eta_leptone_l_w, eta_leptone_t_w, "#eta"),
        "Eta antilepton": (eta_aleptone_l_w, eta_aleptone_t_w, "#eta"),
        "Phi lepton": (phi_leptone_l_w, phi_leptone_t_w, "#phi (rad)"),
        "Phi antilepton": (phi_aleptone_l_w, phi_aleptone_t_w, "#phi (rad)"),
        "Pt W": (pt_w_l, pt_w_t, "#it{p}_{T}^{W} (GeV)"),
        "Eta W": (eta_w_l, eta_w_t, "#eta^{W}"),
        "Phi W": (phi_w_l, phi_w_t, "#phi^{W} (rad)"),
        "Pt balance": (pt_balance_l_w, pt_balance_t_w, "lepton #it{p}_{T} balance"),
        "cos(θ*)": (theta_star_l_w, theta_star_t_w, "cos(#theta^{*})"),
        "cos(θ1)": (theta_1_l_w, theta_1_t_w, "cos(#theta_{1})"),
        "φ1": (phi_1_l_w, phi_1_t_w, "#phi_{1} (rad)"),
        "Δφ": (deltaphi_l_w, deltaphi_t_w, "|#Delta#phi| (rad)"),
        "WH invariant mass": (m_wh_l, m_wh_t, "M_{WH} (GeV)")
    }

    
    
    
    #higgs
    LHE_L = "unweighted_events_L.lhe"
    events_dict = read_file(LHE_L)
    apply_smearing = True
    selected_id = apply_selections(events_dict)
    
    selected_events = {}
    for e_id in selected_id:
        selected_events[e_id] = events_dict[e_id]
    
    fm_h = get_fm_of([25], selected_events, apply_smearing = apply_smearing)
    
    pt_H_L = get_pt_of([25], selected_events, apply_smearing = apply_smearing)
    eta_H_L = get_eta_of([25], selected_events, apply_smearing = apply_smearing)
    phi_H_L = get_phi_of([25], selected_events, apply_smearing = apply_smearing)
    pt_h_l_list = [val for (i, val) in pt_H_L.items()]
    eta_h_l_list = [val for (i, val) in eta_H_L.items()]
    phi_h_l_list = [val for (i, val) in phi_H_L.items()]
    
    LHE_T = "unweighted_events_T.lhe"
    events_dict = read_file(LHE_T)
    apply_smearing = True
    selected_id = apply_selections(events_dict)
    
    print(len(selected_id))
    
    selected_events = {}
    for e_id in selected_id:
        selected_events[e_id] = events_dict[e_id]
    
    fm_h = get_fm_of([25], selected_events, apply_smearing = apply_smearing)
    
    pt_H_T = get_pt_of([25], selected_events, apply_smearing = apply_smearing)
    eta_H_T = get_eta_of([25], selected_events, apply_smearing = apply_smearing)
    phi_H_T = get_phi_of([25], selected_events, apply_smearing = apply_smearing)
    pt_h_t_list = [val for (i, val) in pt_H_T.items()]
    eta_h_t_list = [val for (i, val) in eta_H_T.items()]
    phi_h_t_list = [val for (i, val) in phi_H_T.items()]
    
    ROOT_hist1d(pt_h_l_list, pt_h_t_list, "Higgs transverse momentum", "pt (GeV)", ylim = 0.17)
    ROOT_hist1d(eta_h_l_list, eta_h_t_list, "Higgs pseudorapidity", "eta", nbins = 30, ylim = 0.05)
    ROOT_hist1d(phi_h_l_list, phi_h_t_list, "Higgs phi", "phi (rad)", nbins = 25)
    
    ROOT_hist1d(pt_leptone_l_z, pt_leptone_t_z, "Lepton transverse moomentum", "pt (GeV)", nbins = 40, ylim = 0.13)
    ROOT_hist1d(deltaphi_l_z, deltaphi_t_z, "Delta phi", "#Delta #phi (rad)", nbins = 30)
    ROOT_hist1d(theta_star_l_z, theta_star_t_z, "#theta *", "cos(#theta *)", nbins = 25)
    ROOT_hist1d(theta_1_l_z, theta_1_t_z, "#theta 1", "cos(#theta 1)", nbins = 25)

    
    return 0
    
    for (var, (l_list, t_list, axis_name)) in variables_z.items():
    
        output_dir = "/home/giorgia/MG/Risultati/H > bb/plots_selection"
        os.makedirs(output_dir, exist_ok=True)
        
        N_bins1 = sturges (len (l_list))
        N_bins2 = sturges (len (t_list))
               
        hist1 = ROOT.TH1F("hist1", f"{var} distribution", N_bins1, np.min(l_list), np.max(l_list))
        hist2 = ROOT.TH1F("hist2", f"{var} distribution", N_bins2, np.min(l_list), np.max(t_list))
        
        for value in l_list:
            hist1.Fill(value)
        for value in t_list:
            hist2.Fill(value)
            
        #normalisation
        if hist1.Integral() > 0:    
            hist1.Scale(1.0/hist1.Integral())
        
        if hist2.Integral() > 0:
            hist2.Scale(1.0/hist2.Integral())
            
        steelblue  = ROOT.TColor.GetColor("#4682B4")  #steelblue
        firebrick  = ROOT.TColor.GetColor("#B22222")  #firebrick

        hist1.SetLineColor(steelblue)
        hist2.SetLineColor(firebrick)
        hist1.SetLineWidth(2)
        hist2.SetLineWidth(2)
        
        hist1.SetLineColor(steelblue)
        hist1.SetFillColor(steelblue)
        hist1.SetFillStyle(3002)
        
        hist2.SetLineColor(firebrick)
        hist2.SetFillColor(firebrick)
        hist2.SetFillStyle(3002)
        #hist1.GetXaxis().SetRangeUser(0, 600)
        #hist1.GetXaxis().SetRangeUser(0, 600)
        hist1.GetXaxis().SetTitle(f"{axis_name}")
        hist1.GetYaxis().SetTitle("Event Fraction")
        
        ROOT.gROOT.SetBatch(True)
        
        canvas = ROOT.TCanvas("canvas", f"Selected {var} distribution", 1600, 1200)
        hist1.Draw("HIST")
        hist2.Draw("HIST SAME")
        
        legend = ROOT.TLegend(0.6, 0.8, 0.89, 0.89)
        legend.SetFillStyle(0)                
        legend.SetTextFont(42)                
        legend.SetTextSize(0.03)              
        legend.AddEntry(hist1, "longitudinal polarisation", "f")
        legend.AddEntry(hist2, "transverse polarisation", "f")
        legend.Draw()

        canvas.SetFillColor(0)                
        canvas.SetLeftMargin(0.12)
        canvas.SetBottomMargin(0.08)
        
        ROOT.gStyle.SetOptStat(0)

        filename = os.path.join(output_dir, f"{var}_selected.png")
        
        canvas.SaveAs(filename)
        
    
    
    
    
if __name__ == '__main__':
    main()  
                    
