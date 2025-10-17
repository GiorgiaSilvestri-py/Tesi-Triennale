import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
from scipy.stats import wasserstein_distance, ks_2samp
from utils_nonlhe import *
import sys
import ROOT

            
def main():
    
    LHE_L = "unweighted_events_L.lhe"
    events_dict = read_file(LHE_L)
    apply_smearing = True
    selected_id = apply_selections(events_dict)
    
    print(len(selected_id))
    
    selected_events = {}
    for e_id in selected_id:
        selected_events[e_id] = events_dict[e_id]


    ##pt lepton smearing + selections

    pt_lepton = get_pt_of([11, 13, 15], selected_events, apply_smearing = apply_smearing) 
    l_list = [val for (e, val) in pt_lepton.items()]

    LHE_T = "unweighted_events_T.lhe"
    events_dict = read_file(LHE_T)
    apply_smearing = True
    selected_id = apply_selections(events_dict)
    
    print(len(selected_id))
    
    selected_events = {}
    for e_id in selected_id:
        selected_events[e_id] = events_dict[e_id]

    pt_lepton = get_pt_of([11, 13, 15], selected_events, apply_smearing = apply_smearing) 
    t_list = [val for (e, val) in pt_lepton.items()]

    ROOT_hist1d(l_list, t_list, "pt lepton", "pt (GeV)", ylim = 0.11)



    #pt Higgs smearing + selection

    LHE_L = "unweighted_events_L.lhe"
    events_dict = read_file(LHE_L)
    apply_smearing = True
    selected_id = apply_selections(events_dict)
    
    print(len(selected_id))
    
    selected_events = {}
    for e_id in selected_id:
        selected_events[e_id] = events_dict[e_id]

    
    fm_h = get_fm_of([25], selected_events, apply_smearing = apply_smearing)
    print(fm_h)
    
    pt_H_L = get_pt_of([25], selected_events, apply_smearing = apply_smearing)
    pt_h_l_list = [val for (i, val) in pt_H_L.items()]
    
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
    pt_h_t_list = [val for (i, val) in pt_H_T.items()]
    
    xmin = min(min(pt_h_l_list), min(pt_h_t_list))
    xmax = max(max(pt_h_l_list), max(pt_h_t_list)) 
    N_bins_1 = 4 * sturges (len (pt_h_l_list))
    N_bins_2 = 4 * sturges (len (pt_h_t_list))
    #bin_edges = np.linspace (xmin, xmax, N_bins)

    hist1 = ROOT.TH1F("hist1", "Higgs transverse momentum distribution", 50, 0, 600)
    hist2 = ROOT.TH1F("hist2", "Higgs transverse momentum distribution", 50, 0, 600)
    
    for value in pt_h_l_list:
        hist1.Fill(value)
    for value in pt_h_t_list:
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
    hist1.GetXaxis().SetRangeUser(0, 600)
    hist1.GetXaxis().SetRangeUser(0, 600)
    ymax = max(hist1.GetMaximum(), hist2.GetMaximum()) * 1.1  # margine del 10%
    hist1.SetMaximum(ymax)
    hist2.SetMaximum(ymax)
    #hist1.SetMaximum(0.125)  # Y-axis upper limit
    #hist2.SetMaximum(0.125)  # Y-axis upper limit

    hist1.GetXaxis().SetTitle("pt (GeV)")
    hist1.GetYaxis().SetTitle("Event Fraction")
    
    ROOT.gROOT.SetBatch(True)
    
    canvas = ROOT.TCanvas("canvas", "Selected Higgs transverse momentum distribution", 1600, 1200)
    hist1.Draw("HIST")
    hist2.Draw("HIST SAME")
    ROOT.gPad.RedrawAxis()  # Ensure axes are updated

    
    legend = ROOT.TLegend(0.6, 0.8, 0.89, 0.89)
    legend.SetFillStyle(0)                # trasparente
    legend.SetTextFont(42)                # font leggibile
    legend.SetTextSize(0.03)              # dimensione testo
    legend.AddEntry(hist1, "longitudinal polarisation", "f")
    legend.AddEntry(hist2, "transverse polarisation", "f")
    legend.Draw()

    canvas.SetFillColor(0)                # sfondo bianco
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.08)
    
    ROOT.gStyle.SetOptStat(0)

    
    canvas.SaveAs("pt_Higgs_selected.png")
    
    
    
    
    
    
    return 0
    
    
    
    
    
    
    
    
    
    
    apply_smearing = False
    print("---- BUILDING DICTIONARIES ----")
    
    
    LHE_L = "unweighted_events_L.lhe"
    events_dict = read_file(LHE_L)
    
    fm_h_dictionary = get_fm_of([25], events_dict, apply_smearing = False)
    
    
    
    #invariant mass check
    inv_mass_list = []
    
    for e_id in fm_h_dictionary.keys():
        fm = fm_h_dictionary[e_id]
        inv_mass_list.append(fm.M)
        
    #plot
    hist = ROOT.TH1F("hist", "Higgs decay products invariant mass distribution", 40, np.min(inv_mass_list), np.max(inv_mass_list))
    
    for val in inv_mass_list:
        hist.Fill(val)
        
    #normalisation
    if hist.Integral() > 0:
        hist.Scale(1.0/hist.Integral())
        
    steelblue  = ROOT.TColor.GetColor("#4682B4")  # steelblue
    firebrick  = ROOT.TColor.GetColor("#B22222")  # firebrick

    hist.SetLineColor(steelblue)
    hist.SetLineWidth(2)
    
    hist.SetLineColor(steelblue)
    hist.SetFillColor(steelblue)
    hist.SetFillStyle(3002)
    
    hist.GetXaxis().SetTitle("m (GeV)")
    hist.GetYaxis().SetTitle("Event Fraction")
    
    ROOT.gROOT.SetBatch(True)
    
    canvas = ROOT.TCanvas("canvas", "Higgs decay products invariant mass distribution", 1600, 1200)
    hist.Draw("HIST")
    
    '''
    legend = ROOT.TLegend(0.6, 0.8, 0.89, 0.89)
    legend.SetFillStyle(0)                # trasparente
    legend.SetTextFont(42)                # font leggibile
    legend.SetTextSize(0.03)              # dimensione testo
    legend.AddEntry(hist1, "longitudinal polarisation", "f")
    legend.AddEntry(hist2, "transverse polarisation", "f")
    legend.Draw()
    '''
    canvas.SetFillColor(0)                # sfondo bianco
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.08)

    ROOT.gStyle.SetOptStat(0)  
    
    canvas.SaveAs("H_invariant_mass.png")
   
   
    #transverse momentum
    pt_H_L = get_pt_of([25], events_dict)
    pt_h_l_list = [val for (i, val) in pt_H_L.items()]
    
    LHE_T = "unweighted_events_T.lhe"
    events_dict = read_file(LHE_T)
    
    pt_H_T = get_pt_of([25], events_dict)
    pt_h_t_list = [val for (i, val) in pt_H_T.items()]
    
    hist1 = ROOT.TH1F("hist1", "Higgs transverse momentum distribution", 60, 0, np.max(pt_h_l_list))
    hist2 = ROOT.TH1F("hist2", "Higgs transverse momentum distribution", 50, 0, np.max(pt_h_t_list))
    
    for value in pt_h_l_list:
        hist1.Fill(value)
    for value in pt_h_t_list:
        hist2.Fill(value)
        
    #normalisation
    if hist1.Integral() > 0:
        hist1.Scale(1.0/hist1.Integral())
    
    if hist2.Integral() > 0:
        hist2.Scale(1.0/hist2.Integral())
        
    steelblue  = ROOT.TColor.GetColor("#4682B4")  # steelblue
    firebrick  = ROOT.TColor.GetColor("#B22222")  # firebrick

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
    hist1.GetXaxis().SetRangeUser(0, 600)
    hist1.GetXaxis().SetRangeUser(0, 600)
    hist1.GetXaxis().SetTitle("pt (GeV)")
    hist1.GetYaxis().SetTitle("Event Fraction")
    
    ROOT.gROOT.SetBatch(True)
    
    canvas = ROOT.TCanvas("canvas", "Higgs transverse momentum distribution", 1600, 1200)
    hist1.Draw("HIST")
    hist2.Draw("HIST SAME")
    
    legend = ROOT.TLegend(0.6, 0.8, 0.89, 0.89)
    legend.SetFillStyle(0)                # trasparente
    legend.SetTextFont(42)                # font leggibile
    legend.SetTextSize(0.03)              # dimensione testo
    legend.AddEntry(hist1, "longitudinal polarisation", "f")
    legend.AddEntry(hist2, "transverse polarisation", "f")
    legend.Draw()

    canvas.SetFillColor(0)                # sfondo bianco
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.08)
    
    ROOT.gStyle.SetOptStat(0)

    
    canvas.SaveAs("pt_Higgs_lhe.png")
    
    
    
    
    
    
    
    #smearing
    
    LHE_L = "unweighted_events_L.lhe"
    events_dict = read_file(LHE_L)
    apply_smearing = True
    pt_H_L = get_pt_of([25], events_dict, apply_smearing = apply_smearing)
    pt_h_l_list = [val for (i, val) in pt_H_L.items()]
    
    LHE_T = "unweighted_events_T.lhe"
    events_dict = read_file(LHE_T)
    apply_smearing = True
    pt_H_T = get_pt_of([25], events_dict, apply_smearing = apply_smearing)
    pt_h_t_list = [val for (i, val) in pt_H_T.items()]
    
    hist1 = ROOT.TH1F("hist1", "Higgs transverse momentum distribution", 60, 0, np.max(pt_h_l_list))
    hist2 = ROOT.TH1F("hist2", "Higgs transverse momentum distribution", 50, 0, np.max(pt_h_t_list))
    
    for value in pt_h_l_list:
        hist1.Fill(value)
    for value in pt_h_t_list:
        hist2.Fill(value)
        
    #normalisation
    if hist1.Integral() > 0:
        hist1.Scale(1.0/hist1.Integral())
    
    if hist2.Integral() > 0:
        hist2.Scale(1.0/hist2.Integral())
        
    steelblue  = ROOT.TColor.GetColor("#4682B4")  # steelblue
    firebrick  = ROOT.TColor.GetColor("#B22222")  # firebrick

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
    hist1.GetXaxis().SetRangeUser(0, 600)
    hist1.GetXaxis().SetRangeUser(0, 600)
    hist1.GetXaxis().SetTitle("pt (GeV)")
    hist1.GetYaxis().SetTitle("Event Fraction")
    
    ROOT.gROOT.SetBatch(True)
    
    canvas = ROOT.TCanvas("canvas", "Higgs transverse momentum distribution", 1600, 1200)
    hist1.Draw("HIST")
    hist2.Draw("HIST SAME")
    
    legend = ROOT.TLegend(0.6, 0.8, 0.89, 0.89)
    legend.SetFillStyle(0)                # trasparente
    legend.SetTextFont(42)                # font leggibile
    legend.SetTextSize(0.03)              # dimensione testo
    legend.AddEntry(hist1, "longitudinal polarisation", "f")
    legend.AddEntry(hist2, "transverse polarisation", "f")
    legend.Draw()

    canvas.SetFillColor(0)                # sfondo bianco
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.08)
    
    ROOT.gStyle.SetOptStat(0)

    
    canvas.SaveAs("pt_Higgs_smeared.png")
    
    
    
    return 0
    
    
    
    
    
    
    pt_H_L = get_pt_of([25], events_H_L, apply_smearing = apply_smearing)
    eta_H_L = get_eta_of([25], events_H_L, apply_smearing = apply_smearing)
   
   
    #analisi polarizzazione trasversale
    LHE_T = "unweighted_events_T.lhe"
    events_H_T = read_file(LHE_T)
    #print("N events_T: ", len(events_T))

    pt_H_T = get_pt_of([25], events_H_T, apply_smearing = apply_smearing)
    eta_H_T = get_eta_of([25], events_H_T, apply_smearing = apply_smearing)
    
    pt_h_l_list = [value for (id, value) in pt_H_L.items()]
    pt_h_t_list = [value for (id, value) in pt_H_T.items()]
    
    print("---- BUILDING HISTOGRAMS ----") 
    
    hist1 = ROOT.TH1F("hist1", "Higgs transverse momentum distribution", 40, 0, np.max(pt_h_l_list))
    hist2 = ROOT.TH1F("hist2", "Higgs transverse momentum distribution", 40, 0, np.max(pt_h_t_list))
    
    for value in pt_h_l_list:
        hist1.Fill(value)
    for value in pt_h_t_list:
        hist2.Fill(value)
        
    #normalisation
    if hist1.Integral() > 0:
        hist1.Scale(1.0/hist1.Integral())
    
    if hist2.Integral() > 0:
        hist2.Scale(1.0/hist2.Integral())
        
    steelblue  = ROOT.TColor.GetColor("#4682B4")  # steelblue
    firebrick  = ROOT.TColor.GetColor("#B22222")  # firebrick

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


    
    hist1.GetXaxis().SetTitle("pt (GeV)")
    hist1.GetYaxis().SetTitle("Event Fraction")
    
    ROOT.gROOT.SetBatch(True)
    
    canvas = ROOT.TCanvas("canvas", "Higgs transverse momentum distribution", 1600, 1200)
    hist1.Draw("HIST")
    hist2.Draw("HIST SAME")
    
    legend = ROOT.TLegend(0.6, 0.8, 0.89, 0.89)
    legend.SetFillStyle(0)                # trasparente
    legend.SetTextFont(42)                # font leggibile
    legend.SetTextSize(0.03)              # dimensione testo
    legend.AddEntry(hist1, "longitudinal polarisation", "f")
    legend.AddEntry(hist2, "transverse polarisation", "f")
    legend.Draw()

    canvas.SetFillColor(0)                # sfondo bianco
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.08)
    #canvas.SetTopMargin(0.08)

    ROOT.gStyle.SetOptStat(0)



    
    
    print(ROOT.gROOT.IsBatch())
    
    
    canvas.SaveAs("pt_Higgs2.png")
    
    

        
    
    
    
    #rappresentazione confronto
    data = list(pt_H_L.values())
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)
    
    data = list(eta_H_L.values())
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_bins = int((max(data) - min(data)) / bin_width)
     
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(pt_H_L, bins=nBins, color='royalblue', edgecolor = 'steelblue', stat='probability', ax=ax[0],  label='longitudinal polarisation',  alpha = 0.8)
    sb.histplot(pt_H_T, bins=nBins, color='firebrick', stat='probability', edgecolor = 'firebrick', element='step', linewidth=1.5, alpha = 0.4, ax=ax[0], label='transverse polarisation')
    ax[0].set_xlabel("Pt (GeV)")
    ax[0].set_ylabel("Event fraction")
    ax[0].set_title("Transverse momentum distribution (H)")
    ax[0].legend()
    
    sb.histplot(eta_H_L, bins=n_bins, color='royalblue', edgecolor = 'steelblue', stat='probability', ax=ax[1], label='longitudinal polarisation',  alpha = 0.8)
    sb.histplot(eta_H_T, bins=n_bins, color='firebrick', stat='probability', edgecolor = 'firebrick', element='step', linewidth=1.5,  alpha = 0.4,
             ax=ax[1], label='transverse polarisation')
    ax[1].set_xlabel("η (rad)")
    ax[1].set_ylabel("Event fraction")
    ax[1].set_title("Pseudorapidity distribution (H)")
    ax[1].legend()

    plt.savefig("def_plot/Confronto Higgs.png", dpi=300)
    plt.show()
    
    print("---- COMPUTING METRICS ----")
    #wasserstein distance pt
    distance_z = wasserstein_distance(list(pt_H_L.values()), list(pt_H_T.values()))
    k_s_z = ks_2samp(list(pt_H_L.values()), list(pt_H_T.values()))

    #bin con massima distanza
    height_l, bin_l = np.histogram(list(pt_H_L.values()), bins=40, density=True)
    height_t, bin_t = np.histogram(list(pt_H_T.values()), bins=40, density=True)
    bin_center = [0.5 * (bin_l[i] + bin_l[i+1]) for i in range(len(bin_l) - 1)]
    difference = [abs(x-y) for x,y in zip(height_l, height_t)]
    index_max = np.argmax(difference)
    pt_max = bin_center[index_max]


    print("---- COMPUTING METRICS ----")
    print(f"EMD (pt): {distance_z:.2f} Gev")
    print(f"Kolmogorov-Smirnov stat: {k_s_z.statistic:.3f}, p-value: {k_s_z.pvalue:.3f}")
    print(f"Valore di pt per cui si ha massima differenza: {pt_max:0f}")
    
    #wasserstein distanze eta
    distance_w = wasserstein_distance(list(eta_H_L.values()), list(eta_H_T.values()))
    k_s_w = ks_2samp(list(eta_H_L.values()), list(eta_H_T.values()))
    
    #bin con massima distanza w
    height_l, bin_l = np.histogram(list(eta_H_L.values()), bins=40, density=True)
    height_t, bin_t = np.histogram(list(eta_H_T.values()), bins=40, density=True)
    bin_center = [0.5 * (bin_l[i] + bin_l[i+1]) for i in range(len(bin_l) - 1)]
    difference = [abs(x-y) for x,y in zip(height_l, height_t)]
    index_max = np.argmax(difference)
    eta_max = bin_center[index_max]
    
    print(f"EMD (eta): {distance_w:.2f} rad")
    print(f"Kolmogorov-Smirnov stat: {k_s_w.statistic:.3f}, p-value: {k_s_w.pvalue:.3f}")
    print(f"Valore di eta per cui si ha massima differenza: {eta_max:0f}")
    
    
    
    #confronto massa invariante 
    
    #longitudinale
    LHE_L = "unweighted_events_L.lhe"
    eventsL_dict = read_file(LHE_L)
    
    W_fm_L = get_fm_of([24, -24], eventsL_dict, apply_smearing = apply_smearing)                                 
    Z_fm_L = get_fm_of([23], eventsL_dict, apply_smearing = apply_smearing)                                          
    H_fm_L = get_fm_of([25], eventsL_dict, apply_smearing = apply_smearing)  
    
    #somma dei due e massa invariante
    fm_ZH_L = compute_tot_fm(Z_fm_L, H_fm_L)
    fm_WH_L = compute_tot_fm(W_fm_L, H_fm_L)
    
    M_ZH_l = [zh.M for zh in fm_ZH_L.values() if zh.M > 0]                                        #lista con massa invariante di ZH
    M_WH_l = [wh.M for wh in fm_WH_L.values() if wh.M > 0]                                        #lista con massa invariante di WH
    
   
    #trasversale
    LHE_T = "unweighted_events_T.lhe"
    eventsT_dict = read_file(LHE_T)
    
    W_fm_T = get_fm_of([24, -24], eventsT_dict, apply_smearing = apply_smearing)                                 
    Z_fm_T = get_fm_of([23], eventsT_dict, apply_smearing = apply_smearing)                                          
    H_fm_T = get_fm_of([25], eventsT_dict, apply_smearing = apply_smearing)  
    
    #somma dei due e massa invariante
    fm_ZH_T = compute_tot_fm(Z_fm_T, H_fm_T)
    fm_WH_T = compute_tot_fm(W_fm_T, H_fm_T)
    
    M_ZH_t = [zh.M for zh in fm_ZH_T.values()]                                        #lista con massa invariante di ZH
    M_WH_t = [wh.M for wh in fm_WH_T.values()]                                        #lista con massa invariante di WH
    
    
    
    #rappresentazione sul grafico
    data = M_ZH_l
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    nBins = int((max(data) - min(data)) / bin_width)
    
    data = M_WH_l
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    bin_width = 2 * iqr / (len(data) ** (1/3))
    n_Bins = int((max(data) - min(data)) / bin_width)
    
    print("---- BUILDING HISTOGRAMS ----")
    sb.set(style="whitegrid")
    fig, ax = plt.subplots(1, 2, figsize=(16, 5))
    plt.subplots_adjust(wspace=0.5)
       
    sb.histplot(M_ZH_l, bins=nBins, color='royalblue', edgecolor = 'steelblue', stat='probability', ax=ax[0],  label='longitudinal polarisation',  alpha = 0.8)
    sb.histplot(M_ZH_t, bins=nBins, color='firebrick', stat='probability', edgecolor = 'firebrick', element='step', alpha = 0.4, ax=ax[0], label='transverse polarisation')
    ax[0].set_xlabel("M (GeV)")
    ax[0].set_ylabel("Event fraction")
    ax[0].set_title("ZH invariant mass distribution")
    ax[0].legend()
    
    print(min(M_WH_l), max(M_WH_t))
    
    sb.histplot(M_WH_l, bins=n_Bins, color='royalblue', edgecolor = 'steelblue', stat='probability', ax=ax[0],  label='longitudinal polarisation',  alpha = 0.8)
    sb.histplot(M_WH_t, bins=n_Bins, color='firebrick', stat='probability', edgecolor = 'firebrick', element='step', alpha = 0.4, ax=ax[0], label='transverse polarisation')
    ax[1].set_xlabel("M (GeV)")
    ax[1].set_ylabel("Event fraction")
    ax[1].set_title("WH invariant mass distribution")
    ax[1].legend()

    plt.savefig("def_plot/Confronto massa invariante.png", dpi=300)
    plt.show()
 
    #wasserstein distance zh
    distance_z = wasserstein_distance(M_ZH_l, M_ZH_t)
    k_s_z = ks_2samp(M_ZH_l, M_ZH_t)

    #bin con massima distanza
    height_l, bin_l = np.histogram(M_ZH_l, bins=40, density=True)
    height_t, bin_t = np.histogram(M_ZH_t, bins=40, density=True)
    bin_center = [0.5 * (bin_l[i] + bin_l[i+1]) for i in range(len(bin_l) - 1)]
    difference = [abs(x-y) for x,y in zip(height_l, height_t)]
    index_max = np.argmax(difference)
    m_max = bin_center[index_max]

    print(f"EMD (M_ZH): {distance_z:.2f} Gev")
    print(f"Kolmogorov-Smirnov stat: {k_s_z.statistic:.3f}, p-value: {k_s_z.pvalue:.3f}")
    print(f"Valore di M per cui si ha massima differenza: {m_max:0f}")
    
    #wasserstein distanze wh
    distance_w = wasserstein_distance(M_WH_l, M_WH_t)
    k_s_w = ks_2samp(M_WH_l, M_WH_t)
    
    #bin con massima distanza w
    height_l, bin_l = np.histogram(M_WH_l, bins=40, density=True)
    height_t, bin_t = np.histogram(M_WH_t, bins=40, density=True)
    bin_center = [0.5 * (bin_l[i] + bin_l[i+1]) for i in range(len(bin_l) - 1)]
    difference = [abs(x-y) for x,y in zip(height_l, height_t)]
    index_max = np.argmax(difference)
    m_max = bin_center[index_max]
    
    print(f"EMD (M_WH): {distance_w:.2f} rad")
    print(f"Kolmogorov-Smirnov stat: {k_s_w.statistic:.3f}, p-value: {k_s_w.pvalue:.3f}")
    print(f"Valore di M per cui si ha massima differenza: {m_max:0f}")

if __name__ == '__main__':
    main()  
                    
