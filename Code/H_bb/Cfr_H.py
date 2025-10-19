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

    ROOT_hist1d(l_list, t_list, "pt lepton", "pt (GeV)")

    


    
    fm_h = get_fm_of([25], selected_events, apply_smearing = apply_smearing)
    
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

    hist1 = ROOT.TH1F("hist1", "Higgs transverse momentum distribution", N_bins_1, xmin, xmax)
    hist2 = ROOT.TH1F("hist2", "Higgs transverse momentum distribution", N_bins_2, xmin, xmax)
    
    for value in pt_h_l_list:
        hist1.Fill(value)
    for value in pt_h_t_list:
        hist2.Fill(value)
        
    #normalisation
    if hist1.Integral() > 0:
        hist1.Scale(1.0/hist1.Integral())
    
    if hist2.Integral() > 0:
        hist2.Scale(1.0/hist2.Integral())
        
    steelblue  = ROOT.TColor.GetColor("#4682B4")
    firebrick  = ROOT.TColor.GetColor("#B22222")

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
    hist1.SetMaximum(ymax)
    hist2.SetMaximum(ymax)
    #hist1.SetMaximum(0.125)  #Y-axis upper limit
    #hist2.SetMaximum(0.125)  #Y-axis upper limit

    hist1.GetXaxis().SetTitle("pt (GeV)")
    hist1.GetYaxis().SetTitle("Event Fraction")
    
    ROOT.gROOT.SetBatch(True)
    
    canvas = ROOT.TCanvas("canvas", "Selected Higgs transverse momentum distribution", 1600, 1200)
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

    
    canvas.SaveAs("pt_Higgs_selected.png")
    
    
    
    
    
    
    
    
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
        
    steelblue  = ROOT.TColor.GetColor("#4682B4")  
    firebrick  = ROOT.TColor.GetColor("#B22222")  

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
    legend.SetFillStyle(0)                
    legend.SetTextFont(42)                
    legend.SetTextSize(0.03)              
    legend.AddEntry(hist1, "longitudinal polarisation", "f")
    legend.AddEntry(hist2, "transverse polarisation", "f")
    legend.Draw()
    '''
    canvas.SetFillColor(0)                
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
        
    steelblue  = ROOT.TColor.GetColor("#4682B4")  
    firebrick  = ROOT.TColor.GetColor("#B22222")  

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
        
    steelblue  = ROOT.TColor.GetColor("#4682B4")  
    firebrick  = ROOT.TColor.GetColor("#B22222")  

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

    
    canvas.SaveAs("pt_Higgs_smeared.png")
    
    
    
    
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
        
    steelblue  = ROOT.TColor.GetColor("#4682B4")
    firebrick  = ROOT.TColor.GetColor("#B22222")

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
    legend.SetFillStyle(0)                
    legend.SetTextFont(42)                
    legend.SetTextSize(0.03)              
    legend.AddEntry(hist1, "longitudinal polarisation", "f")
    legend.AddEntry(hist2, "transverse polarisation", "f")
    legend.Draw()

    canvas.SetFillColor(0)                
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.08)
    #canvas.SetTopMargin(0.08)

    ROOT.gStyle.SetOptStat(0)
    canvas.SaveAs("pt_Higgs2.png")
    
    

if __name__ == '__main__':
    main()  
                    
