import matplotlib.pyplot as plt
import math
import numpy as np
import vector
import seaborn as sb
import sys
from ROOT import TGraph, TF1, TCanvas, gStyle
from array import array

            
def main():


    x_coord = [30, 40, 50, 60, 80, 100, 200, 300, 400, 700, 1000, 2000]
    y_coord = [0.225, 0.194, 0.174, 0.16,  0.138, 0.125, 0.092, 0.08, 0.07, 0.061, 0.058, 0.048]
    
    #hyperbole
    graph = TGraph(len(x_coord), array('d', x_coord), array('d', y_coord))
    graph.SetTitle("hyperbolic fit; p_{T}^{jet} [GeV]; Jet energy resolution")
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(4)

    f_hyper = TF1("f_hyper", "[0]/x + [1]", 0, 2000)
    f_hyper.SetParameters(1.0, 0.1)  
    f_hyper.SetLineColor(2)

    graph.Fit(f_hyper, "R")

    c = TCanvas("c", "Fit iperbolico", 1400, 1000)
    gStyle.SetOptFit(1111)  
    #c.SetLogx()
    graph.Draw("AP")
    c.Update()
    c.SaveAs("hyperbolic_jet_resolution_fit.png")  
    
     
    
if __name__ == '__main__':
    main()  
                    
