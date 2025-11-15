from array import array
import os
import sys
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

    # region Load df

    df_ZLH = ROOT.RDataFrame("tree_ZLH", "complete_ZLH.root")
    df_WLH = ROOT.RDataFrame("tree_WLH", "complete_WLH.root")
    df_ZTH = ROOT.RDataFrame("tree_ZTH", "complete_ZTH.root")
    df_WTH = ROOT.RDataFrame("tree_WTH", "complete_WTH.root")
    
    # endregion
    
    #------------------------------------------------------------------------------------------------------------------------------------

    # region Define Cuts
    
    #DEFINE CUTS - VLH
    
    df_ZLH = df_ZLH\
            .Define("cut0_ZLH", '(event_type == "e") || (event_type == "mu")')\
            .Define("cut1_ZLH", 'pt_leading > leading_th && pt_subleading > subleading_th')\
            .Define("cut2_ZLH", '(event_type == "e" && R1 > 1.0 && R2 > 1.0) || (event_type == "mu" && R1 > 0.5 && R2 > 0.5)')
    
    df_WLH = df_WLH\
            .Define("cut0_WLH", 'event_type == "e/mu"')\
            .Define("cut1_WLH", 'pt_leading > leading_th && pt_subleading > subleading_th')\
            .Define("cut2_WLH", 'event_type == "e/mu" && R1 > 1.0 && R2 > 1.0')\
            .Define("cut3_WLH", 'pt_MET > 45')
        
    
    #DEFINE CUTS - VTH
    
    df_ZTH = df_ZTH\
            .Define("cut0_ZTH", '(event_type == "e") || (event_type == "mu")')\
            .Define("cut1_ZTH", 'pt_leading > leading_th && pt_subleading > subleading_th')\
            .Define("cut2_ZTH", '(event_type == "e" && R1 > 1.0 && R2 > 1.0) || (event_type == "mu" && R1 > 0.5 && R2 > 0.5)')
    
    df_WTH = df_WTH\
            .Define("cut0_WTH", 'event_type == "e/mu"')\
            .Define("cut1_WTH", 'pt_leading > leading_th && pt_subleading > subleading_th')\
            .Define("cut2_WTH", 'event_type == "e/mu" && R1 > 1.0 && R2 > 1.0')\
            .Define("cut3_WTH", 'pt_MET > 45')

    # endregion
    
    #--------------------------------------------------------------------------------------------------------------------------------------------------------

    # region Channel Evts df [tot & cumulative]

    #########################################################################################################################################################
    #CHANNEL DATAFRAME
    #########################################################################################################################################################

    N_tot_ZLH    = df_ZLH.Count().GetValue()
    N_lep_ZLH    = df_ZLH.Filter("cut0_ZLH").Count().GetValue()
    N_cut1_ZLH   = df_ZLH.Filter("cut0_ZLH && cut1_ZLH").Count().GetValue()
    N_cut2_ZLH   = df_ZLH.Filter("cut0_ZLH && cut2_ZLH").Count().GetValue()
    N_cut12_ZLH  = df_ZLH.Filter("cut0_ZLH && cut1_ZLH && cut2_ZLH").Count().GetValue()

    N_tot_WLH    = df_WLH.Count().GetValue()
    N_lep_WLH    = df_WLH.Filter("cut0_WLH").Count().GetValue()
    N_cut1_WLH   = df_WLH.Filter("cut0_WLH && cut1_WLH").Count().GetValue()
    N_cut2_WLH   = df_WLH.Filter("cut0_WLH && cut2_WLH").Count().GetValue()
    N_cut3_WLH   = df_WLH.Filter("cut0_WLH && cut3_WLH").Count().GetValue()
    N_cut123_WLH = df_WLH.Filter("cut0_WLH && cut1_WLH && cut2_WLH && cut3_WLH").Count().GetValue()

    N_tot_ZTH = df_ZTH.Count().GetValue()
    N_lep_ZTH    = df_ZTH.Filter("cut0_ZTH").Count().GetValue()
    N_cut1_ZTH   = df_ZTH.Filter("cut0_ZTH && cut1_ZTH").Count().GetValue()
    N_cut2_ZTH   = df_ZTH.Filter("cut0_ZTH && cut2_ZTH").Count().GetValue()
    N_cut12_ZTH  = df_ZTH.Filter("cut0_ZTH && cut1_ZTH && cut2_ZTH").Count().GetValue()

    N_tot_WTH = df_WTH.Count().GetValue()
    N_lep_WTH    = df_WTH.Filter("cut0_WTH").Count().GetValue()
    N_cut1_WTH   = df_WTH.Filter("cut0_WTH && cut1_WTH").Count().GetValue()
    N_cut2_WTH   = df_WTH.Filter("cut0_WTH && cut2_WTH").Count().GetValue()
    N_cut3_WTH   = df_WTH.Filter("cut0_WTH && cut3_WTH").Count().GetValue()
    N_cut123_WTH = df_WTH.Filter("cut0_WTH && cut1_WTH && cut2_WTH && cut3_WTH").Count().GetValue()


    data = [
            ["Total events",            N_tot_ZLH,   N_tot_ZTH,   N_tot_WLH,    N_tot_WTH],
            ["Events with e/mu",        N_lep_ZLH,   N_lep_ZTH,   N_lep_WLH,    N_lep_WTH],
            ["Events after cut1 only",  N_cut1_ZLH,  N_cut1_ZTH,  N_cut1_WLH,   N_cut1_WTH],
            ["Events after cut2 only",  N_cut2_ZLH,  N_cut2_ZTH,  N_cut2_WLH,   N_cut2_WTH],
            ["Events after cut3 only",  "-",         "-",         N_cut3_WLH,   N_cut3_WTH],
            ["Eventi selezionati",      N_cut12_ZLH, N_cut12_ZTH, N_cut123_WLH, N_cut123_WTH]
            ]


    df_channel = pd.DataFrame(data, columns=["", "ZLH", "ZTH", "WLH", "WTH"])

    print("\n" + "="*58)
    print("{:^60}".format("H > γγ"))
    print("="*58 + "\n")
    print(tabulate(df_channel.values.tolist(), headers=df_channel.columns.tolist(), tablefmt="grid"))


    #--------------------------------------------------------------------------------------------------------------------------------------------------------


    #########################################################################################################################################################
    #CHANNEL CUMULATIVE DATAFRAME
    #########################################################################################################################################################

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

    print("\n" + "="*53)
    print("{:^60}".format("H > γγ [cumulative values]"))
    print("="*53 + "\n")
    print(tabulate(df_cumulative.values.tolist(), headers=df_cumulative.columns.tolist(), tablefmt="grid"))

    # endregion


    #--------------------------------------------------------------------------------------------------------------------------------------------------------


    # region Efficiencies and LHC evts number

   
    print("Computing efficiencies")
    
    tot_eff_list_channel  = []
    part_eff_list_channel = []
    lhc_values_channel    = []

    xsection_L = 282            #fb
    BR         = 0.00270
    xsection_T = 252            #fb
    lum        = 100            #1/fb 


    for ch in ["ZLH", "ZTH", "WLH", "WTH"]:

        if ch == "ZLH" or ch == "WLH":
            xsec = xsection_L
        else:
            xsec = xsection_T
     
        #total efficiencies
        epsilon     = df_channel.iloc[0][ch] / 50000
        epsilon_0   = df_channel.iloc[1][ch] / 50000
        epsilon_1   = df_channel.iloc[2][ch] / df_channel.iloc[1][ch]
        epsilon_2   = df_channel.iloc[3][ch] / df_channel.iloc[1][ch]

        #partial efficiencies
        part_epsilon_0  = df_cumulative.iloc[1][ch] / df_cumulative.iloc[0][ch]
        part_epsilon_1  = df_cumulative.iloc[2][ch] / df_cumulative.iloc[1][ch]
        part_epsilon_2  = df_cumulative.iloc[3][ch] / df_cumulative.iloc[2][ch]

        #LHC events number
        N_tot_LHC  = xsec * lum
        N_tot_LHC  = N_tot_LHC  * epsilon                             #ZLH events number
        N_tot_diph = N_tot_LHC  * BR

        N_em_LHC   = N_tot_diph * part_epsilon_0
        N_cut1_LHC = N_em_LHC   * part_epsilon_1
        N_cut2_LHC = N_cut1_LHC * part_epsilon_2

        if df_channel.iloc[4][ch] == "-":
            epsilon_3      = -1
            part_epsilon_3 = -1
            N_cut3_LHC     = -1
        else:
            epsilon_3      = df_channel.iloc[4][ch]    / df_channel.iloc[2][ch]
            part_epsilon_3 = df_cumulative.iloc[4][ch] / df_cumulative.iloc[3][ch]
            N_cut3_LHC     = N_cut2_LHC * part_epsilon_3

        epsilon_tot = df_channel.iloc[5][ch] / df_channel.iloc[0][ch]

        #LHC efficiencies
        kin_epsilon  = N_cut2_LHC / N_em_LHC
        chan_epsilon = N_cut2_LHC / N_tot_diph
        tot_epsilon  = N_cut2_LHC / N_tot_LHC

        #add channel computed values to list
        lhc_values_channel.append(
                                  [
                                   f"{N_tot_LHC:.0f}", f"{N_tot_LHC:.0f}", f"{N_tot_diph:.0f}", f"{N_em_LHC:.0f}", f"{N_cut1_LHC:.0f}", f"{N_cut2_LHC:.0f}", f"{N_cut3_LHC:.0f}",
                                   f"{kin_epsilon:.2f}", f"{chan_epsilon:.2f}", f"{tot_epsilon:.4f}"
                                  ]
                                 )

        tot_eff_list_channel.append(
                                    [f"{epsilon:.2f}", f"{epsilon_0:.2f}", f"{epsilon_1:.2f}", f"{epsilon_2:.2f}", f"{epsilon_3:.2f}", f"{epsilon_tot:.2f}"]
                                   )

        part_eff_list_channel.append(
                                    [f"{epsilon:.2f}", f"{part_epsilon_0:.2f}", f"{part_epsilon_1:.2f}", f"{part_epsilon_2:.2f}", f"{part_epsilon_3:.2f}", f"{epsilon_tot:.2f}"]
                                    )

        
    
    tot_eff_array_channel  = np.array(tot_eff_list_channel)
    part_eff_array_channel = np.array(part_eff_list_channel)
    lhc_values_array       = np.array(lhc_values_channel)


    # endregion


    #########################################################################################################################################################
    #CHANNEL CUMULATIVE and TOTAL EFFICIENCIES
    #########################################################################################################################################################

    # region Cumulative and Total efficiencies

    columns = [
               "",
               "ZLH cumulative eff", "ZLH total eff", 
               "ZTH cumulative eff", "ZTH total eff",
               "WLH cumulative eff", "WLH total eff",
               "WTH cumulative eff", "WTH total eff"
               ]
    
    rows = []

    for eps in range(6):    #efficiencies

        epsilon_row = []

        #build index column
        if eps == 0:
            epsilon_row.append("epsilon channel")
        elif eps == 5:
            epsilon_row.append("epsilon tot")
        else:
            epsilon_row.append(f"epsilon_{eps - 1}")


        for ch in range(4):  #channels

            partial = part_eff_array_channel[ch][eps]               
            total   = tot_eff_array_channel[ch][eps]
            epsilon_row.append(partial)
            epsilon_row.append(total)

        rows.append(epsilon_row)


    df_part_tot = pd.DataFrame(rows, columns = columns)

    print("\n" + "="*183)
    print("{:^183}".format("Cumulative and Total efficiencies"))
    print("="*183 + "\n")

    print(tabulate(df_part_tot.values.tolist(), headers=df_part_tot.columns.tolist(), tablefmt="grid"))

    # endregion

    # region LHC values

    #############################################################################################
    #CHANNEL LHC VALUES
    #############################################################################################

    columns = ["", "ZLH", "ZTH", "WLH", "WTH"]
    row_index = ["Total V{i}H events", 
                "Total ch events", 
                "Total γγ events", 
                "Events with e/mu", 
                "Events after cut1", 
                "Events after cut1, cut2", 
                "Events after cut1, cut2, cut3",
                "Kinematic cuts efficiency",
                "Channel efficiency",
                "Total efficiency"
                ]
            
    rows = []

    for i in range(len(lhc_values_array[1])):

        single_row = []

        #build index column
        single_row.append(row_index[i])

        for ch in range(4):
            single_row.append(lhc_values_array[ch][i]) 

        rows.append(single_row)


    df_LHC = pd.DataFrame(rows, columns = columns)

    print("\n" + "="*82)
    print("{:^82}".format("H > γγ (LHC)"))
    print("="*82 + "\n")

    print(tabulate(df_LHC.values.tolist(), headers=df_LHC.columns.tolist(), tablefmt="grid"))

    # endregion

    #------------------------------------------------------------------------------------------------------------------------

    # region Ordering efficiencies

    print("Ordering cuts")

    cut_list = ["cut1", "cut2", "cut3"]
    channels = ["ZLH", "ZTH", "WLH", "WTH"]
    
    df_list = [df_ZLH, df_ZTH, df_WLH, df_WTH]

    i = 0                                                                               #variable for loop over df_list
    for ch in channels:

        col = f"{ch} total eff"

        #re-arranging
        values     = df_part_tot.loc[2:4, col].copy()
        num_values = [float(v) for v in values if v != -1]

        ordered_val = sorted(num_values)
        print("Ordered total efficiencies:", ordered_val)
        sorted_indices = sorted(range(len(num_values)), key=lambda n: num_values[n])

        ordered_cut_list = [cut_list[n] for n in sorted_indices]
        print(f"Ordered cuts for {col}:", ordered_cut_list)

        j = 0
        for m in range(2, 5): 

            if df_part_tot.at[m, col] == -1:
                continue

            df_part_tot.at[m, col] = ordered_val[j]
            j += 1


        col  = f"{ch} cumulative eff"
        cut0 = "cut0" + f"_{ch}"
        chan_ordered_cut_list = [ordered_cut_list[n] + f"_{ch}" for n in range(len(ordered_cut_list))]

        print(df_list[i])

        #re-define cut order
        df_0 = df_list[i].Filter(cut0)
        df_1 = df_0.Filter(chan_ordered_cut_list[0])
        df_2 = df_1.Filter(chan_ordered_cut_list[1])

        #re-define partial efficiencies
        part_1 = df_1.Count().GetValue() / df_0.Count().GetValue()
        part_2 = df_2.Count().GetValue() / df_1.Count().GetValue()

        #skip if cut doesn't exist
        if df_part_tot.at[4, col] != -1:
            df_3   = df_2.Filter(chan_ordered_cut_list[2])
            part_3 = df_3.Count().GetValue() / df_2.Count().GetValue()

            #overwrite df
            df_part_tot.at[4, col] = part_3

        df_part_tot.at[2, col] = part_1
        df_part_tot.at[3, col] = part_2

        i += 1



    print("\n" + "="*120)
    print("{:^120}".format("Cumulative and Total efficiencies [ordered]"))
    print("="*120 + "\n")

    columns  = [
                "",
                "ZLH cumulative eff", "ZLH total eff", 
                "ZTH cumulative eff", "ZTH total eff",
                "WLH cumulative eff", "WLH total eff",
                "WTH cumulative eff", "WTH total eff"
                ]

    print(tabulate(df_part_tot.values.tolist(), headers=df_part_tot.columns.tolist(), tablefmt="grid"))


    # endregion


if __name__ == '__main__':
    main()