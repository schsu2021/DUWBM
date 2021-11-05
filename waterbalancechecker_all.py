# -*- coding: utf-8 -*-
"""
Created on Mon May  4 23:35:51 2020

@author: ShuChen Hsu
"""

import pandas as pd
import numpy as np

def WaterBalanceCheckall(result_all, param_all, output, forcing):
    """
    Args:
        result_all (dataframe): list of results of each grid
        #outdf (dataframe): result at model area outlet 
        param_all (dataframe): list of parameters of each grid
        #IWU(float): daily indoor water consumption, same for each grids for now [L]
    Return:
        df (dataframe): In, Out, dS, WB, In_p, In_ir, totalE
    Compute the water balance of the cluster.
    Input: P * (AR + RTop * ART + AP + APer + SWSop * ASWS + WWSop * AWWS)
          IR [mm AR/AP/APer], LD [mm AGW]
          inflow from upstream (Rw_out, Rs_out, BF_out) [L]
    Output: Ea, ET, Q_seep, Rw, Rs, BF_out
    dS: RSTt, PSTt, PerSTt, ht [mm]
        RTt, SWSt, WWSt [L]
        GWt, GWabt [m-SL]
    
    """
    number = len(result_all)   # number of grids
    dataL = np.shape(result_all[0])[0]   # length of time
    lst = [
        {
            "1_In": np.nan, "1_Out":np.nan, "1_dS":np.nan, "1_WB1":np.nan, "1_WB2":np.nan,
            "2_In_p": np.nan,"2_In_sewer": np.nan,"2_In_imp":np.nan,"2_In_sup": np.nan, "2_IR": np.nan, 
            '2_Out_subsurface':np.nan,"2_totalE": np.nan,"2_Out_ET": np.nan,"2_Out_subsurface": np.nan,
            "2_Out_sewer":np.nan,"2_Out_sup": np.nan
            
        }
    ]
    
    for t in range(1,dataL):
        In = In_p = In_sewer = In_imp = In_sup = IR = 0
        Out = Out_e_rpper = Out_e_st = Out_ET = Out_sewer = Out_subsurface = Out_sup = totalE = 0
        dS = dS_imp = dS_per = dS_rt = dS_uz = dS_gw = dS_s = 0
        
        for i in range(0, number):
            result = result_all[i]
            param = param_all[i]
            #Natural
            In_p = In_p + (forcing.P[t]*(param["AR"]+param["AP"]+param["APer"]+
                                        param["RTop"]*param["allART"]+param["SWSop"]*param["ASWS"]))
            #In_in = In_in + result.LD[t]*param["AGW"]
            In_sewer = In_sewer + (result.Rs_up[t] + result.Rw_up[t]) # + result.BF_up[t])
            In_imp = In_imp + result.Import[t]
            In_sup = In_sup + result.cWWSsup[t] + result.SWSsup[t]
            IR = IR + (result.IR_per[t]*param["APer"] + result.IR_p[t]*param["AP"] 
                             + result.IR_r[t]*param["AR"])
            Out_e_rpper = Out_e_rpper + (result.Ea_r[t]*param["AR"] + result.Ea_p[t]*param["AP"]
                                         + result.Ea_per[t]*param["APer"])
            Out_e_st = Out_e_st + result.Ea_rt[t] + result.Ea_sws[t]# + param["IWU"]
            Out_ET = Out_ET + result.ET[t] * param["AUZ"] 
            Out_subsurface = Out_subsurface + (result.Q_seep[t] + result.BF_out[t]) # + result.BF_down[t])
            Out_sewer = Out_sewer + result.Rw[t] + result.Rs[t]
            Out_sup = Out_sup + result.cWWSuse[t] + result.SWSuse[t]
            totalE = totalE + ((result.Ea_r[t]*param["AR"] + result.Ea_p[t]*param["AP"] + result.Ea_per[t]*param["APer"]) + 
                               result.Ea_rt[t] + result.ET[t] * param["AUZ"] + result.Ea_sws[t])
            dS_imp = dS_imp + (result.RSTt[t]-result.RSTt[t-1])*param["AR"] + (result.PSTt[t]-result.PSTt[t-1])*param["AP"]
            dS_per = dS_per + (result.PerSTt[t]-result.PerSTt[t-1])*param["APer"]
            dS_rt = dS_rt + result.RTt[t]-result.RTt[t-1]
            dS_uz = dS_uz + (result.ht[t]-result.ht[t-1])*param["AUZ"]
            dS_gw = dS_gw + ((result.sc_gw[t] * (result.GWt[t-1] + result.GWabt[t-1] / result.sc_gw[t] - (
                    result.GWt[t] + result.GWabt[t] / result.sc_gw[t])) * 1000) * param["AGW"])
            dS_s = dS_s + ((result.SWSt[t]-result.SWSt[t-1]) 
                           + (result.WWSt[t]-result.WWSt[t-1]) 
                           + (result.cWWSt[t]-result.cWWSt[t-1]))
            #dS_ow = dS_ow - (result.OWt[t] - result.OWt[t-1])*1000*param["AOW"]
        
        In =  In_p + In_imp  #+ In_ir, LD: included in Import and internal water supply
        Out = Out_e_rpper + Out_e_st + Out_ET + Out_subsurface + output.Rs_out[t] + output.Rw_out[t]
        dS = dS_imp + dS_per + dS_rt + dS_uz + dS_gw + dS_s
        
        WB1 = In - Out - dS
        WB2 = (In + In_sewer + In_sup) - (totalE + Out_subsurface + Out_sewer + Out_sup) - dS
        lst.append({"1_In": In, "1_Out":Out, "1_dS":dS, "1_WB1":WB1, "1_WB2":WB2,
                    "2_In_p": In_p,"2_In_sewer": In_sewer,"2_In_imp": In_imp,"2_In_sup": In_sup, "2_IR": IR, 
                    '2_Out_subsurface':Out_subsurface,"2_totalE": totalE,"2_Out_ET": Out_ET,"2_Out_subsurface": Out_subsurface,
                    "2_Out_sewer":Out_sewer,"2_Out_sup": Out_sup})
    
    df = pd.DataFrame(lst)
    df2 = df.set_index(result_all[0].index)
    return df2
    