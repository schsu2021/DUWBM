# -*- coding: utf-8 -*-
"""
Created on Sat May  2 11:45:18 2020

@author: ShuChen Hsu
"""

import pandas as pd
import numpy as np

def WaterBalanceCheck(result_all, param_all, forcing):
    """
    Args:
        result_all (dataframe): list of results of each grid
        param_all (dataframe): list of parameters of each grid
        #IWU(float): daily indoor water consumption, same for each grids for now [L]
    
    Compute the water balance of the cluster.
    Input: P * (AR + RTop * ART + AP + APer + SWSop * ASWS + WWSop * AWWS)
          IR [mm AR/AP/APer], LD [mm AGW]
          inflow from upstream (Rw_out, Rs_out, BF_out) [L]
    Output: Ea, ET, Q_seep, Rw, Rs, BF_out, BF_down
    dS: RSTt, PSTt, PerSTt, ht [mm]
        RTt, SWSt, WWSt [L]
        GWt, GWabt [m-SL]
    
    """
    number = len(result_all)   # number of grids
    wblist = []
    names = locals()
    for i in range(1, number+1):
        result = result_all[i-1]  #(dataframe) result of grid i
        param = param_all[i-1]  #(dataframe) result of grid i
        dataL = np.shape(result)[0]
        names['lst%s' % i] = [
            {
                "In": np.nan,
                "Out": np.nan,
                "dS": np.nan,
                "WB": np.nan,
                "WBuse": np.nan
            }
        ]
        for t in range(1,dataL):
            In_p = (forcing.P[t]*(param["AR"]+param["AP"]+param["APer"]+param["RTop"]*param["allART"]+
                                 param["SWSop"]*param["ASWS"]))
            In_ir = (result.IR_per[t]*param["APer"] + result.IR_p[t]*param["AP"] 
                     + result.IR_r[t]*param["AR"])
            In_sewer = result.Rs_up[t] + result.Rw_up[t]  # + result.LD[t]*param["AGW"]
            In_sup = result.Import[t] + result.cWWSsup[t] + result.SWSsup[t]
            In =  In_p + In_sewer + In_sup #+ In_ir
            
            Out_e_rpper = (result.Ea_r[t]*param["AR"] + result.Ea_p[t]*param["AP"] + result.Ea_per[t]*param["APer"])
            #Out_e_rt = result.Ea_rt[t]
            Out_ET = result.ET[t] * param["AUZ"]
            Out_storage = result.Ea_rt[t] + result.Ea_sws[t]
            Out_sewer = result.Rw[t] + result.Rs[t]
            Out_subsurface = result.Q_seep[t] + result.BF_out[t] # + result.ET[t] * param["AUZ"] + result.BF_down[t]
            #Out_ow = result.Ea_ow[t]*param["AOW"] + result.Qout[t]
            Out_sup = result.cWWSuse[t] + result.SWSuse[t]
            Out = Out_e_rpper + Out_ET + Out_storage + Out_sewer + Out_subsurface + Out_sup
            
            dS_imp = ((result.RSTt[t]-result.RSTt[t-1])*param["AR"] 
                      + (result.PSTt[t]-result.PSTt[t-1])*param["AP"])
            dS_per = (result.PerSTt[t]-result.PerSTt[t-1])*param["APer"]
            dS_rt = result.RTt[t]-result.RTt[t-1]
            dS_uz = (result.ht[t]-result.ht[t-1])*param["AUZ"]
            dS_gw = ((result.sc_gw[t] * (result.GWt[t-1] + result.GWabt[t-1] / result.sc_gw[t] - (
                result.GWt[t] + result.GWabt[t] / result.sc_gw[t])) * 1000) * param["AGW"])
            dS_s = ((result.SWSt[t]-result.SWSt[t-1]) 
                    + (result.WWSt[t]-result.WWSt[t-1]) 
                    + (result.cWWSt[t]-result.cWWSt[t-1]))
            dS = dS_imp + dS_per + dS_rt + dS_uz + dS_gw + dS_s
            
            WB = In - Out - dS
            WBuse = ((In_sup - result.LD[t]*param["AGW"]) + result.RTuse[t] + result.SSGuse[t] + result.WWSuse[t] - In_ir - param["IWU"])
            names['lst%s' % i].append({"In": In, "Out":Out, "dS":dS, "WB":WB, "WBuse":WBuse})
        
        names['df%s' % i] = pd.DataFrame(names['lst%s' % i])
        names['df%s' % i] = names['df%s' % i].set_index(result_all[0].index)
        wblist.append(names['df%s' % i])
    
    return wblist