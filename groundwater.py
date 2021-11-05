# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:27:58 2020

@author: ShuChen Hsu
"""
import numpy as np
from DUWBM.selector import soil_selector
from DUWBM.gwlcalculator import gwlcalc

class Groundwater:
    """
    Args: 
        AGW(float): area of groundwater storage [m^2] ### assume = AR + AP + APer
        AUZ(float): area of unsaturated zone [m^2]
        LR(float): leakage rate [%]
        GWpret = GW0(float): initial groundwater level (at t=0) [m-SL]
        GWabpret = 0 (float): groundwater level above surface level at previous time step [m-SL]
        seep_def (int): seepage to deep groundwater defined as either constant downward flux or dynamic computed \
        flux which is determined by head difference and resistance [0=flux; 1=level]
        w (float): drainage resistance from groundwater to open water [d]
        vc (float): vertical flow resistance from shallow groundwater to deep groundwater [d]
        IRC(float): infiltration store recession constant ratio [-]
        h_dgw (float): predefined hydraulic head of deep groundwater [m-SL]
        down_seep (float): predefined constant downward flux from shallow groundwater to deep groundwater [mm/d]
        soiltype (int): soil type
        croptype (int): crop type
        dt(float): time step [day]
        soil_prm(dataframe): soil parameter matrix dependent on soil type and crop type
        OW (boolean): if there is open water in this grid
    """
    def __init__(self, dict_param, soilmatrix, etmatrix):
        self.AGW = dict_param["AGW"]
        self.AR = dict_param["AR"]
        self.AP = dict_param["AP"]
        self.APer = dict_param["APer"]
        self.AUZ = dict_param["AUZ"]
        self.LR = dict_param["LR"]/100   # [%] --> [-]
        self.GWpret = dict_param["GW0"]
        self.GWabpret = 0
        self.seep_def = dict_param["seep_def"]
        self.w = dict_param["w"]
        self.vc = dict_param["vc"]
        self.IRC = dict_param["IRC"]
        self.h_dgw = dict_param["h_dgw"]
        self.down_seep = dict_param["down_seep"]
        self.soiltype = dict_param["soiltype"]
        self.croptype = dict_param["croptype"]
        self.dt = dict_param["dt"]
        self.soil_prm = soil_selector(soilmatrix, etmatrix, self.soiltype, self.croptype)
        
    def sol(self, UZ_per, Inf_p, OW_pret=0, IWU=0, IR_r=0, IR_p=0, IR_per=0):
        """
        Calculate the states and fluxes on grounodwater during current time step.
        
        Args: 
            UZ_per(float): percolation from unsaturated zone to groundwater [mm AUZ]
            Inf_p(float): infiltration from pavement to groundwater (only happen when PST = PIL) [mm AP]
            OW_pret(float): open water level at pevious time step [m-SL]
            BF_up(float): baseflow from upstream grid [L]
            IWU(float): indoor water usage (L)
            IR_r(float): irrigation on roof area [mm]
            IR_p(float): irrigation on pavement area [mm]
            IR_per(float): irrigation on pervious area [mm]
                        
        Returns: (dictionary): A dictionary of computed states and fluxes of grounodwater during current time step
            IRtot(float): total irrigation [L]
            LD(float): leakage depth of the current time step [mm AGW]
            Inflow_gw (float): percolation from unsaturated zone and infiltration from road to groundwater [mm AGW]
            sc_gw (float): storage coefficient of groundwater of the current time step [-]
            GWt0 (float): groundwater storage level in the beginning of the time step [m-SL]
            Q_seep (float): seepage from GW to deep GW [mm]  --> Q_seepL[L]
            BF_out (float): baseflow to open water (if negative: open water feeding GW) [mm]  --> BF*AGW[L]
            INFS(float): infiltration to wastewater pipe [mm], = max(0.0, (3-GWt)*IRC)
            GWt (float): groundwater storage level at current time step [m-SL]
            GWabt (float): groundwater level above surface level at current time step [m-SL]
            WB_gw (float): water balance in groundwater storage [L]
        """
        IRtot = IR_r * self.AR + IR_p * self.AP + IR_per * self.APer
        #LD (leakage depth = leakage rate * total input water (total IR + IWU))
        LD = (((IRtot + IWU) * self.LR)/(1-self.LR)) / self.AGW
        #INFS = np.sqrt(ht0 - UZ_per) * self.IRC
        Inflow_gw = (UZ_per * self.AUZ + Inf_p * self.AP) / self.AGW + LD  #+ BF_up

        gwl_sol = gwlcalc(self.GWpret)
        gw_up = gwl_sol[0]
        gw_low = gwl_sol[1]
        id1 = gwl_sol[2]
        id2 = gwl_sol[3]

        if self.GWpret < 10.0:
            sc_gw = self.soil_prm[id2]["stor_coef"] + (self.GWpret - gw_low) / (gw_up - gw_low) * (
                    self.soil_prm[id1]["stor_coef"] - self.soil_prm[id2]["stor_coef"])
        else:
            sc_gw = self.soil_prm[29]["stor_coef"]

        #calculate the GWt: Inflow_gw, BF, Q_seep, INFS
        if self.seep_def > 0.5:
            GWt0 = - ((Inflow_gw/1000 * self.w * self.vc - self.h_dgw * self.w 
                       - OW_pret * self.vc - 3.0 * self.IRC * self.w * self.vc) 
                      / (self.w + self.vc +self.IRC * self.w * self.vc)
                     + ( - (self.GWpret + self.GWabpret)
                        - (Inflow_gw/1000 * self.w * self.vc - self.h_dgw * self.w 
                           - OW_pret * self.vc - 3.0 * self.IRC * self.w * self.vc) 
                        / (self.w + self.vc + self.IRC * self.w * self.vc))
                     * np.exp(-self.dt * (self.w + self.vc + self.IRC * self.w * self.vc) 
                              / (sc_gw * self.w * self.vc)))
            Q_seep = 1000 * (self.h_dgw - 0.5 * (GWt0 + self.GWpret + self.GWabpret)) / self.vc * self.dt
            INFS = 1000 * (3.0 - 0.5 * (GWt0 + self.GWpret + self.GWabpret)) * self.IRC * self.dt
        else:
            GWt0 = - (((Inflow_gw - self.down_seep) /1000 * self.w 
                       - 3.0 * self.IRC * self.w - OW_pret) / (1 + self.IRC * self.w)
                     + ( - (self.GWpret + self.GWabpret) 
                        - ((Inflow_gw - self.down_seep) /1000 * self.w 
                           - 3.0 * self.IRC * self.w - OW_pret) / (1 + self.IRC * self.w))
                     * np.exp(-self.dt * (1 + self.IRC * self.w) / (sc_gw * self.w)))
            Q_seep = self.dt * self.down_seep
            INFS = 1000 * (3.0 - 0.5 * (GWt0 + self.GWpret + self.GWabpret)) * self.IRC * self.dt
        # BF = In - Out - dS
        BF = Inflow_gw - Q_seep - INFS - sc_gw * (self.GWpret + self.GWabpret / sc_gw - GWt0) * 1000
        #If (self.GWpret + self.GWabpret/sc_gw) > (Inflow_gw - Q_seep - BF)/(1000*sc_gw) --> GWt>0, GWabt=0
        GWt = max(0, self.GWpret + self.GWabpret / sc_gw - (Inflow_gw - Q_seep - INFS - BF) / (1000 * sc_gw))
        #If (self.GWpret + self.GWabpret/sc_gw) < (Inflow_gw - Q_seep - BF)/(1000*sc_gw) --> GWt=0, GWabt<0
        GWabt = -1 * max(0, (0 - (self.GWpret + self.GWabpret / sc_gw 
                                  - (Inflow_gw - Q_seep - INFS - BF) / (1000 * sc_gw)))) * sc_gw
        WB_gw = (Inflow_gw - Q_seep - INFS - BF 
                 - sc_gw * (self.GWpret + self.GWabpret / sc_gw - (GWt + GWabt / sc_gw)) * 1000) * self.AGW

        #if self.OW > 0.5:
        BF_out = BF * self.AGW
            #BF_down = 0.0
        #else:
            #BF_out = 0
            #BF_down = BF * self.AGW
        
        Q_seepL = Q_seep * self.AGW
        self.GWpret = GWt
        self.GWabpret = GWabt
        
        return {
            "IRtot": IRtot,
            "LD": LD,
            "Inflow_gw":Inflow_gw,
            "sc_gw":sc_gw,
            "GWt0":GWt0,
            "Q_seep": Q_seepL,
            "INFS": INFS,
            "BF_out": BF_out,
            "GWt": GWt,
            "GWabt": GWabt,
            "WB_gw": WB_gw,
        } 
        