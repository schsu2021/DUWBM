# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 10:30:31 2020

@author: ShuChen Hsu
"""
from DUWBM.selector import soil_selector

class Pervious:
    """
    calculate the status(storage) and flow with given area and storage capacity.
    Inflow: Rainfall(P), irrigation(IR_per), inflow from roof and pavement(NEAR_r, NEAR_p)
    Outflow: Evaporation(Ea_per), overflow from pervious interception(EXC), 
             infiltration to GW (Inf_per)
    Args: 
        APer(float): area of pervious area [m^2]
        AR(float): area of roof [m^2]
        AP(float): area of pavement [m^2]
        dt(float):time step [d]
        PerIL(float): pervious area maximum initial loss (storage capacity) [mm]
        PerST0(float): pervious area initial storage(at t=t0) [mm]
        Infilc_per: predefined infiltration capacity of pervious area [mm/d]
        soiltype(int): soil type
        croptype(int): crop type
        UZc(float): maximum moisture content of root zone [mm]
        UZKsat(float): predefined saturated permeability of soil [mm/d]
    """
    
    def __init__(self, dict_param, soilmatrix, etmatrix):
        self.APer = dict_param["APer"]
        self.AR = dict_param["AR"]
        self.AP = dict_param["AP"]
        self.dt = dict_param["dt"]
        self.PerIL = dict_param["PerIL"]
        self.PerSTpret = dict_param["PerST0"]
        self.Infilc_per = dict_param["Infilc_per"]
        self.soiltype = dict_param["soiltype"]
        self.croptype = dict_param["croptype"]
        #select soil parameter matrix [at gwl=0] depend on soil type and crop type
        self.soil_prm = soil_selector(soilmatrix, etmatrix, self.soiltype, self.croptype)[0]
        self.UZc = self.soil_prm["moist_cont_eq_rz[mm]"]
        self.UZKsat = 10 * self.soil_prm["k_sat"]
        
    def sol(self, P, Ep, NEAR_r, NEAR_p, UZpret, IR_per=0):
        """
        Calculate the states and fluxes on paved area during current time step.
        
        Args: 
            P(float): Precipitation of the time step [mm]
            Ep(float): Potential evaporation of the time step [mm]
            NEAR_r(float): non-effective runoff of roof [mm AR]
            NEAR_p(float): non-effective runoff of pavement [mm AP]
            UZpret(float): moisture content in unsaturated zone of previous time step [mm]
            IR_per(float): irrigation on pervious area [mm]
                        
        Returns: (dictionary): A dictionary of computed states and fluxes of pervious area durning current time step
            Inflow_per(float): inflow from impervious area [mm (APer)]
            PerSTt0(float): pervious area interception storage after inflow, rainfall and irrigation at the beginning of the current time step [mm]
            Infilc_pera(float): actual infiltration capacity on pervious area[mm]
            timef_per(float): time factor, i.e. part of time step that interception storage on pervious area is available [mm]
            Ea_per(float): evaporation from interception storage on pervious area during the current time step [mm]
            Inf_per(float): infiltration from pervious area to unsaturaed zone [mm]
            PerSTt(float): pervious area interception storage level at the end of the current time step [mm]
            EXC(float): overflow from pervious interception [mm]
            WB_per(float): water balance of pervious area [L]
        """
        if self.APer ==0:
            Inflow_per = PerSTt0 = Infilc_pera = timef_per = Ea_per = Inf_per = PerSTt = EXC = WB_per = 0.0
        else:
            Inflow_per = (NEAR_r * self.AR + NEAR_p * self.AP) /self.APer
            PerSTt0 = max(0.0, self.PerSTpret + P + Inflow_per + IR_per)  #Not limited by GIL
            Infilc_pera = min(self.dt * self.Infilc_per, 
                           self.UZc - UZpret + min(self.UZc - UZpret, self.dt * self.UZKsat))
            timef_per = min(1.0, PerSTt0 / (Ep + Infilc_pera))
            Ea_per = timef_per * Ep
            Inf_per = timef_per * Infilc_pera
            PerSTt = min(self.PerIL, max(0.0, PerSTt0 - Ea_per - Inf_per))
            EXC = max(0.0, (P + Inflow_per + IR_per) - (Ea_per + Inf_per) - (PerSTt - self.PerSTpret))
            WB_per = ((P + Inflow_per + IR_per) - (Ea_per + Inf_per + EXC) - (PerSTt - self.PerSTpret)) * self.APer
        
        self.PerSTpret = PerSTt
       
        return {
            "Inflow_per":Inflow_per,
            "PerSTt0":PerSTt0,
            "Infilc_pera":Infilc_pera,
            "timef_per":timef_per,
            "Ea_per":Ea_per,
            "Inf_per":Inf_per,
            "PerSTt":PerSTt,
            "EXC":EXC,
            "WB_per":WB_per,
        }