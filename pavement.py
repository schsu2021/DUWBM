# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 10:03:44 2020

@author: ShuChen Hsu
"""
class Pavement:
    """
    calculate the status(storage) and flow with given area and storage capacity.
    Inflow: Rainfall(P), irrigation(IR_p), inflow from raintank(IRUN_rp)
    Outflow: Evaporation(Ea_p), 
             impervious surface runoff(IRUN_p), which fows to stormwater system,
             non-effective surface runoff(NEAR_p), which flows to pervious area
    Args: 
        AP(float): area of pavement area [m^2]
        EPA(float): effective pavement area [%]
        PIL(float): pavement area maximum initial loss (storage capacity) [mm]
        PST0(float): pavement area initial storage(at t=0) [mm]
        Infilc_p(float): infiltration (to GW) capacity of pavement [mm/d]
        dt(float): time step [day]
    """
    
    def __init__(self, dict_param):
        self.AP = dict_param["AP"]
        # if there is no pervious area, non-effective runoff also flows to stormwater system
        if dict_param["APer"] == 0:
            self.EPA = 1.0
        else:
            self.EPA = dict_param["EPA"]/100  #[%] --> [-]
        self.PIL = dict_param["PIL"]
        self.PSTpret = dict_param["PST0"]
        self.Infilc_p = dict_param["Infilc_p"]
        self.dt = dict_param["dt"]
        
    def sol(self, P, Ep, IRUN_rp, IR_p=0):
        """
        Calculate the states and fluxes on pavement area during current time step.
        
        Args: 
            P(float): Precipitation of the time step [mm]
            Ep(float): Potential evaporation of the time step [mm]
            IR_p(float): irrigation on paved area [mm]
            IRUN_rp(float): effective impervious surface runoff from raintank to pavement [L]
            
        Returns: (dictionary): A dictionary of computed states and fluxes of paved area during current time step
            Inflow_rp(float): A part of effective impervious surface runoff inflow to pavement [L--> mm AP] 
            PSTt0(float): pavement area interception storage after rainfall and irrigation at the beginning of the current time step [mm]
            Ea_p(float): evaporation from interception storage on pavement area during the current time step [mm]
            PSTt(float): pavement area interception storage level at the end of the current time step [mm]
            Inf_p(float): infiltration to groundwater (only happen when PST = PIL) [mm]
            IRUN_p(float): effective impervious surface runoff from pavement area [mm]
            NEAR_p(float): non-effective runoff from pavement area [mm]
            WB_p(float): water balance of pavement area [L]
        """
        if self.AP ==0:
            Inflow_rp = PSTt0 = Ea_p = PSTt = Inf_p = IRUN_p = NEAR_p = WB_p = 0.0
        else:
            Inflow_rp = IRUN_rp / self.AP
            PSTt0 = min(self.PIL, max(0.0, self.PSTpret + P + IR_p + Inflow_rp))
            Ea_p = min(Ep, PSTt0)
            PSTt = PSTt0 - Ea_p  #didn't abstract Inf_p because it is calculated as part of overflow
            Inf_p = max(0.0, 
                        min(P + IR_p + Inflow_rp - (self.PIL - self.PSTpret),
                            self.Infilc_p * self.dt))
            # (PIL - PSTpret) = available space in pavement, if P < (PIL - PSTpret) --> Inf_p = 0
            IRUN_p = self.EPA * max(0.0,
                                    (P + IR_p + Inflow_rp - Ea_p - Inf_p - (PSTt - self.PSTpret)))
            NEAR_p = max(0.0, 
                         (P + IR_p + Inflow_rp - Ea_p - Inf_p - (PSTt - self.PSTpret) - IRUN_p))

            WB_p = ((P + IR_p + Inflow_rp) - (Ea_p + Inf_p + IRUN_p + NEAR_p) - (PSTt - self.PSTpret)) * self.AP
        
        self.PSTpret = PSTt
        
        return {
            "Inflow_rp": Inflow_rp,
            "PSTt0": PSTt0,
            "Ea_p": Ea_p,
            "PSTt": PSTt,
            "Inf_p":  Inf_p,
            "IRUN_p": IRUN_p,
            "NEAR_p": NEAR_p,
            "WB_p": WB_p
        }
