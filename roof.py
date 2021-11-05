# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 10:00:57 2020

@author: ShuChen Hsu
"""

class Roof:
    """
    calculate the status(storage) and flow with given area and storage capacity.
    Inflow: Rainfall(P), irrigation(IR_r)
    Outflow: Evaporation(Ea_r), 
             impervious surface runoff(IRUN_r), which fows to RainTank,
             non-effective surface runoff(NEAR_r), which flows to pervious area

    Args: 
        AR(float): area of roof [m^2]
        ERA(float): effective roof area (connected with gutter) [%]
        RIL(float): Roof area maximum initial loss (storage capacity) [mm]
        RST0(float): roof initial storage(t0) [mm]
        dt(float): time step [day]
    """
    def __init__(self, dict_param):
        self.AR = dict_param["AR"]
        self.APer = dict_param["APer"] 
        # if there is no pervious area, non-effective runoff also flows to gutter
        if dict_param["APer"] == 0:
            self.ERA = 1.0
        else:
            self.ERA = dict_param["ERA"]/100  #[%] --> [-]
        self.RIL = dict_param["RIL"]
        self.RSTpret = dict_param["RST0"]
        self.dt = dict_param["dt"]

    def sol(self, P, Ep, IR_r=0):
        """
        Calculate the states and fluxes on roof during current time step.

        Args: 
            P(float): Precipitation of the time step [mm]
            Ep(float): Potential evaporation of the time step [mm]
            IR_r(float): irrigation on roof area[mm]

        Returns: (dictionary): A dictionary of computed states and fluxes of roof during current time step
            IR_r(float): irrigaiton [mm AR]
            RSTt0(float): roof interception storage after rainfall and irrigation at the beginning of the current time step [mm AR]
            Ea_r(float): evaporation from interception storage on roof during the current time step [mm AR]
            RSTt(float): roof interception storage level at the end of the current time step [mm AR]
            IRUN_r(float): effective impervious surface runoff (collected to rain tank(if exist),  storm sewer and pavement) [mm AR]
            NEAR_r(float): non-effective runoff flows to pavement and pervious area [mm AR]
            WB_r(float): water balance of roof [L]= In - Out - dS
            
        """
        
        if self.AR ==0:
            RSTt0 = Ea_r = RSTt = IRUN_r = NEAR_r = WB_r = 0.0
        else:
            RSTt0 = min(self.RIL, max(0.0, self.RSTpret + P + IR_r))
            Ea_r = min(Ep, RSTt0)
            RSTt = RSTt0 - Ea_r

            IRUN_r = self.ERA * max(0.0, (P + IR_r - Ea_r - (RSTt - self.RSTpret)))
            NEAR_r = max(0.0, (P + IR_r - Ea_r - (RSTt - self.RSTpret) - IRUN_r))

            WB_r = ((P + IR_r) - (Ea_r + IRUN_r + NEAR_r) - (RSTt - self.RSTpret)) * self.AR

        self.RSTpret = RSTt

        return {
            "IR_r": IR_r,
            "RSTt0": RSTt0,
            "Ea_r": Ea_r,
            "RSTt": RSTt,
            "IRUN_r": IRUN_r,
            "NEAR_r": NEAR_r,
            "WB_r": WB_r
        }