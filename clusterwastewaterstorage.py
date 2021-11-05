# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 14:22:51 2020

@author: ShuChen Hsu
"""

class ClusterWasteWaterStorage:
    """
    Args:
        AcWWS(float): Area of WWS [m2]
        cWWSc(float): WWS storage capacity [L]
        cWWS0(float): initial storage of WWS [L]
        AGW (float): area of groundwater storage = total area [m^2] 
    """
    
    def __init__(self, dict_param):
        self.AGW = dict_param["AGW"]
        self.AcWWS = dict_param["AcWWS"]
        self.cWWSc = dict_param["cWWSc"]
        #self.cWWSpret = dict_param["cWWS0"]
        
    def sol(self, cWWSpret, WWSsp, INFS, ISI, Rw_up=0):
        """
        Calculate the states and fluxes on wastewater storage during current time step.
        
        Args:
            cWWSpret(float): input the updated raintank after WaterReuse
            WWSsp(float):outflow from onsite wastewater storage [L]
            INFS: infiltration from UZ to wastewater pipe [mm AGW]
            ISI(float): inflow of stormwater to wastewater sewer system[L]
            Rw_up(float): wastewater sewer system flow from upstream grid [L]
            
        Return:(dictionary): A dictionary of computed states and fluxes of wastewater storage during current time step
            In_cwws(float): inflow to WWS [L]
            cWWSt(float): WWS storage in the end of the time step[L]
            WB_cwws(float): water balance of WWS[L]
            Rw_up(float): wastewater sewer system flow from upstream grid [L]
            Rw(float): wastewater sewer system inflow [L]
        """
        RUNw = WWSsp + INFS * self.AGW + ISI + Rw_up
        
        if self.cWWSc ==0:
            In_cwws = cWWSt = WB_cwws = 0.0
            Rw = RUNw
        else:
            In_cwws = RUNw
            cWWSt = min(self.cWWSc, cWWSpret + In_cwws)
            Rw = max(0.0, In_cwws - (cWWSt - cWWSpret))  #EXC_cwws
            WB_cwws = In_cwws - Rw - (cWWSt - cWWSpret)
        
        #self.cWWSpret = cWWSt
        
        return{
            "In_cwws": In_cwws,
            "cWWSt": cWWSt,
            "WB_cwws": WB_cwws,
            "Rw_up": Rw_up,
            "Rw": Rw,
        }