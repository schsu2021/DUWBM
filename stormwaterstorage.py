# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 14:05:39 2020

@author: ShuChen Hsu
"""
class StormWaterStorage:
    """
    Args:
        AP, Aper (float): area [m^2]
        SWSop(boolean): is the stormwater storage (SWS) open for P and E? [-]
        ASWS(float): Area of SWS [m2]
        SWSc(float): SWS storage capacity [L]
        SWSff(float): predifined first flush of SWS [L]
        perI(float): percentage of runoff become inflow to wastewater system [%] 
    """
    
    def __init__(self, dict_param):
        self.AP = dict_param["AP"]
        self.APer = dict_param["APer"]
        self.SWSop = dict_param["SWSop"]
        self.ASWS = dict_param["ASWS"]
        self.SWSc = dict_param["SWSc"]
        self.SWSff = dict_param["SWSff"]
        self.perI = dict_param["perI"]/100   # [%] --> [-]
        
    def sol(self, P, Ep, SWSpret, IRUN_rrun, IRUN_p, EXC, Rs_up=0):
        """
        Calculate the states and fluxes on stormwater storage during current time step.
        
        Args:
            P(float): Precipitation of the time step [mm]
            Ep(float): Potential evaporation of the time step [mm]
            SWSpret(float): the storage volume in the end of previous time step [L]
            IRUN_rrun(float): Outflow from roof-raintank system [L]
            IRUN_p(float): effective impervious surface runoff from pavement area [mm AP]
            EXC(float): overflow from pervious area interception [mm APer]
            Rs_up(float): stormwater sewer system flow from upstream grid [L]
            
        Return: (dictionary): A dictionary of computed states and fluxes of stormwater storage during current time step
            RUN_tot(float): total runoff [L]
            ISI(float): inflow of stormwater to wastewater sewer system [L]
            RUN(float): runoff to stormwater system [L]
            SWSff_a (float): actual first flush [L]
            Rs_up (float): stormwater sewer system flow from upstream grid [L]
            In_sws(float): inflow to SWS, RUN - first flush + Precipitation if SWSop=1 [L]
            SWSt0(float): storage volume in the beginning of the time step [L]
            Ea_sws(float): evaporation volume if SWSop=1 [L]
            SWSt(float): SWS storage in the end of the time step[L]
            EXC_sws(float): overflow from SWS [L]
            Rs(float): stormwater sewer system inflow [L]
            WB_sws(float): water balance of SWS[L]
        """
        RUN_tot = IRUN_rrun + IRUN_p * self.AP + EXC * self.APer + Rs_up
        ISI = self.perI * RUN_tot
        RUN = RUN_tot - ISI
            
        if self.SWSc ==0:
            SWSff_a = In_sws = SWSt0 = Ea_sws = SWSt = EXC_sws = 0.0
            Rs = RUN
        else:
            SWSff_a = min(RUN, self.SWSff)
            # self.RTop = 1 (open for rainfall)
            In_sws = RUN - SWSff_a + self.SWSop * P * self.ASWS  # + Rs_up
            SWSt0 = min(self.SWSc, max(0.0, SWSpret + In_sws))
            Ea_sws = self.SWSop * min(Ep * self.ASWS, SWSt0)
            SWSt = SWSt0 - Ea_sws
            EXC_sws = max(0.0, In_sws - Ea_sws - (SWSt - SWSpret))
            Rs = SWSff_a + EXC_sws
        
        WB_sws = In_sws - (Ea_sws + EXC_sws) - (SWSt - SWSpret)
        
        #self.SWSpret = SWSt
        
        return{
            "RUN_tot": RUN_tot,
            "ISI": ISI,
            "RUN": RUN,
            "SWSff_a": SWSff_a,
            "Rs_up": Rs_up,
            "In_sws": In_sws,
            "SWSt0": SWSt0,
            "Ea_sws": Ea_sws,
            "SWSt": SWSt,
            "EXC_sws": EXC_sws,
            "Rs": Rs,
            "WB_sws": WB_sws,
        }