# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 10:43:11 2020

@author: ShuChen Hsu
"""
class RainTank:
    """
    Args:
        RTop(boolean): is the raintank open for P and E? [-]
        ART(float): area of one Rain tank [m2]
        RTc(float): raintank storage capacity of one Rain tank [L]
        RTff(float): predifined first flush of one raintank [L]
        RT0(float): initial storage of raintank [L]
        ERA_out(float): effective roof area (connected with gutter) flow to pavement [%]
        pRT(float): % of houses(AR) installing Raintank [%]
            (NH * pRT = number of raintank)
        NH(float): Number of houses in the cell [-]
    """
    
    def __init__(self, dict_param):
        self.RTop = dict_param["RTop"]
        self.pRT = dict_param["pRT"]/100              #[%] --> [-]
        #self.Aoccu = dict_param["Aoccu"]
        self.NH = dict_param["NH"]
        self.ART = dict_param["ART"] * self.NH * self.pRT
        self.RTc = dict_param["RTc"] *self.NH * self.pRT
        self.RTff = dict_param["RTff"] * self.NH * self.pRT
        #self.RTpret = dict_param["RT0"]
        if dict_param["AP"] == 0:
            self.ERA_out = 1.0
        else:
            self.ERA_out = dict_param["ERA_out"]/100  #[%] --> [-]
        self.AR = dict_param["AR"]
        
    def sol(self, P, Ep, RTpret, IRUN_r):
        """
        Calculate the states and fluxes on raintank during current time step.
        
        Args:
            P(float): Precipitation of the time step [mm]
            Ep(float): Potential evaporation of the time step [mm]
            RTpret(float): raintank storage in the end of previous time step[L]
            IRUN_r(float): effective impervious surface runoff (collected to rain tank(if exist),  storm sewer and pavement) [mm AR]
            
        Return: (a dictionary)
            RTff_a (float): actual first flush [L]
            In_rt(float): inflow to raintank [L]
            RTt0(float): storage volume in the beginning of the time step [L]
            Ea_rt(float): evaporation volume if RTop=1 [L]
            RTt(float): raintank storage in the end of the time step[L]
            EXC_rt(float): overflow from raintank [L]
            Out_r(float): outflow from roof-raintank system [L]
            WB_rt(float): water balance of raintank[L]
            IRUN_rrun(float): effective impervious surface runoff to storm sewer [L]
            IRUN_rp(float): effective impervious surface runoff to pavement [L]
        """
        if self.RTc ==0:
            RTff_a = In_rt = RTt0 = Ea_rt = RTt = EXC_rt = WB_rt = 0.0
            Out_r = IRUN_r * self.AR
        else:
            RTff_a = min(IRUN_r * self.AR * self.pRT, self.RTff)
            # self.RTop = 1 (open for rainfall)
            In_rt = IRUN_r * self.AR * self.pRT - RTff_a + self.RTop * P * self.ART
            RTt0 = min(self.RTc, max(0.0, RTpret + In_rt))
            Ea_rt = self.RTop * min(Ep * self.ART, RTt0)
            RTt = RTt0 - Ea_rt
            EXC_rt = max(0.0, In_rt - Ea_rt - (RTt - RTpret))
            Out_r = RTff_a + EXC_rt + IRUN_r * self.AR * (1.0 - self.pRT)
            WB_rt = In_rt - (Ea_rt + EXC_rt) - (RTt - RTpret)
        
        IRUN_rrun = self.ERA_out * Out_r
        IRUN_rp = Out_r - IRUN_rrun
        #self.RTpret = RTt
        
        return{
            "RTff_a": RTff_a,
            "In_rt": In_rt,
            "RTt0": RTt0,
            "Ea_rt": Ea_rt,
            "RTt": RTt,
            "EXC_rt": EXC_rt,
            "Out_r": Out_r,
            "WB_rt": WB_rt,
            "IRUN_rrun": IRUN_rrun,
            "IRUN_rp": IRUN_rp
        }