# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 09:36:04 2020

@author: ShuChen Hsu
"""

#import numpy as np

class CellReuse:
    """
    Calculate the total IWU (indoor water use) according to occupancy and number of houses in the cell.
    
    Args:
        NH (float): nuber of houses in this cell
        demand(series): the water demand of different indoor water use per block [L/day/block]
        setreuse(dataframe): setting of the supply and use priority (boolean numbers)
        AWWS(float): Area of onsite WWS [m2]
        WWSc(float): onsite WWS storage capacity [L]
        WWS0(float): initial storage of onsite WWS [L]
        IWU (float): daily indoor water use
    """
    
    def __init__(self, demand, setreuse, dict_param):
        self.NH = dict_param["NH"]
        #demand_pop = demand * self.NH
        demand_pop = demand
        #supply of SSG
        self.SSGS = (setreuse.KforSSG * demand_pop["K"].values[0]
                      + setreuse.BforSSG * demand_pop["B"].values[0]
                      + setreuse.LforSSG * demand_pop["L"].values[0])
        #Demand from rain tank (Which water use can raintank suppply)
        self.RTD0 = (setreuse.RTforK * demand_pop["K"].values[0]
                      + setreuse.RTforB * demand_pop["B"].values[0]
                      + setreuse.RTforL * demand_pop["L"].values[0])
        self.K0 = demand_pop["K"].values[0]
        self.B0 = demand_pop["B"].values[0]
        self.L0 = demand_pop["L"].values[0]
        self.T0 = demand_pop["T"].values[0]
        self.IWU = dict_param["IWU"]
        self.setreuse = setreuse
        self.pRT = dict_param["pRT"]/100.0
        self.AR = dict_param["AR"]
        self.AP = dict_param["AP"]
        self.APer = dict_param["APer"]
        self.AGW = dict_param["AGW"]
        self.RTc = dict_param["RTc"] *self.NH * self.pRT
        self.AWWS = dict_param["AWWS"] * self.NH
        self.WWSc = dict_param["WWSc"] * self.NH
        self.WWSpret = dict_param["WWS0"] * self.NH
        
    def sol(self, IR_r, IR_p, IR_per, RTt, WB_rt, LD):
        """
        Simulate the supply(S)/ demand(D)/ use(use)/ deficit(def) of SSG, wws, and raintank.
        Update the water level and water balance of raintank.
        
        Args:
            IR_r, IR_p, IR_per: irrigation demand on the surface [mm r/p/per]
            RTt: rain tank amount in this time step from the RainTank object [L]
            WB_rt: water balance check of RT (form RainTank object) [L]
            LD: leakage rate [mm AGW]
        Return:
            RTt, RTtend: real water content in RT after fill/spill/consumption [L]
             (the RTt from RainTank will be replace by the value derived here at every time step)
            Import(float): required imported water = KBL1 + T2 + IR3 + LD [L]
        """
        IRt = IR_r*self.AR + IR_p*self.AP + IR_per*self.APer
        # IR0 = IR(t)
        if self.SSGS ==0:
            SSGD = SSGsp = SSGuse = SSGdef = 0
            IR1 = IRt
        else:
            SSGD = IRt #input from hydrological scheme
            SSGuse = min(self.SSGS, SSGD)  #SSGS = SSG supply/ availble amount
            SSGsp = max(0, self.SSGS - SSGuse)
            SSGdef = max(SSGD - SSGuse,0)
            IR1 = SSGdef
        checkssg = (IRt - IR1 == SSGuse)
        
        # onsite wastewater storage
        WWSin = self.T0 + SSGsp + ((1-self.setreuse.KforSSG) * self.K0
                                   + (1-self.setreuse.BforSSG) * self.B0
                                   + (1-self.setreuse.LforSSG) * self.L0)
        if self.WWSc ==0:
            WWSD = WWSS = WWSuse = WWSt = WWSdef = WB_wws = 0.0
            WWSsp = WWSin
            T1 = self.T0
            IR2 = IR1
        else:
            WWSD = self.T0 * self.setreuse.WWSforT + IR1 * self.setreuse.WWSforIR
            if WWSD ==0:
                WWSS = WWSuse = WWSdef = WB_wws = 0
                WWSt = min(self.WWSpret + WWSS, self.WWSc)
                WWSsp = max(self.WWSpret + WWSS- self.WWSc, 0)
                T1 = self.T0
                IR2 = IR1
            else:
                WWSS = WWSin
                WWSt0 = min(self.WWSpret + WWSS, self.WWSc)
                WWSuse = min(WWSt0, WWSD)
                WWSt = WWSt0 - WWSuse
                WWSdef = max(WWSD - WWSuse, 0)
                WWSsp = max(self.WWSpret + WWSS- self.WWSc, 0)
                WB_wws = WWSS - (WWSsp + WWSuse) - (WWSt - self.WWSpret)
                if (self.setreuse.WWSforT ==0) & (self.setreuse.WWSforIR !=0):
                    T1 = self.T0
                    IR2 = IR1 - WWSuse
                elif (self.setreuse.WWSforT !=0) & (self.setreuse.WWSforIR ==0):
                    T1 = self.T0 - WWSuse
                    IR2 = IR1
                else:
                    if WWSuse < self.T0:
                        T1 = self.T0 - WWSuse
                        IR2 = IR1
                    else:
                        T1 = 0
                        IR2 = IR1 - (WWSuse - self.T0)
        self.WWSpret = WWSt
        checkwws = ((self.T0 - T1) + (IR1 - IR2) == WWSuse)
                
        # raintank
        
        if self.RTc ==0:
            RTD = RTuse = RTdef = 0
            RTtend = RTt
            KBL1 = (self.K0 + self.B0 + self.L0)
            T2 = T1
            IR3 = IR2
        else:
            RTD = self.RTD0 + T1 * self.setreuse.RTforT + IR2 * self.setreuse.RTforIR
            if RTD == 0:
                RTuse = RTdef = 0
                RTtend = RTt
                KBL1 = (self.K0 + self.B0 + self.L0)
                T2 = T1
                IR3 = IR2
            else:
                #RTavailable = RTt
                RTuse = min(RTt, RTD)
                RTtend = RTt - RTuse
                RTdef = max(RTD - RTuse, 0)
                #RTuse: sth left after satisfying RTD0 /all for KBL 
                KBL1 = (self.K0 + self.B0 + self.L0)-min(self.RTD0, RTuse)
                if self.setreuse.RTforT == 0:
                    T2 = T1
                    IR3 = IR2 - max(RTuse - self.RTD0, 0)
                else:
                    T2 = T1 - min(max(RTuse - self.RTD0,0), T1)
                    IR3 = IR2 - max(RTuse - self.RTD0 - T1,0)
        checkrt = ((self.K0 + self.B0 + self.L0) - KBL1) + (T2-T1) + (IR3 - IR2) == RTuse
        newWB_rt = WB_rt - RTuse + (RTt-RTtend)
        Import = KBL1 + T2 + IR3 + LD*self.AGW
        return {
                "SSGD": SSGD, "SSGuse": SSGuse, "SSGsp":SSGsp, "SSGdef": SSGdef,"WWSD":WWSD,
                "WWSS": WWSS, "WWSuse": WWSuse, "WWSt": WWSt, "WWSdef": WWSdef, "WWSsp":WWSsp,"WB_wws":WB_wws,
                "RTD":RTD, "RTuse":RTuse, "RTtend": RTtend, "RTdef": RTdef,"RTt":RTtend, "Import": Import,
                "KBL": KBL1, "T": T2, "IR": IR3, "checkssg":checkssg, "checkwws": checkwws, "checkrt":checkrt,
                "WB_rt":newWB_rt, 'cWWSuse':0, 'SWSuse':0, 'cWWSsup':0, 'SWSsup':0 
            }