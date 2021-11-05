# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:14:23 2020

@author: ShuChen Hsu
"""

from DUWBM.selector import soil_selector, et_selector
from DUWBM.gwlcalculator import gwlcalc

class UnsaturatedZone:
    """
    Args:
        AUZ (float): area of unsaturated zone = APer [m^2]
        dt(float): time step [day]
        soiltype (int): soil type
        croptype (int): crop type
        h0 (float): initial (volumetric) moisture content of soil in the root zone(assume = h2) (at t=0) [mm]
        h3l (float): equilibrium moisture content of soil in root zone, transpiration (E_pot ≤ 1 mm/d) reduction starts
        h3h (float): equilibrium moisture content of soil in root zone, transpiration (E_pot ≥ 5 mm/d) reduction starts
        h1 (float): equilibrium moisture content of soil in root zone, groundwater level at surface level, i.e. complete saturation
        h2 (float): equilibrium moisture content of soil in root zone, groundwater level at bottom root zone, i.e. field capacity
        h4 (float): equilibrium moisture content of soil in root zone, transpiration = 0, i.e. permanent wilting point
        
    """
    def __init__(self, dict_param, soilmatrix, etmatrix):
        self.AUZ = dict_param["AUZ"]
        self.dt = dict_param["dt"]
        self.hpret = dict_param["h0"]
        self.soiltype = dict_param["soiltype"]
        self.croptype = dict_param["croptype"]
        #select soil parameter matrix depend on soil type and crop type
        self.soil_prm = soil_selector(soilmatrix, etmatrix, self.soiltype, self.croptype)
        self.UZKsat = 10 * self.soil_prm[0]["k_sat"]
        # self.et_prm (dataframe): a matrix of root zone - related parameters
        self.et_prm = et_selector(etmatrix, self.soiltype, self.croptype)
        self.h3l = self.et_prm["theta_h3l_mm"].values[0]
        self.h3h = self.et_prm["theta_h3h_mm"].values[0]
        self.h1 = self.et_prm["theta_h1_mm"].values[0]
        self.h2 = self.et_prm["theta_h2_mm"].values[0]
        self.h4 = self.et_prm["theta_h4_mm"].values[0]
        
    def sol(self, Eref, Inf_per, GWpret):
        """
        Calculate the states and fluxes on unsaturated area during current time step.
        
        Args: 
            Eref(float): Reference evaporation of the time step [mm]
            Inf_per(float): infiltration from pervious area[mm APer]
            GWpret(float): groundwater level at previous time step [m-SL]
                        
        Returns: (dictionary): A dictionary of computed states and fluxes of unsaturated zone durning current time step
            h3_uz(float): equilibrium moisture content of soil in root zone where transpiration reduction starts [mm]
            alpha(float): transpiration reduction factor during current time step [-]
            ET(float): actual transpiration from unsaturated zone [mm]
            gw_up(float): groundwater level in the lookup table (upper one) [mm]
            gw_low(float): groundwater level in the lookup table (lower one) [mm]
            heq(float): equilibrium moisture content of soil of the current time step [mm]
            Cmax(float): maximum capillary rise rate [mm/d]
            ht0(float): moisture content of root zone at the beginning of the time step [mm]
            UZ_per(float): percolation from UZ to GW [mm] (+:percolation, -:capillary rise)
            ht(float): moisture content of root zone at the end of the time step [mm]
            WB_uz(float): water balance in unsaturated zone [L]
        """
        #Inf_per = Inf_per * self.APer / self.AUZ (APer = AUZ)
        if self.AUZ ==0:
            h3_uz = alpha = ET = gw_up = gw_low = heq = Cmax = UZ_per = WB_uz = 0.0
            ht0 = ht = self.hpret
        else:
            if Eref < 1.0:
                h3_uz = self.h3l
            elif Eref > 5.0:
                h3_uz = self.h3h
            else:
                h3_uz = (Eref - 1.0)/(5.0 - 1.0) * (self.h3h - self.h3l) + self.h3l

            if self.hpret + Inf_per > self.h1:
                alpha = 0.0
            elif self.hpret + Inf_per > self.h2:
                alpha = 1 - (self.hpret + Inf_per - self.h2)/(self.h1 - self.h2)
            elif self.hpret + Inf_per > h3_uz:
                alpha = 1.0
            elif self.hpret + Inf_per > self.h4:
                alpha = (self.hpret + Inf_per - self.h4)/(h3_uz - self.h4)
                #alpha = 1 - (self.hpret + Inf_per - self.h4)/(h3_uz - self.h4)
            else:
                alpha = 0.0

            ET = alpha * Eref
        
            gwl_sol = gwlcalc(GWpret)
            gw_up = gwl_sol[0]
            gw_low = gwl_sol[1]
            id1 = gwl_sol[2]
            id2 = gwl_sol[3]
            # heq and Cmax of the (soiltype+croptype) depend on groundwater level 
            if GWpret < 10.0:
                # for GW level not deeper than 10 m: interpolation
                heq = self.soil_prm[id1]["moist_cont_eq_rz[mm]"] + (GWpret - gw_up)/(gw_low - gw_up) * (
                    self.soil_prm[id2]["moist_cont_eq_rz[mm]"] - self.soil_prm[id1]["moist_cont_eq_rz[mm]"])
    
                Cmax = self.soil_prm[id1]["capris_max[mm/d]"] + (GWpret - gw_up)/(gw_low - gw_up) * (
                    self.soil_prm[id2]["capris_max[mm/d]"] - self.soil_prm[id1]["capris_max[mm/d]"])
            else:
                # for GW level deeper than 10 meter: use the values at 10 m-SL
                heq = self.soil_prm[29]["moist_cont_eq_rz[mm]"]
                Cmax = self.soil_prm[29]["capris_max[mm/d]"]
    
            ht0 = self.hpret + Inf_per - ET
            if ht0 > heq:
                UZ_per = min(ht0 - heq, self.dt * self.UZKsat)
            else:
                UZ_per = -1 * min((heq-ht0), self.dt * Cmax)
    
            #INFS = np.sqrt(ht0 - UZ_per) * self.IRC
            ht = ht0 - UZ_per #- INFS
            WB_uz = (Inf_per - (ET + UZ_per) - (ht - self.hpret)) * self.AUZ   # (ET + UZ_per + INFS)
        
        self.hpret = ht
        
        return {
            "h3_uz": h3_uz,
            "alpha": alpha,
            "ET": ET,
            "gw_up": gw_up,
            "gw_low": gw_low,
            "heq": heq,
            "Cmax": Cmax,
            "ht0": ht0,
            "UZ_per": UZ_per,
            "ht": ht,
            "WB_uz": WB_uz,
        } 