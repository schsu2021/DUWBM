# -*- coding: utf-8 -*-
"""
Created on Sat May  2 10:48:24 2020

@author: ShuChen Hsu
"""

import numpy as np
import pandas as pd
#from random import seed
from random import choice
from simpledbf import Dbf5
from DUWBM.selector import soil_selector
from DUWBM.gwlcalculator import gwlcalc
from DUWBM.findorderUB import findordererUB
from DUWBM.nearest import NearestDownOW
from DUWBM.roof import Roof
from DUWBM.raintank import RainTank
from DUWBM.pavement import Pavement
from DUWBM.pervious import Pervious
from DUWBM.unsaturatedzone import UnsaturatedZone
from DUWBM.groundwater import Groundwater
from DUWBM.stormwaterstorage import StormWaterStorage
from DUWBM.cellreuse import CellReuse
from DUWBM.clusterwastewaterstorage import ClusterWasteWaterStorage
from DUWBM.waterbalancechecker import WaterBalanceCheck
from DUWBM.waterbalancechecker_all import WaterBalanceCheckall


class NewModel(object):
    """
    Running only the block scale to check water balance
    """
    
    def __init__(self, dict_param, soilmatrix, etmatrix, demand, setreuse):
        self.param = dict_param
        self.roof = Roof(self.param)
        self.raintank = RainTank(self.param)
        self.pavement = Pavement(self.param)
        self.pervious = Pervious(self.param, soilmatrix, etmatrix)
        self.unsaturatedzone = UnsaturatedZone(self.param, soilmatrix, etmatrix)
        self.groundwater = Groundwater(self.param, soilmatrix, etmatrix)
        self.stormws = StormWaterStorage(self.param)
        self.reuse = CellReuse(demand, setreuse, dict_param)
        self.cwastews = ClusterWasteWaterStorage(self.param)
        #self.openwater = OpenWater(self.param)
        
    def __iter__(self):
        return self
    
    def __next__(self, P, Ep, OWpret, IR, Eref, lst_pret, lstup):
        """
        Calculate storages, fluxes at current time step.
        """
        try:
            r_sol = self.roof.sol(P=P, Ep=Ep, IR_r=IR["IR_r"])
            rt_sol = self.raintank.sol(P=P, Ep=Ep, RTpret = lst_pret["RTt"], IRUN_r=r_sol["IRUN_r"])
            p_sol = self.pavement.sol(P=P, Ep=Ep, IRUN_rp=rt_sol["IRUN_rp"], IR_p=IR["IR_p"])
            per_sol = self.pervious.sol(P=P, Ep=Ep, NEAR_r=r_sol["NEAR_r"], NEAR_p=p_sol["NEAR_p"],
                                        UZpret=lst_pret["ht"], IR_per=IR["IR_per"])
            uz_sol = self.unsaturatedzone.sol(Eref=Eref, Inf_per=per_sol["Inf_per"],GWpret=lst_pret["GWt"])
            gw_sol = self.groundwater.sol(UZ_per=uz_sol["UZ_per"], Inf_p=p_sol["Inf_p"],OW_pret=OWpret,
                                          IWU=self.reuse.IWU, IR_r=IR["IR_r"], IR_p=IR["IR_p"], IR_per=IR["IR_per"])
            sws_sol = self.stormws.sol(P=P, Ep=Ep, SWSpret=lst_pret["SWSt"], IRUN_rrun=rt_sol["IRUN_rrun"],
                                       IRUN_p=p_sol["IRUN_p"], EXC=per_sol["EXC"], Rs_up=lstup["Rs"])
            reuse_sol = self.reuse.sol(IR_r=IR["IR_r"], IR_p=IR["IR_p"], IR_per=IR["IR_per"], 
                                       RTt=rt_sol["RTt"], WB_rt=rt_sol["WB_rt"], LD=gw_sol["LD"])
            wws_sol = self.cwastews.sol(cWWSpret=lst_pret["cWWSt"], WWSsp=reuse_sol["WWSsp"], INFS=gw_sol["INFS"],
                                        ISI=sws_sol["ISI"],Rw_up=lstup["Rw"])
            #ow_sol = self.openwater.sol(P=P, Ep=Ep, Rs=sws_sol["Rs"], BF=gw_sol["BF"])
            merged_dict = {**r_sol, **rt_sol, **p_sol, **per_sol,
                           **uz_sol, **gw_sol, **sws_sol, **reuse_sol, **wws_sol}   #OrderedDict()
        except IndexError:
            raise StopIteration
        return merged_dict
    
def read_param2(UrbanB, par, size, altwater, groundwater, soilfile, etfile, direction, w_k=5.0, w_D=0.13):
    """
    Args:
        UrbanB(string): UrbanBEATS parameters data file name.
        par(string): calibrated parameters data file name. When SA, par(df)
        size(float): cell size (side length) [m]
        altwater(string): the dimension of alternative water device 
                        (raintank, onsite WWS, SWS, cluster-scale WWS) data file name.
        groundwater (string): initial groundwater level (in the begining of simulation) csv file name 
                        (ID, GW0 [m-SL], h_dgw [m-SL])
        soilfile (string): soilparameter data file name.
        etfile (string): etparameter data file name.
        direction(int): Direction that is considered for neighborhoods. (options: 4, 6, 8)
    Return:
        paramlist()
    """
    if type(par)==str:
        par = pd.read_csv(par,header=None,index_col = 0)
    dbf = Dbf5(UrbanB, codec='utf-8')
    df = dbf.to_dataframe()
    df.set_index('BlockID', inplace=True)
    ##
    downOWids, downOWds = NearestDownOW(df, direction, size)
    ##
    alternativew = pd.read_csv(altwater)
    #alternativew.loc[0] = np.zeros(len(alternativew.iloc[0,:]))
    alternativew.loc[len(alternativew)] = np.zeros(len(alternativew.columns))  #new
    alternativew.set_index('id', inplace=True)                                   #new
    GW = pd.read_csv(groundwater) #[m-SL]
    GW.set_index('BlockID', inplace=True)
    soildata = pd.read_csv(soilfile, header = 0)
    etdata = pd.read_csv(etfile, header =0)
    paramlist = []
    names = locals()
    i = 1
    OWi = 0
    for idarea in df.index:
        cellid = idarea
        if direction==6:
            Atotal = df.Active[idarea]*(1.5*np.sqrt(3)*(size**2))
        else:
            Atotal = df.Active[idarea]*(size**2)
        AR = df.Blk_RoofsA[idarea]
        AP = df.Blk_TIA[idarea] - AR
        APer = Atotal - AR - AP
        NH = df.ResHouses[idarea] + df.HDRFlats[idarea]
        Blk_IWU = df.WD_In[idarea]*1000.0  # kL/d/block --> L/day/block
        if np.isnan(NH):
            NH = 0
        Aoccu = ((df.HouseOccup[idarea]*df.ResHouses[idarea]
                 +df.HDROccup[idarea]*df.HDRFlats[idarea])/max(NH,1.0))
        if np.isnan(Aoccu):
            Aoccu = 0
        #----------------------------------------------
        if GW.loc[idarea,'gw0mSL'] > 20:
            Qseep = 0
            GW0 = max(GW.loc[idarea,'gw0mSL'],0)
            h_dgw = GW.loc[idarea,'gwmmSL']
        else:
            Qseep = par.loc["down_seep",i] 
            GW0 = max(GW.loc[idarea,'gw0mSL'],0)
            h_dgw = max(GW.loc[idarea,'gwmmSL'],0)
        if (df.pLU_WAT[idarea]>0.0001): # [NEW:] & (GW.loc[idarea,'gw0mSL']<12.0)
            w = par.loc['w',i]
            #GW0 = 0.5
            #h_dgw = 0.8
        else:
            w = downOWds[OWi]**2/(8*w_k*w_D)
            #w = par.loc['w_deep',i]
        OWi = OWi+1
        #--------------------------------------------
        soil_prm = soil_selector(soilmatrix=soildata, etmatrix=etdata,a=par.loc["soiltype",i], b=par.loc["croptype",i])  #soiltype, croptype
        h0 = soil_prm[gwlcalc(GW.loc[idarea,'gw0mSL'])[2]]["moist_cont_eq_rz[mm]"]   #need correction
        if idarea in alternativew.index:
            altw = alternativew.loc[idarea,:]
        else:
            altw = alternativew.iloc[-1,:]
        if df.Blk_TIF[idarea] < 0.05*0.01:
            perI = 5
        else:
            perI = par.loc["perI",i]
        names['par%s' % idarea] = {
            "cellid":cellid,
            "dt": par.loc["dt",i],
            "soiltype": par.loc["soiltype",i],
            "croptype": par.loc["croptype",i],
            "AR": AR,
            "ERA": par.loc["ERA",i],
            "RIL": par.loc["RIL",i],
            "RST0": 0,
            "RTop": altw.RTop,
            "ART": altw.ART,
            "RTc": altw.RTc,
            "RTff": altw.RTff,
            "RT0": altw.RT0,
            "ERA_out": par.loc["ERA_out",i],
            "AP": AP,
            "EPA": par.loc["EPA",i],
            "PIL": par.loc["PIL",i],
            "PST0": 0,
            "Infilc_p": par.loc["Infilc_p",i],
            "APer": APer,
            "PerIL": par.loc["PerIL",i],
            "PerST0":0,
            "Infilc_per": par.loc["Infilc_per",i],
            "LR": par.loc["LR",i],
            "IRC": par.loc["IRC",i],
            "GW0": GW0,
            "seep_def": par.loc["seep_def",i],
            "w": w,
            "vc": par.loc["vc",i],
            "h_dgw": h_dgw,
            "down_seep": Qseep,
            "ASWS": altw.ASWS,
            "SWSop": altw.SWSop,
            "SWSc": altw.SWSc,
            "SWS0": altw.SWS0,
            "SWSff": altw.SWSff,
            "perI": perI,
            "AWWS": altw.AWWS,
            "WWSc": altw.WWSc,
            "WWS0": altw.WWS0,
            "pRT": altw.pRT,
            "Aoccu": Aoccu,
            "NH": NH,
            "AcWWS": altw.AcWWS,
            "cWWSc": altw.cWWSc,
            "cWWS0": altw.cWWS0,
            "h0": h0,
            "AUZ": APer,
            "AGW": Atotal,
            "IWU":Blk_IWU,
            "Elev":df.AvgElev[idarea],
            "OWele":520.0,
            "Blk_WDIrr": df.WD_Out[idarea]*365.0,
            "Population": df.Population[idarea]
        } #OW": par.iloc[39,i],
        paramlist.append(names['par%s' % idarea])
        if len(par.iloc[0,:])>2:
            i = i+1
            
    largepath = []
    for mid in df.index[0:]:
        nei = df.Neighbours[mid].split(',')
        if df.downID[mid] > 0.0:
            d = df.downID[mid]
        else:
            d = 0.0
        path0 = [mid, d]
        path = []
        # check if the downstream cell of the neighboring cells = main cell 
        # Boolean(df.downID[indint] == mid)
        for ind in nei:
            indint = int(ind)
            path.append(indint*(df.downID[indint]==mid))
        #sort(reverse):from large to small (put zeros to back)
        path.sort(reverse=True)
        path0.extend(path)  #[mid, d, up1, up2, ...]
        largepath.append(path0)

    Lpathdf = pd.DataFrame(largepath)
    Lpathdf = Lpathdf.fillna(0) #not every cell has all neighbors
    Lpathdf.set_index(0, inplace=True)
    if direction == 8:
        Lpathdf.columns = ["d", "u1", "u2", "u3", "u4","u5","u6","u7","u8"]
    elif direction == 6:
        Lpathdf.columns = ["d", "u1", "u2", "u3", "u4","u5","u6"]
    else:
        Lpathdf.columns = ["d", "u1", "u2", "u3", "u4"]
        
    return paramlist, Lpathdf


def read_forcing(PEpfile, OWfile, yr, deltat):
    """
    Raed daily precipitation data and yearly irrigation data.
    Transfer the irrigation data into daily data. (Yearly irrigation amount is distributed to each day based on\
    inverse of potential evaporation)
    Args:
        PEpfile(string):"Water balance data/forcing.csv"
        OWfile(string): "Water balance data/forcingOW.csv"
        yr(int): the year of data you would like to model for
        deltat(float): dict_param["dt"]
    return: 
        forcing2(dataframe): Daily data with first row = 0.
            daily data of P[mm], Ep[mm], OWt[m-SL], Ind, indnorm, IR_r, IR_p, IR_per[mm/d]
    """
    
    PEp = pd.read_csv(PEpfile, header = 0)
    PEp['Date'] = pd.to_datetime(PEp['Date'], format='%Y-%m-%d', errors='coerce')
    PEp.set_index('Date', inplace=True)
    OWmSL = pd.read_csv(OWfile, header = 0)
    OWmSL['Date'] = pd.to_datetime(OWmSL['Date'], format='%Y-%m-%d', errors='coerce')
    OWmSL.set_index('Date', inplace=True)
    forcing = pd.concat([PEp['P'],PEp['Ep'],OWmSL['OWt']], axis=1, join='inner')  #,OWElev-OWdfmean['OWt']
    forcing.columns=['P', 'Ep', 'OWt']
    
    # unit: mm/day --> mm
    forcing["P"] = forcing["P"] * deltat
    forcing["Ep"] = forcing["Ep"] * deltat
    forcing["Ind"] = forcing["Ep"]
    yearly = forcing["Ind"].groupby([forcing.index.year]).sum()
    indnorm = []
    for i in range(len(forcing)):
        indnorm.append(forcing["Ind"][i] / yearly[forcing.index.year[i]])

    forcing["indnorm"] = indnorm
    forcing["IR_r"] = 0.0              #IR_r
    forcing["IR_p"] = 0.0              #IR_p
    forcing["IR_per"] = 1.0            #IR_per
    
    #PEp = pd.read_csv(PEpfile, header = 0)
    forcing2 = pd.DataFrame(np.array([np.zeros(len(forcing.iloc[0,:]))]),
                            columns=forcing.columns).append(forcing, ignore_index=False)
    forcing2['Date'] = pd.to_datetime(forcing2.index, format='%Y-%m-%d', errors='coerce')
    forcing2.set_index('Date', inplace=True)
    #forcing2.index.name = 'Date'
    #forcingdata.insert(0,"Date",forcingdata.index)
    #forcing2018 = forcing2.drop(index = forcing2[(forcing2.index.year!=2018)].index)
    forcing1yr = forcing2.drop(index = forcing2.iloc[1:,:][(forcing2.iloc[1:,:].index.year!=yr)].index)
    return forcing1yr


def running(UrbanB, cali, size, direction, altwater, reusesetting, groundwater, PEpfile, OWfile, yr, demandfile, soilfile, etfile, detail=True):
    """
    
    Args:
        UrbanB(string): UrbanBEATS output file name (.dbf)
        cali(string):calibrated parameters file name, when SA, cali(dataframe)
        size(float): cell size (side length) [m]
        direction(int): Direction that is considered for neighborhoods. (options: 4, 6, 8)
        altwater(string): The dimension of alternative water device
                         (raintank, onsite WWS, SWS, cluster-scale WWS) data file name.
        reusesetting(string): The reusesetting file name ('wateruse/settingreuse.csv')
        groundwater (string): Initial groundwater level (in the begining of simulation) file name
        PEpfile(string): Climate data (Date/P/Ep) file name. (.csv)
        OWfile (string): Boundary condition open water level(Date/OW [m-SL]) file name. (.csv)
        yr(int):year of data that is going to be used
        demandfile(string): The percentage of the four water uses (K, B, T, L) [%]
        soilfile(string): Soilparameter data file name. ("in_soilparameter_new.csv")
        etfile(string): Etparameter data file name. ("in_etparameter.csv")
        detail(boolean): return df or not
        
        Return:
        df(list): list of results of each grid
        dict_param_all(list): list of parameters of each grid
        output(dataframe): result at model area outlet
        forcing (dataframe): summary of the forcing data (column: P/ Ep/ OWt/ Ind/ indnorm)
        #wbeach(dataframe): waterbalance check of each grid at each time step
        #WBasone (dataframe): water balance check at each time step if taking the whole model as one
    """
    #dict_param_all (list): item i-1(dict_param_all[i-1], type = dict) for model i (mi)
    dict_param_all, path = read_param2(UrbanB, cali, size, altwater, groundwater, soilfile, etfile, direction,1.5, 0.20)
    setreuseall = pd.read_csv(reusesetting, header = None, index_col = 0)
    set_len = len(setreuseall.iloc[0,:])
    demand_wu = pd.read_csv(demandfile, header = 0, index_col = 0)  #demand_of different water use
    soilmatrix = pd.read_csv(soilfile, header = 0)
    etmatrix = pd.read_csv(etfile, header =0)
    
    alter = pd.read_csv(altwater)
    cWWSind = alter.id[np.where(alter['cWWSc']!=0)[0]] #block id having cWWS
    SWSind = alter.id[np.where(alter['SWSc']!=0)[0]]   #block id having SWS
    
    forcing = read_forcing(PEpfile, OWfile, yr, dict_param_all[0]["dt"])
    Date = forcing.index
    P = forcing["P"]
    Ep = forcing["Ep"]
    Eref = forcing["Ep"]
    OWtg = forcing["OWt"]
    
    order = findordererUB(path,direction)  #calculation order from upstream to downstream
    number = len(path) #the number of grids
    iters = np.shape(P)[0]
    
    lst0 = [{"Rs":0.0, "Rw":0.0}]   # for zero upstream
    lst0.append({"Rs":0.0, "Rw":0.0})
    names = locals()
    for i in range(1, number+1):
        cellid = dict_param_all[i-1]["cellid"]
        #create IRi considering the surface area of mi
        names['IR%s' % cellid] = forcing[["IR_r","IR_p","IR_per"]]
#         if dict_param_all[i-1]["AR"] == 0:
#             names['IR%s' % cellid] = names['IR%s' % cellid].assign(IR_r = 0)
#         if dict_param_all[i-1]["AP"] == 0:
#             names['IR%s' % cellid] = names['IR%s' % cellid].assign(IR_p = 0)
        names['IR%s' % cellid] = names['IR%s' % cellid].assign(IR_r = 0)
        names['IR%s' % cellid] = names['IR%s' % cellid].assign(IR_p = 0)
        if dict_param_all[i-1]["APer"] == 0:
            names['IR%s' % cellid] = names['IR%s' % cellid].assign(IR_per = 0)
        else:
            #currently only feasible for whole year data, assume all block irrigation goes to Per
            projectedIR = dict_param_all[i-1]['Blk_WDIrr']*1000/dict_param_all[i-1]['APer']*forcing["indnorm"]
              #[kL/year]-->mm/day
            names['IR%s' % cellid] = names['IR%s' % cellid].assign(IR_per = projectedIR)
        #create mi for grid i, and the "answer sheet lsti"
        #demand_per = demand(demandfile, dict_param_all[i-1]['Aoccu'])  #demand per house
        demand_per = dict_param_all[i-1]["IWU"]/100.0 * demand_wu  # daily demand per block
        if set_len>1:
            s = i
        else:
            s = 1
        setreuse = setreuseall[s]
        names['m%s' % cellid] = NewModel(dict_param_all[i-1], soilmatrix, etmatrix, demand_per, setreuse)
        dict_param_all[i-1]['allART'] = names['m%s' % cellid].raintank.ART
        dict_param_all[i-1]['allRTc'] = names['m%s' % cellid].raintank.RTc
        dict_param_all[i-1]['allRTff'] = names['m%s' % cellid].raintank.RTff
        dict_param_all[i-1]['UZc'] = names['m%s' % cellid].pervious.UZc
        dict_param_all[i-1]['UZKsat'] = names['m%s' % cellid].pervious.UZKsat
        dict_param_all[i-1]['SSGS'] = names['m%s' % cellid].reuse.SSGS
        dict_param_all[i-1]['RTD0'] = names['m%s' % cellid].reuse.RTD0
        dict_param_all[i-1]['K0'] = names['m%s' % cellid].reuse.K0
        dict_param_all[i-1]['B0'] = names['m%s' % cellid].reuse.B0
        dict_param_all[i-1]['L0'] = names['m%s' % cellid].reuse.L0
        dict_param_all[i-1]['T0'] = names['m%s' % cellid].reuse.T0
        #dict_param_all[i-1]['IWU'] = names['m%s' % cellid].reuse.IWU
        dict_param_all[i-1]['allAWWS'] = names['m%s' % cellid].reuse.AWWS
        dict_param_all[i-1]['allWWSc'] = names['m%s' % cellid].reuse.WWSc
        
        names['lst%s' % cellid] =[
                {
                        "IR_r": np.nan,
                        "RSTt0": np.nan,
                        "Ea_r": np.nan,
                        "RSTt": dict_param_all[i-1]["RST0"],
                        "IRUN_r": np.nan,
                        "NEAR_r": np.nan,
                        "WB_r": np.nan,
                        "RTff_a": np.nan,
                        "In_rt": np.nan,
                        "RTt0": np.nan,
                        "Ea_rt": np.nan,
                        "RTt": dict_param_all[i-1]["RT0"],
                        "EXC_rt": np.nan,
                        "Out_r": np.nan,
                        "WB_rt": np.nan,
                        "IRUN_rrun": np.nan,
                        "IRUN_rp": np.nan,
                        "Inflow_rp": np.nan,
                        "IR_p": np.nan,
                        "PSTt0": np.nan,
                        "Ea_p": np.nan,
                        "PSTt": dict_param_all[i-1]["PST0"],
                        "Inf_p": np.nan,
                        "IRUN_p": np.nan,
                        "NEAR_p": np.nan,
                        "WB_p": np.nan,
                        "Inflow_per": np.nan,
                        "IR_per": np.nan,
                        "PerSTt0": np.nan,
                        "Infilc_pera": np.nan,
                        "timef_per": np.nan,
                        "Ea_per": np.nan,
                        "Inf_per": np.nan,
                        "PerSTt": dict_param_all[i-1]["PerST0"],
                        "EXC": np.nan,
                        "WB_per": np.nan,
                        "IRtot":np.nan,
                        "LD": np.nan,
                        "h3_uz": np.nan,
                        "alpha": np.nan,
                        "ET": np.nan,
                        "gw_up": np.nan,
                        "gw_low": np.nan,
                        "heq": np.nan,
                        "Cmax": np.nan,
                        "ht0": np.nan,
                        "UZ_per": np.nan,
                        "INFS": np.nan,
                        "ht": dict_param_all[i-1]["h0"],
                        "WB_uz": np.nan,
                        "Inflow_gw": np.nan,
                        "sc_gw": np.nan,
                        "GWt0": np.nan,
                        "Q_seep": 0.0,
                        "BF_out": 0.0,
                        "GWt": dict_param_all[i-1]["GW0"],
                        "GWabt": 0.0,
                        "WB_gw": np.nan,
                        "RUN_tot": np.nan,
                        "ISI": np.nan,
                        "RUN": np.nan,
                        "SWSff_a": np.nan,
                        "Rs_up": 0.0,
                        "In_sws": np.nan,
                        "SWSt0": np.nan,
                        "Ea_sws": np.nan,
                        "SWSt": dict_param_all[i-1]["SWS0"],
                        "EXC_sws": np.nan,
                        "Rs": 0.0,
                        "WB_sws": np.nan,
                        "In_cwws": np.nan,
                        "cWWSt": dict_param_all[i-1]["cWWS0"],
                        "WB_cwws": np.nan,
                        "Rw_up": 0.0,
                        "Rw": 0.0,
                    "SSGD": np.nan, "SSGuse": np.nan, "SSGsp":np.nan, "SSGdef": np.nan,"WWSD":np.nan,
                    "WWSS": np.nan, "WWSuse": np.nan, "WWSt": dict_param_all[i-1]["WWS0"],
                    "WWSdef": np.nan, "WWSsp":np.nan,"WB_wws":np.nan,"RTD":np.nan, "RTuse":np.nan,
                    "RTtend": np.nan, "RTdef": np.nan,"Import": np.nan,"KBL": np.nan,
                    "T": np.nan, "IR": np.nan, "checkssg":np.nan, "checkwws": np.nan, "checkrt":np.nan,
                    'cWWSuse':np.nan, 'SWSuse':np.nan, 'cWWSsup':np.nan, 'SWSsup':np.nan 
                    }
                ]
        
    output = [{"Rs_out":np.nan,"Rw_out":np.nan,"BFout_all":np.nan,
               "Qseep_all":np.nan, "Import_all":np.nan, "totalET": np.nan, "totalE": np.nan}]
    #flowpath
    df_param_all = pd.DataFrame(dict_param_all)
    df_param_all.set_index('cellid', inplace=True)
    
    print("Data loaded.")
    for t in range(1,iters):
        Qseep_all = 0.0
        BFout_all = 0.0
        Rs_out = 0.0  # possible multiple outlet
        Rw_out = 0.0  # possible multiple outlet
        Import_all = 0.0
        totalET = 0.0
        totalE = 0.0
        for i in order:  #each model in the same time step
            i = int(i)
            u = path.loc[i][1:]  #upstream grid no.
            d = path.d[i]
            Rsup = 0
            Rwup = 0
            for upup in u:
                Rsup = Rsup + names['lst%s' % int(upup)][t]["Rs"]
                Rwup = Rwup + names['lst%s' % int(upup)][t]["Rw"]
            lstup = {"Rs":Rsup, "Rw":Rwup}
            #print(Rsup, Rwup)
            names['lst%s' % i].append(names['m%s' % i].__next__(P.iloc[t], Ep.iloc[t], OWtg.iloc[t], names['IR%s' % i].iloc[t,:], 
                                                                Eref.iloc[t], names['lst%s' % i][t-1],lstup))
            if d == 0.0:
                Rs_out = Rs_out + names['lst%s' % i][t]["Rs"]
                Rw_out = Rw_out + names['lst%s' % i][t]["Rw"]
                #BFdown = names['lst%s' % i][t]["BF_down"]
            
            BFout_all = BFout_all + names['lst%s' % i][t]["BF_out"]  # + BFdown
            Qseep_all = Qseep_all + names['lst%s' % i][t]["Q_seep"]
            Import_all = Import_all + names['lst%s' % i][t]["Import"]
            E = ((names['lst%s' % i][t]["Ea_r"]*df_param_all.AR[i] + names['lst%s' % i][t]["Ea_p"]*df_param_all.AP[i] 
                  + names['lst%s' % i][t]["Ea_per"]*df_param_all.APer[i]) + names['lst%s' % i][t]["Ea_rt"] 
                 + names['lst%s' % i][t]["ET"] * df_param_all.AUZ[i] + names['lst%s' % i][t]["Ea_sws"])
            totalE = totalE + E
            totalET = totalET + names['lst%s' % i][t]["ET"] * df_param_all.AUZ[i]
        for w in cWWSind:
            forcWWS = [k for k in order]
            while names['lst%s' % w][t]["cWWSt"]>0:
                select = choice(forcWWS)
                setreuse = setreuseall[select]
                use1 = min(names['lst%s' % w][t]["cWWSt"],
                           names['lst%s' % select][t]["T"]*setreuse.cWWSforT)
                names['lst%s' % select][t]["T"] = names['lst%s' % select][t]["T"] - use1
                use2 = min(names['lst%s' % w][t]["cWWSt"]-use1,
                           names['lst%s' % select][t]["IR"]*setreuse.cWWSforIR)
                names['lst%s' % select][t]["IR"] = names['lst%s' % select][t]["IR"] - use2
                names['lst%s' % w][t]["cWWSt"] = names['lst%s' % w][t]["cWWSt"] - use1 - use2
                names['lst%s' % w][t]["cWWSuse"] = names['lst%s' % w][t]["cWWSuse"] + use1 + use2
                names['lst%s' % select][t]["cWWSsup"] = names['lst%s' % select][t]["cWWSsup"] + use1 + use2
                names['lst%s' % select][t]["Import"] = names['lst%s' % select][t]["Import"] - use1 - use2
                forcWWS.remove(select)
                if len(forcWWS) ==0:
                    break
            Import_all = Import_all - names['lst%s' % w][t]["cWWSuse"]
            
        for s in SWSind:
            forSWS = [k for k in range(1, number+1)]
            while names['lst%s' % s][t]["SWSt"]>0:
                select = choice(forSWS)
                setreuse = setreuseall[select]
                use3 = min(names['lst%s' % s][t]["SWSt"],
                           names['lst%s' % select][t]["T"]*setreuse.SWSforT)
                names['lst%s' % select][t]["T"] = names['lst%s' % select][t]["T"] - use3
                use4 = min(names['lst%s' % s][t]["SWSt"]-use3,
                           names['lst%s' % select][t]["IR"]*setreuse.SWSforIR)
                names['lst%s' % select][t]["IR"] = names['lst%s' % select][t]["IR"] - use4
                names['lst%s' % s][t]["SWSt"] = names['lst%s' % s][t]["SWSt"] - use3 - use4
                names['lst%s' % w][t]["SWSuse"] = names['lst%s' % w][t]["SWSuse"] + use3 + use4
                names['lst%s' % select][t]["SWSsup"] = names['lst%s' % select][t]["SWSsup"] + use3 + use4
                names['lst%s' % select][t]["Import"] = names['lst%s' % select][t]["Import"] - use3 - use4
                print(select)
                forSWS.remove(select)
                if len(forSWS) ==0:
                    break
            Import_all = Import_all - names['lst%s' % s][t]["SWSuse"]
          
        output.append({"Rs_out": Rs_out,
                       "Rw_out": Rw_out,
                       "BFout_all": BFout_all,
                       "Qseep_all": Qseep_all,
                       "Import_all": Import_all,
                       "totalET": totalET,
                       "totalE": totalE
                      })
        lst0.append({"Rs":0.0, "Rw":0.0})
    
    print("Simulation done. Processing results.")
    if detail:
        df = []
        for i in range(1, number+1):
            cellid = dict_param_all[i-1]["cellid"]
            names['df%s' % i] = pd.DataFrame(names['lst%s' % cellid])
            names['df%s' % i].insert(0, "Date", Date)
            names['df%s' % i]['Date'] = pd.to_datetime(names['df%s' % i]['Date'],
                                                       format='%Y-%m-%d', errors='coerce')
            names['df%s' % i].set_index('Date', inplace=True)
            df.append(names['df%s' % i])
        print("List of df built.")
    else:
        df=[]
    output = pd.DataFrame(output)
    output.insert(0, "Date", Date)
    output.set_index('Date', inplace=True)
    
    #wbeach = WaterBalanceCheck(result_all=df, param_all=dict_param_all) #, wbeach
    #WBasone = WaterBalanceCheckall(result_all=df, param_all=dict_param_all,output=output)
    
    return df, dict_param_all, output, forcing  #, wbeach, WBasone