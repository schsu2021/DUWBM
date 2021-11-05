# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 00:26:57 2020

@author: ShuChen Hsu
"""

import pandas as pd
from DUWBM.exploreQQ import exploreQ


def ExportForGIS(flux, trans, itsA, df, param_all, period, saveat, filename, yrsum=False):
    """
    Args:
        flux(str): name of the flux in model result (df) you are interested in
        trans(boolean): want to transfer from mm to L ot not (yes:1, no:0)
        itsA(string): the corresponding area variable name for the transfer
        df(list of dataframe)
        param_all(list of dictionary)
        period1(list of date in str): The dates you would like to export, needed if yrsum=False
        saveat(str): the folder name
        filename(str)
        yrsum(Boolean): true for yearly sum (transA and period input doesn't matter)
    Return:
        (Export filename.csv to the folder "saveat")
        cellID_QTsT(dataframe): row:time, column: cells, useful when you would like to plot one cell
    """
    
    # Change the datatype of the datetime index to string so QGIS can read properly
    cellID_QTsT = pd.DataFrame([])
    if yrsum:
        for i in range(len(flux)):
            fluxi = flux[i]
            itsAi = itsA[i]
            if type(itsAi)==str:
                transi = 1
            else:
                transi = 0
            # input trans doesn't matter
            Qi= exploreQ(fluxi, transi, itsAi, df, param_all)
            yrQi = Qi.groupby([Qi.index.year]).sum()
            cellID_QTsT[fluxi] = yrQi.T.iloc[0:,0]
        cellID_QTsT['cellid'] = yrQi.T.index
        exportpath = saveat + filename + '.csv'
        cellID_QTsT.to_csv(exportpath)
    else:
        cellID_QTsT = exploreQ(flux, trans, itsA, df, param_all)
        Date = ['0']
        for day in range(1,len(cellID_QTsT)):
            Date.append(str(cellID_QTsT.index[day].date()))
        cellID_QTsT['strDate'] = Date
        cellID_QTsT.set_index('strDate', inplace=True)
        cellID_QTs = cellID_QTsT.T
        cellID_QTsE = cellID_QTs[period]
        exportpath = saveat + filename + '.csv'
        cellID_QTsE.to_csv(exportpath)
    return cellID_QTsT