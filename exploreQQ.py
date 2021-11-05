# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:28:44 2020

@author: ShuChen Hsu
"""


import pandas as pd

def exploreQ(flux, trans, itsA, df, param_all):
    """
    Args:
        df(list of dataframe)
        flux(string): name of the flux in model result (df) you are interested in
        trans(boolean): want to transfer from mm to L ot not (yes:1, no:0)
        itsA(string): the corresponding area variable name for the transfer
        df(list of dataframe)
        param_all(list of dictionary)
    Return:
        exportQ(dataframe): row:time, column: cells
    """
    column = [param_all[0]['cellid'], param_all[1]['cellid']]
    #column = [0,1]
    if trans:
        Q0 = pd.DataFrame(df[0][flux]) * param_all[0][itsA]
        Q1 = pd.DataFrame(df[1][flux]) * param_all[1][itsA]
    else:
        Q0 = pd.DataFrame(df[0][flux])
        Q1 = pd.DataFrame(df[1][flux])
    exportQ = Q0.join(Q1, lsuffix='0', rsuffix='1')
    for i in range(2,len(df)):
        if trans:
            mmtoL = param_all[i][itsA]
        else:
            mmtoL = 1.0
        Q = pd.DataFrame(df[i][flux]) * mmtoL
        exportQ = exportQ.join(Q, rsuffix=str(i))
        column.append(param_all[i]['cellid'])
        #column.append(i)
    
    exportQ.columns = column
    return exportQ