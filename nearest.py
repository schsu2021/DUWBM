# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 10:20:20 2020

@author: ShuChen Hsu
"""

import numpy as np

def NearestDownOW(df, direction, size):
    """
    Derive the nearest downward block having open water.
    Args:
        df (dataframe): index = cellid, important column:downID, CentreX, CentreY
        direction (int): number of neighbor (6,4,8)
        size (float): side length of one cell
    Return:
        downOW (list): the index of the nearest downward block having open water [-]
        downOWd (list): the distance bewteen the block and its downOW [m]
        
    """
#     direction=8
#     size = 300.0
    downOW = []
    downOWd = []
    for idarea in df.index:
        upid = idarea
        dcc = [0, 0, 0]   #1, sqrt(2), sqrt(3)
        downid = df.downID[upid]
        if downid>min(df.index):
            #print(upid, downid, dcc) 
            while df.pLU_WAT[downid]<0.0001:
                if direction == 6:
                    dcc[2] = dcc[2] + 1
                else:
                    ds = (df.CentreX[upid]!=df.CentreX[downid]) & (df.CentreY[upid]!=df.CentreY[downid])  #Boolean:1 for sqrt(2), 0 for 1
                    dcc[0] = dcc[0] + 1*(1-ds)
                    dcc[1] = dcc[1] + 1*ds
                    #print(ds)
                upid = downid
                downid = df.downID[downid]
                if downid<min(df.index):
                    break
        if downid>min(df.index):
            if direction == 6:
                dcc[2] = dcc[2] + 1
            else:
                ds = (df.CentreX[upid]!=df.CentreX[downid]) & (df.CentreY[upid]!=df.CentreY[downid])  #Boolean:1 for sqrt(2), 0 for 1
                dcc[0] = dcc[0] + 1*(1-ds)
                dcc[1] = dcc[1] + 1*ds
            #print(upid, downid, dcc)

        downOW.append(downid)
        downOWd.append((dcc[0]+np.sqrt(2)*dcc[1]+np.sqrt(3)*dcc[2])*size)
        
    return downOW, downOWd