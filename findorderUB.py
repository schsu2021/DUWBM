# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 09:42:01 2020

@author: ShuChen Hsu
"""

import pandas as pd
import numpy as np

def findordererUB(pathindex,direction):
    """
    Derive the calculation order
    Args:
    pathindex(dataframe): flow path data derived from read_param2 ["d", "u1", "u2", "u3", "u4",...]
    direction (int): number of neighbor (6,4,8)
    """
    #pathindex = pd.read_csv(path, header=None)
    #pathindex.index = list(range(1,len(pathindex)+1))
    #pathindex.columns = ["d", "u1", "u2", "u3", "u4"]
    d1 = np.where(pathindex.u1==0)  # u1=0: most upstream cell
    order = pathindex.index[d1]   #d1[0] +1    # calculation order
    order = np.insert(order, 0, 0)
    more = np.zeros(1)  # The cells that have more than one upstream
    #direction = len(pathindex.iloc[0,1:])

    for i in pathindex.iloc[d1]["d"]:     # for the downstream cells of the most upstream cells
        if i ==0:
            break
        if pathindex.loc[i,"u2"] ==0:     # if they only have one upstream cell
            order = np.append(order, i)   # it can be calculated
            down= pathindex.loc[i,"d"]    # down = the downstream of these cells
            if down !=0:
                if pathindex.loc[down,"u2"] !=0:   # if the u2 of these down cells != 0
                        more = np.append(more, down)   # the downstrean of these down cells will not be track
                while pathindex.loc[down,"u2"] ==0:    # if yes: keep track the ones with only one upstream 
                    order = np.append(order, down)
                    down = pathindex.loc[down,"d"]
                    if down == 0:
                        break
                    if pathindex.loc[down,"u2"] !=0:
                        more = np.append(more, down)   # until the downstream requires more than one inflow
        else:
            more = np.append(more, i)
    more = np.delete(more,0)
    more = np.unique(more)

    new = np.zeros(1)
    for mul in more:    # among the cells having more than one upstream, if their upstream is already available in the "order"
        if direction == 8:
            available = ((pathindex.loc[mul,"u1"] in order) and (pathindex.loc[mul,"u2"] in order)
                         and (pathindex.loc[mul,"u3"] in order) and (pathindex.loc[mul,"u4"] in order)
                        and (pathindex.loc[mul,"u5"] in order) and (pathindex.loc[mul,"u6"] in order)
                        and (pathindex.loc[mul,"u7"] in order) and (pathindex.loc[mul,"u8"] in order))
        elif direction == 6:
            available = ((pathindex.loc[mul,"u1"] in order) and (pathindex.loc[mul,"u2"] in order)
                     and (pathindex.loc[mul,"u3"] in order) and (pathindex.loc[mul,"u4"] in order)
                    and (pathindex.loc[mul,"u5"] in order) and (pathindex.loc[mul,"u6"] in order))
        else:
            available = ((pathindex.loc[mul,"u1"] in order) and (pathindex.loc[mul,"u2"] in order)
                     and (pathindex.loc[mul,"u3"] in order) and (pathindex.loc[mul,"u4"] in order))
        if available:                     # they can be calculated in the next order
            new = np.append(new,mul)      # we need to track the downstream of these new included cells
            order = np.append(order, mul) # include them in the "order"
            more = np.delete(more,np.where(more==mul))     # remove them from the waiting list
    new = np.delete(new, 0)
    #print(order, new, more)

    while len(order) < len(pathindex)+1:
        for i in pathindex.loc[new, "d"]:
            if i in order:
                continue
            if pathindex.loc[i,"u2"] ==0:
                order = np.append(order, i)
                down= pathindex.loc[i,"d"]
                if down == 0:
                        break
                if pathindex.loc[down,"u2"] !=0:
                        more = np.append(more, down)
                while pathindex.loc[down,"u2"] ==0:
                    order = np.append(order, down)
                    down = pathindex.loc[down,"d"]
                    if down == 0:
                        break
                    if pathindex.loc[down,"u2"] !=0:
                        more = np.append(more, down)
            else:
                more = np.append(more, i)
            #print(order, new, more)

        more = np.unique(more)
        new = np.zeros(1)
        for mul in more:
            if direction == 8:
                available = ((pathindex.loc[mul,"u1"] in order) and (pathindex.loc[mul,"u2"] in order)
                             and (pathindex.loc[mul,"u3"] in order) and (pathindex.loc[mul,"u4"] in order)
                             and (pathindex.loc[mul,"u5"] in order) and (pathindex.loc[mul,"u6"] in order)
                             and (pathindex.loc[mul,"u7"] in order) and (pathindex.loc[mul,"u8"] in order))
            elif direction == 6:
                available = ((pathindex.loc[mul,"u1"] in order) and (pathindex.loc[mul,"u2"] in order)
                             and (pathindex.loc[mul,"u3"] in order) and (pathindex.loc[mul,"u4"] in order)
                             and (pathindex.loc[mul,"u5"] in order) and (pathindex.loc[mul,"u6"] in order))
            else:
                available = ((pathindex.loc[mul,"u1"] in order) and (pathindex.loc[mul,"u2"] in order)
                         and (pathindex.loc[mul,"u3"] in order) and (pathindex.loc[mul,"u4"] in order))
            if available:
                order = np.append(order, mul)
                new = np.append(new,mul)
                more = np.delete(more,np.where(more==mul))
        new = np.delete(new, 0)
        #print(order, new, more)
    order = np.delete(order, 0)
    return order