# -*- coding: utf-8 -*-
"""
Created on Sat May  2 10:40:55 2020

@author: urbanwb package
"""

import math


def gwlcalc(x):
    """
    calculates the groundwater up and groundwater low and the corresponding indexes for referencing in the database.

    Args:
        x (float): groundwater level

    Returns:
        (float): First value in predefined table above (below) groundwater level and its corresponding index:

        * **gwl_up** -- First value in predefined table above groundwater level at the end of previous time step [m-SL]
        * **gwl_low** -- First value in predefined table below groundwater level at the end of previous time step [m-SL]
        * **index** -- index of gwl_up value in the database.
        * **index2** -- index of gwl_low value in the database.
    """

    gwl_up = float(x)
    if 0.0 <= gwl_up <= 2.5:
        gwl_up = math.floor(gwl_up * 10.0) / 10.0
        index = int(gwl_up * 10.0)
    elif gwl_up < 3.0:
        gwl_up = 2.5
        index = 25
    elif gwl_up < 5.0:
        gwl_up = int(gwl_up)
        index = 23 + gwl_up
    elif gwl_up < 10.0:
        gwl_up = 5.0
        index = 28
    else:
        gwl_up = 10.0
        index = 29

    if gwl_up < 2.5:
        gwl_low = round(gwl_up + 0.1, 2)
        index2 = index + 1
    elif gwl_up < 3.0:
        gwl_low = 3.0
        index2 = index + 1
    elif gwl_up < 4.0:
        gwl_low = 4.0
        index2 = index + 1
    elif gwl_up < 5.0:
        gwl_low = 5.0
        index2 = index + 1
    else:
        gwl_low = 10.0
        index2 = 29

    return gwl_up, gwl_low, index, index2