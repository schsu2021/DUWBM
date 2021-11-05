# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:24:47 2020

@author: modified based on urbanwb package
"""

import pandas as pd
import DUWBM

#path = DUWBM.DUWBMdir / ".." / "input"

#soilmatrix = pd.read_csv(path / "in_soilparameter_new.csv", header = 0)
#etmatrix = pd.read_csv(path / "in_etparameter.csv", header =0)


def et_selector(etmatrix, a, b):
    """
    defines moisture content - related parameters based on given soil type(a) and crop type(b).
    etmatrix(dataframe)
    etmatrix = pd.read_csv("in_etparameter.csv", header =0)
    """
    # a --- soil type
    # b --- crop type
    sol = etmatrix.loc[(etmatrix.soil_type == int(a)) & (etmatrix.crop_type == int(b))]
    return sol


def soil_selector(soilmatrix, etmatrix, a, b):
    """
    returns a database of soil parameters namely equilibrium moisture content, maximum capillary rise,
    storage coefficient, saturated permeability and unsaturated permeability based on given soil type, crop type
    soilmatrix = pd.read_csv("in_soilparameter_new.csv", header = 0)
    etmatrix = pd.read_csv("in_etparameter.csv", header =0)
    """
    # a --- soil type
    # b --- crop type

    rootzone_thickness = 100 * et_selector(etmatrix, a, b)["th_rz_m"].values
    soil_prm = soilmatrix.loc[
        (soilmatrix.soil_type == int(a)) & (soilmatrix.th_rz == int(rootzone_thickness))
    ]
    # convert data frame to dictionary (list) for quick lookup.
    soil_prm = soil_prm.to_dict(orient="Records")

    return soil_prm
