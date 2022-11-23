#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import sys
import math
import gc
import numpy as np

import pcraster as pcr

def hit_rate(array1, array2):
    """
    calculate the hit rate based upon 2 boolean maps. (i.e. where are both 1)
    """
    # count the number of cells that are flooded in both array1 and 2
    idx_both = np.sum(np.logical_and(array1, array2))
    idx_1 = np.sum(array1)
    return float(idx_both)/float(idx_1)

def false_alarm_rate(array1, array2):
    """
    calculate the false alarm rate based upon 2 boolean maps. (i.e. amount of cells where array2 is True but array1 False)
    """
    # count the number of cells that are flooded in both array1 and 2
    idx_2_only = np.sum(np.logical_and(array2, array1!=1))
    idx_2_total = np.sum(array2)
   
    return float(idx_2_only)/float(idx_2_total)

def critical_success(array1, array2):
    """
    calculate the critical success rate based upon 2 boolean maps.
    """
    idx_both = np.sum(np.logical_and(array1, array2))
    idx_either = np.sum(np.logical_or(array1, array2))
    return float(idx_both)/float(idx_either)

# source for hit_rate: https://github.com/edwinkost/glofris_inundation_benchmarking


def main():

    # note all maps should be in the WGS84 with the resolution 30 arcsec
    
    # open your reference:
    reference_raster = pcr.readmap("/scratch/depfg/sutan101/data/gde_australia/example/aquatic_boolean.map")
    
    # open groundwater model output maps
    # - variables that would be considered:
    groundwater_depth     = pcr.readmap("/scratch/depfg/otoo0001/pcrglobwb_gmglob_30sec_demo/australia_scenario_2/steady-state_only/states/groundwaterDepthLayer2_1958-01-01.ini.masked.map")
    groundwater_discharge = pcr.readmap("/scratch/depfg/otoo0001/pcrglobwb_gmglob_30sec_demo/australia_scenario_2/steady-state_only/maps/baseflow_1958-01-01.ini.map")
    
    # pixels that are classified as groundwater dependent ecosystems based on the model
    gde_based_on_model = pcr.ifthen(groundwater_depth < 10.0, pcr.boolean(1.0))
    # ~ gde_based_on_model = pcr.cover(gde_based_on_model, pcr.ifthen(groundwater_discharge lt XX.X, pcr.boolean(1.0)))
    
    # convert the reference and model classification array to nummpy
    np_reference_raster   = pcr.pcr2numpy(reference_raster, 0)
    np_gde_based_on_model = pcr.pcr2numpy(gde_based_on_model, 0)
    
    hit_rate_score = hit_rate(np_reference_raster, np_gde_based_on_model)
    print("test")
    print(hit_rate_score)   
