#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import sys
import math
import gc
import numpy as np

import pcraster as pcr

import virtualOS as vos


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
    
    # clone map
    clone_map = "/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/cloneMaps/australia_30sec.map"

    # output directory
    output_directory = "/scratch/depfg/sutan101/test_nicole/"
    cleanOutputDir = False
    if cleanOutputDir:
        try: 
            shutil.rmtree(output_directory)
        except: 
            pass # for new outputDir (not exist yet)
    try: 
        os.makedirs(output_directory)
    except: 
        pass # for new outputDir (not exist yet)
    
    # temporary directory
    tmp_directory = output_directory + "/tmp/"
    if os.path.exists(tmp_directory): shutil.rmtree(tmp_directory)
    os.makedirs(tmp_directory)
    

    # the REFERENCE map
    # ~ reference_raster  = pcr.readmap("/scratch/depfg/sutan101/data/gde_australia/example/aquatic_boolean.map")
    print("reading the reference map")
    reference_raster_file = "/scratch/depfg/otoo0001/comparison_matrix/shapefiles/aquatic/aquatic_boolean.map"
    reference_raster      = vos.readPCRmapClone(v = reference_raster_file, \
                                                cloneMapFileName = clone_map, \
                                                tmpDir = tmp_directory)
    reference_raster      = pcr.ifthen(reference_raster, reference_raster)

    
    # read the aquatic extent based on the extents of PCR-GLOBWB rivers, lakes and reservoirs
    print("set the aquatic extent")
    # - river width (m) file based on PCR-GLOBWB
    river_width_file = "/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/routing/channel_properties/version_2021-02-XX/maps_covered_with_zero/bankfull_width_channel_parameters_30sec_february_2021_global_covered_with_zero.nc"
    river_width = vos.readPCRmapClone(v = river_width_file, \
                                      cloneMapFileName = clone_map, \
                                      tmpDir = tmp_directory)
    # - threshold (m), below it, we assume as 
    river_width_threshold = 50. 
    # - extent of rivers
    river_extent = pcr.ifthen(river_width > river_width_threshold, pcr.boolean(1.0)) 
    # - pcr-globwb water body types
    water_body_type_file = "/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/routing/surface_water_bodies/version_2020-05-XX/maps/waterBodyTyp_lakes_and_reservoirs_30sec_global_2019_version_202005XX.map"
    water_body_type = vos.readPCRmapClone(v = water_body_type_file, \
                                          cloneMapFileName = clone_map, \
                                          tmpDir = tmp_directory)
    lake_and_reservoir_extent = pcr.ifthen(water_body_type > 0.0, pcr.boolean(1.0))
    # - aquatic extent                                     
    aquatic_extent = pcr.cover(river_extent, lake_and_reservoir_extent)
    # - make sure that the reference_raster is included in aquatic_exent
    aquatic_extent = pcr.cover(aquatic_extent, reference_raster)
    aquatic_extent = pcr.ifthen(aquatic_extent, aquatic_extent)
    
    # open groundwater model output maps
    # - variables that would be considered:
    groundwater_depth     = pcr.readmap("/scratch/depfg/otoo0001/pcrglobwb_gmglob_30sec_demo/australia_scenario_2/steady-state_only/states/groundwaterDepthLayer2_1958-01-01.ini.masked.map")
    groundwater_discharge = pcr.readmap("/scratch/depfg/otoo0001/pcrglobwb_gmglob_30sec_demo/australia_scenario_2/steady-state_only/maps/baseflow_1958-01-01.ini.map")
    
    # pixels that are classified as groundwater dependent ecosystems based on the model
    gde_based_on_model_depth        = pcr.ifthen(groundwater_depth < 10.0, pcr.boolean(1.0))
    gde_based_on_model_gw_discharge = pcr.ifthen(groundwater_discharge > 0.0005, pcr.boolean(1.0))
    gde_based_on_model              = pcr.cover(gde_based_on_model_depth, gde_based_on_model_gw_discharge)
    
    # focus our analysis on the aquatic_extent only
    reference_raster   = pcr.ifthen(aquatic_extent, reference_raster)
    gde_based_on_model = pcr.ifthen(aquatic_extent, gde_based_on_model)

    # convert the reference and model classification array to nummpy
    np_reference_raster   = pcr.pcr2numpy(reference_raster, 0)
    np_gde_based_on_model = pcr.pcr2numpy(gde_based_on_model, 0)
    
    hit_rate_score = hit_rate(np_reference_raster, np_gde_based_on_model)
    print("hit rate")
    print(hit_rate_score)   

    false_alarm_rate_score = false_alarm_rate(np_reference_raster, np_gde_based_on_model)
    print("false_alarm_rate")
    print(false_alarm_rate_score)   

    critical_success_index = critical_success(np_reference_raster, np_gde_based_on_model)
    print("critical_success_index")
    print(critical_success_index)   
        
if __name__ == '__main__':
    sys.exit(main())
