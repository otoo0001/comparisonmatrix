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

# clone map
    clone_map = "/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/cloneMaps/australia_30sec.map"

# output directory
    output_directory = "/scratch/depfg/otoo0001/test_nicole_phreatophytes_hitrate/new_hitrates/"
#   output_directory = "/scratch/depfg/sutan101/test_nicole_with_plot/"
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
        
    try: 
        os.makedirs(output_directory)
    except: 
        pass # for new outputDir (not exist yet)    
        
    
    
    # output pcraster map file name
    pcraster_out_filename = "contigency_phreatophytes_fan.map"
    pcraster_out_filename = output_directory + "/" + pcraster_out_filename
    
    
    # temporary directory
    tmp_directory = output_directory + "/tmp/"
    if os.path.exists(tmp_directory): shutil.rmtree(tmp_directory)
    os.makedirs(tmp_directory)

    

    # note all maps should be in the WGS84 with the resolution 30 arcsec
    
    # open your reference:
    reference_raster_file ="/scratch/depfg/otoo0001/comparison_matrix/shapefiles/terrestrial/terrestrial_boolean.map"
    reference_raster      = vos.readPCRmapClone(v = reference_raster_file, \
                                                cloneMapFileName = clone_map, \
                                                tmpDir = tmp_directory)
   
    reference_raster      = pcr.ifthen(reference_raster, reference_raster) 
    
    
    

    # read the aquatic extent based on the extents of PCR-GLOBWB rivers, lakes and reservoirs
    
    # - river width (m) file based on PCR-GLOBWB
    river_width_file = "/scratch/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/routing/channel_properties/version_2021-02-XX/maps_covered_with_zero/bankfull_width_channel_parameters_30sec_february_2021_global_covered_with_zero.nc"
    river_width = vos.readPCRmapClone(v = river_width_file, \
                                      cloneMapFileName = clone_map, \
                                      tmpDir = tmp_directory)
    # - threshold (m), below it, we assume as 
    river_width_threshold = 30. 
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
    aquatic_extent = pcr.cover(aquatic_extent, pcr.boolean(0.0))

    
  
    # Create non-aquatic extent
    print("set the non-aquatic extent")
    non_aquatic_extent = pcr.ifthenelse(aquatic_extent, pcr.boolean(0.0), pcr.boolean(1.0))


      # read the natural vegetation extent based on the bare cover and crop cover from copernicus https://land.copernicus.eu/global/
    print ("set the natural vegetation extent")
    
    # - natural vegetation file based on copernicus land cover data 
    bare_cover_fraction_file = pcr.readmap ("/scratch/depfg/otoo0001/comparison_matrix/shapefiles/landcover/barecover/Bare_cover_30sec_Australia.map")
    
#   bare_cover_fraction = vos.readPCRmapClone(v = bare_cover_fraction_file, \
#                                             cloneMapFileName = clone_map, \
#                                             tmpDir = tmp_directory)
#    
    
    
    
    # extent of vegetated areas 
    non_bare_areas = pcr.ifthen(bare_cover_fraction_file < 50 , pcr.boolean(1.0)) 
    
    
    
    # extent of non-cropped areas 
    crop_area_fraction_file = pcr.readmap ("/scratch/depfg/otoo0001/comparison_matrix/shapefiles/landcover/cropcover/Crop_cover_frac_30secs_Australia.map")
    
    non_cropped_areas = pcr.ifthen(crop_area_fraction_file < 50 , pcr.boolean(1.0))
    

    # - natural vegetaion extent                                     
    
    # - make sure that the reference_raster is included in natural vegetation exent
    
    natural_vegetation_extent = pcr.cover(non_bare_areas,non_cropped_areas )
    
    natural_vegetation_extent = pcr.ifthen(natural_vegetation_extent, reference_raster)
    natural_vegetation_extent = pcr.ifthen(non_aquatic_extent, natural_vegetation_extent)
#    natural_vegetation_extent = vegetated_areas

   


    
    # open groundwater model output maps
    # - variables that would be considered:
    groundwater_depth     = pcr.readmap("/scratch/depfg/otoo0001/pcrglobwb_gmglob_30sec_demo/runs_144/35/steady-state_only/states/groundwaterDepthLayer1_1958-01-01.ini.masked.map")
    
    rooting_depth          = pcr.readmap("/scratch/depfg/otoo0001/comparison_matrix/shapefiles/landcover/root_depth/maxroot_depth_Australi30seca.map")
       
#    rooting_depth_map = "/scratch/depfg/otoo0001/comparison_matrix/shapefiles/landcover/root_depth/maxroot_Australia.nc"
#    rooting_depth = vos.readPCRmapClone(v = rooting_depth_map, \
#                                      cloneMapFileName = clone_map, \
#                                      tmpDir = tmp_directory) 
#       

#    rooting_depth_map = "/depfg/sutan101/data/pcrglobwb_input_arise/develop/global_30sec/landSurface/landCover/naturalVegetationAndRainFedCrops_version_2021-02-XX/maxrootdepth_all.nc"
#    rooting_depth = vos.readPCRmapClone(v = rooting_depth_map, \
#                                      cloneMapFileName = clone_map, \
#                                      tmpDir = tmp_directory)
#                                      
    pcr.aguila(rooting_depth)
     
                                      
    # pixels that are classified as groundwater dependent ecosystems based on the model
    gde_based_on_model      = pcr.ifthen(groundwater_depth <  rooting_depth, pcr.boolean(1.0))
    
    
    
#   gde_based_on_model_gw_discharge = pcr.ifthen(groundwater_discharge > 0.00, pcr.boolean(1.0))




    
    
    
  
    
    # limit the evaluation extent to natural vegetation only
    reference_raster   = pcr.ifthen(natural_vegetation_extent, reference_raster)
    gde_based_on_model = pcr.ifthen(natural_vegetation_extent, gde_based_on_model) 
 
 
    
    
    
    
    # convert the reference and model classification array to nummpy
    np_reference_raster   = pcr.pcr2numpy(reference_raster, 0)
    np_gde_based_on_model = pcr.pcr2numpy(gde_based_on_model, 0)
    
    hit_rate_score = hit_rate(np_reference_raster, np_gde_based_on_model)
    print("test_Hitrate_terrestrial")
    print(hit_rate_score)   
    
    
    false_alarm_rate_score = false_alarm_rate(np_reference_raster, np_gde_based_on_model)
    print("test_Falsealarm_terrestrial")
    print(false_alarm_rate_score)
    
    critical_success_index = critical_success(np_reference_raster, np_gde_based_on_model)
    print("critical_success_index")
    print(critical_success_index)
    
    # plot evaluation map
    print("plot evaluation map")
    model_only = pcr.ifthen(gde_based_on_model , pcr.ifthenelse(pcr.cover(reference_raster, pcr.boolean(0.0)), pcr.boolean(0.0), pcr.boolean(1.0)))
    model_only = pcr.ifthen(model_only, model_only)
    refer_only = pcr.ifthen(reference_raster, pcr.ifthenelse(pcr.cover(gde_based_on_model , pcr.boolean(0.0)) , pcr.boolean(0.0), pcr.boolean(1.0)))
    refer_only = pcr.ifthen(refer_only, refer_only)
    both_agree = pcr.ifthen(gde_based_on_model , pcr.ifthenelse(pcr.cover(reference_raster, pcr.boolean(0.0)), pcr.boolean(1.0), pcr.boolean(0.0)))
    both_agree = pcr.ifthen(both_agree, both_agree)
    
    # plot in values
    refer_only_class = pcr.ifthen(refer_only, pcr.scalar(-1.0))
    both_agree_class = pcr.ifthen(both_agree, pcr.scalar( 0.0))
    model_only_class = pcr.ifthen(model_only, pcr.scalar( 1.0))
    all_classes      = pcr.cover(refer_only_class, both_agree_class)
    all_classes      = pcr.cover(all_classes     , model_only_class)
    
    pcr.aguila(all_classes)
    
    
    # plot and save the map in a file
#    pcr.report(all_classes, pcraster_out_filename)
 
        
if __name__ == '__main__':
    sys.exit(main())

   