# -*- coding: utf-8 -*-
"""
Brain creator 

Created on Thu Jun 15 09:34:00 2023

@author: grife

This script runs the full execution required to create a 3D brain model from a
freesurfer aseg file.

Input file given by path + fileIn
Preferences for model creation are defined in ConfigFile class
Output is a folder Models created in current working directory

main - the main function of the script
"""

from BrainHexMesh import BrainHexMesh
from point_cloud.PointCloud import PointCloud
from config.Config import ConfigFile
from voxel_data import Preprocessor
from mesh.refinement import Refiner
from mesh.PostProcessor import *
import os, sys

# configFilePath = sys.argv[1]
configFilePath = sys.argv[1]

## Preferences are defined in ConfigFile
config = ConfigFile(configFilePath)

brainCreator = BrainHexMesh(config)

# Writes configuration preferences to output location
config.open_config_file()

# Gets aseg data as 3D of segmentation labels in voxels
data = brainCreator.import_data()

# Homogenizes data label according to Materials label as defined in Config class
data = brainCreator.homogenize_data(data)

# Pre-processes data to ensure valid mesh base on config setting:
preprocessor = Preprocessor.PreProcessorFactory.get_preprocessor(config, data)
assert isinstance(preprocessor, Preprocessor.IPreprocessor)
data_new = preprocessor.preprocess_data()

# Creates point cloud from voxel data
pointCloud = PointCloud()
pointCloud.create_point_cloud_from_voxel(data_new)

# Creates mesh from point cloud
mesh = brainCreator.make_mesh(pointCloud.pcd) 
brainCreator.clean_mesh(mesh)

# Moves mesh to the center of the corpus callosum
mesh.center_mesh(251)

# Wrapping of post-processing operations (operation selection defined in config file)
postProcessor = PostProcessor(config, mesh)
if config.get('atrophy'):
    # Creates a concentration profile in the brain to eb used for degree fo atrophy calculations
    postProcessor = ApplyAtrophyConcentration(postProcessor)

if config.get('smooth'):
    #global smoothing
    postProcessor = SmoothMesh(postProcessor, config.get('co_effs'), config.get('iterations'))

    if config.get('add_csf'):
        labels = config.material_labels.get_homogenized_labels_map()
        label_for_csf = labels.get("CSF", -1)
        postProcessor = SmoothMesh(postProcessor,
                                   config.get('co_effs'), config.get('iterations'), excluded_regions=[label_for_csf])

    # Smooths mesh as defined in config file
    regions = config.get("smooth_regions")
    region_coeffs = config.get("region_co_effs")
    region_iterations = config.get("region_iterations")
    count = 0
    for r in regions:
        labels = config.material_labels.get_homogenized_labels_map()
        label_for_region = labels.get(r, -1)
        if label_for_region != -1:
            excluded_materials = []
            labels.pop(r)
            for e in labels.values():
                if not excluded_materials.count(e):
                    excluded_materials.append(e)
            postProcessor = SmoothMesh(postProcessor, region_coeffs[count], region_iterations[count], excluded_regions=excluded_materials)
        count += 1


# Removes regions specified at unused in the materials label description
all_labels = config.material_labels.get_homogenized_labels_map()
label_for_unused = all_labels.get("Unused")
if label_for_unused != 0:
    postProcessor = RemoveRegion(postProcessor, label_for_unused)

# Creates blank spaces in areas of the ventricles if specified
all_labels = config.material_labels.get_homogenized_labels_map()
label_for_ventricles = all_labels.get("Ventricles", all_labels.get("Ventricle", -1))
if label_for_ventricles != -1:
    postProcessor = RemoveRegion(postProcessor, label_for_ventricles)

# Creates boundary elements on specified regions
for count, e_num in enumerate(config.get('boundary_element_numbers')):
    postProcessor = CreateBoundaryElements(postProcessor, e_num,
                                           boundary_test_fx=config.get('boundary_tests')[count],
                                           excluded_regions=config.get('excluded_regions')[count])

if config.get('refine'):
    # Refines mesh in ways defined in config file
    if config.get('refine.point'):
        postProcessor = RefineMesh(postProcessor, 'point')
    if config.get('refine.bounding_box'):
        postProcessor = RefineMesh(postProcessor, 'bounding_box',)
    if config.get('refine.elements'):
        postProcessor = RefineMesh(postProcessor, 'elements')
# Performs are post-processing steps
postProcessor.post_process()

# Write mesh to file (ucd, vtk or abq inp as specified in config file)
brainCreator.write_to_file(mesh)  

# Close config file write out  
config.close_config_file()


