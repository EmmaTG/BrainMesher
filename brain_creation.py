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
import os


pathIn = os.getcwd()
# pathIn = "./"
fileIn = '/mri/aseg_tumor.mgz'

pathOut = "/".join([pathIn, "Models"])
if not os.path.exists(pathOut):
    os.mkdir(pathOut)

fileOut = "aseg_tumor"

brainCreator = BrainHexMesh()

## Preferences are defined in ConfigFile
config = ConfigFile(pathIn, fileIn, pathOut, fileOut)
brainCreator.setConfig(config)

# Writes configuarion preferences to output location
config.openConfigFile()

# Gets aseg data as 3D of segmentation labels in voxels
data = brainCreator.import_data()

# Homogenizes data label according to Materials label as defined in Config class
data = brainCreator.homogenize_data(data)

# Pre-processes data to ensure valid mesh:
preprocessor = Preprocessor.PreProcessorFactory.get_preprocessor(config, data)
assert isinstance(preprocessor, Preprocessor.IPreprocessor)
data_new = preprocessor.preprocess_data()

# Creates point cloud from voxel data
pointCloud = PointCloud()
pointCloud.create_point_cloud_from_voxel(data_new)

# Creates mesh from pointcloud
mesh = brainCreator.make_mesh(pointCloud.pcd) 
brainCreator.clean_mesh(mesh)

# Moves mesh to the center of the corpus callosum
mesh.center_mesh(251)

postProcessor = PostProcessor(config, mesh)
if config.atrophy:
    postProcessor = ApplyAtrophyConcentration(postProcessor)
if config.Smooth:
    postProcessor = SmoothMesh(postProcessor, [0.6, -0.4], 2)
if config.Refine:
    postProcessor = RefineMesh(postProcessor, 'elements', elements=list(mesh.elements.keys())[:10])

all_labels = config.material_labels.get_homogenized_labels_map()
label_for_unused = all_labels.get("Unused")
postProcessor = RemoveRegion(postProcessor, label_for_unused)
if not config.ventricles:
    all_labels = config.material_labels.get_homogenized_labels_map()
    label_for_ventricles = all_labels.get("Ventricles",-1)
    if label_for_ventricles != -1:
        postProcessor = RemoveRegion(postProcessor, label_for_ventricles)


postProcessor.post_process()


# if config.atrophy:
#     brainCreator.apply_atrophy_concentration(mesh)

# Write mesh to file (ucd, vtk or abq inp as specified in config file)
brainCreator.write_to_file(mesh)  

# Close config file write out  
config.closeConfigFile()


