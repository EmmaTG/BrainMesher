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
import os


# pathIn = os.getcwd()
pathIn = "./"
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
pointCloud.create_point_cloud_from_voxel(data) 

# Creates mesh from pointcloud
mesh = brainCreator.make_mesh(pointCloud.pcd) 
brainCreator.clean_mesh(mesh)

# Moves mesh to the center of the corpus callosum
mesh.center_mesh(251)

# Removes elements associated with a region to be excluded as defined in config file
all_labels = config.material_labels.get_homogenized_labels_map()
label_for_unused = all_labels.get("Unused")
mesh.remove_region(label_for_unused)

### Optional local refinement
meshRefiner = Refiner.Refiner(mesh)
bounds = [245, 255, 187, 195, 153, 165]
meshRefiner.refine_within_region(bounds)
config.writeToConfig("Refined within bounds", ",".join([str(x) for x in bounds]))
element_to_refine = list(mesh.elements.keys())[:10]
meshRefiner.refine_elements(element_to_refine)
config.writeToConfig("Refined elements", ",".join([str(x) for x in element_to_refine]))
point = [272, 190, 192]
radius = 5
meshRefiner.refine_around_point(point, radius)
config.writeToConfig("Refined around point", ",".join([str(x) for x in point]))
config.writeToConfig("Refined with radius", str(radius))

# Laplacian smoothing         
if config.Smooth:
    brainCreator.smooth_mesh(mesh)

if config.atrophy:
    brainCreator.apply_atrophy_concentration(mesh)

# Write mesh to file (ucd, vtk or abq inp as specified in config file)
brainCreator.write_to_file(mesh)  

# Close config file write out  
config.closeConfigFile()


