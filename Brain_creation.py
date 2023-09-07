# -*- coding: utf-8 -*-
"""
Brain creator 

Created on Thu Jun 15 09:34:00 2023

@author: grife

This script runs the full execution required to create a 3D brain model from a
freesurfer aseg file.

Input file given by path + fileIn
Preferences for model creation are defined in ConigFile class
Output is writen to same path as input

main - the main function of the script
"""

from BrainHexMesh import BrainHexMesh
from PointCloud import PointCloud
from Config import ConfigFile
import VoxelDataUtils as bm 

# def main():
path = "C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\PythonScripts\\BrainMesher"
fileIn = '\\mri\\aseg_tumor.mgz'

## Preferences are defined in ConfigFile
config = ConfigFile(path,fileIn)

brainCreator = BrainHexMesh()
brainCreator.setConfig(config)
# Writes configuarion preferences to output location
config.openConfigFile()

# Gets aseg data as 3D of segmentation labels in voxels
data = brainCreator.import_data() 
data = brainCreator.homogenize_data(data, unusedLabel="Unused") 

# Pre-processes data to ensure valid mesh, options: lesion edemicTissue ventricles
data = brainCreator.preprocessFactory(data, "basic")

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

# Laplacian smoothing         
if config.Smooth:
    brainCreator.smooth_mesh(mesh)       

# Write mesh to file (ucd, vtk or abq inp as specified in config file)
brainCreator.write_to_file(mesh)  

# Close config file write out  
config.closeConfigFile()


