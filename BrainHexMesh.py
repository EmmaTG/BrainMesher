"""
Created on Wed May 10 09:11:56 2023

@author: grife
"""

import numpy as np
# import mne
# import pooch
import nibabel
from BrainModel import BrainModel
from Material_Label import Material_Label
from Maze_Solver import Maze_Solver
from PointCloud import PointCloud
from Mesh import Mesh

# Step 1: Using freesurfer and 'recon-all' create mri outputs. Ensure aseg.mgz is created.
t1_file = 'aseg.mgz'
t1 = nibabel.load(t1_file)
# t1.orthoview()
data = np.asarray(t1.dataobj)

# data_test = data[80:120, 80:120, 80:120]
# data = data_test
# data = np.ones((6,6,6))*2
# data[2,0,2] = 0
# data[2,1,2] = 0
# data[2,2,2] = 0

# Step 2: Determine segmentation of brain model via labels map
material_labels  = Material_Label()
material_labels.addLabelToMap('BrainStem', 16)
material_labels.addLabelToMap('GreyMatter', [3,42]) # Left, Right
material_labels.addLabelToMap('WhiteMatter' , [2,41,77]); # Left, Right, WM-hypointensities
material_labels.addLabelToMap('Corpuscallosum' , [251,252,253,254,255]); # CC_Posterior, CC_Mid_Posterior, CC_Central, CC_Mid_Anterior, CC_Anterior
material_labels.addLabelToMap('BasalGanglia' , [11,50,12,51,13,52,26,58,62,30]); # Caudate(L&R), Putamen(L&R), Palladium(L&R), Accumbens Area(L&R), vessel(L&R)
material_labels.addLabelToMap('Cerebellum' , [7,46,8,47]); # WM(L&R), GM(L&R)
material_labels.addLabelToMap('Thalamus' , [10,49,28,60]); # Thalamus(L&R), Ventral DC(L&R)
material_labels.addLabelToMap('Hippocampus' , [17,53]); # Left, Right
material_labels.addLabelToMap('Amygdala' , [18,54]); # Left, Right
material_labels.addLabelToMap('Lesion' , [25,57]); # Left, Right
material_labels.addLabelToMap('CSF' , [24]); # Left, Right

# Unused labels (will be set to 0)
# material_labels.addLabelToMap('Left-Lateral-Ventricle' , [4]);
# material_labels.addLabelToMap('Right-Lateral-Ventricle' , [43]);
# material_labels.addLabelToMap('Left-Inf-Lat-Vent' , [5]);
# material_labels.addLabelToMap('Right-Inf-Lat-Vent ' , [44]);
# material_labels.addLabelToMap('3rd-Ventricle' , [14]);
# material_labels.addLabelToMap('4th-Ventricle' , [15]);
# material_labels.addLabelToMap('Left-choroid-plexus' , [31]);
# material_labels.addLabelToMap('Right-choroid-plexus' , [63]);
# material_labels.addLabelToMap('Optic-Chiasm' , [85]);

# Homogenize labels
data = material_labels.homogenize_material_labels(data);

brainModel = BrainModel()
# Coarsen the brain model
print("########## Coarsening data ##########")
voxel_size = 2
data = brainModel.coarsen(voxel_size, data)

# Clean image removing isolated pixels and small holes
print("########## Performing cleaning operations on the data ##########")
brainModel.clean_mesh_data(data);
brainModel.two_d_cleaning(data);

# Remove empty rows/columns and plains from 3D array
brainModel.trim_mesh(data)

# Find and fill voids within model
print("########## Removing voids from data ##########")
solver = Maze_Solver(data);
data = solver.find_and_fill_voids();

# # Create CSF layer around GM
# print("########## Adding layers of CSF ##########")
# brainModel.add_CSF(data,layers=1);
          
# Create point cloud
print("########## Creating point cloud from data ##########")
pointCloud_full = PointCloud();
pc = pointCloud_full.create_point_cloud_of_data(data);
# # Data for visualization
# pointCloud.view_slice(0, 50); # View slice of point cloud about chosen axis
# pointCloud.view_point_cloud(); # View full 3D point cloud

# Create mesh from point cloud
print("########## Creating mesh from point cloud ##########")
mesh = Mesh(pointCloud_full.pcd,voxel_size)
mesh.create_edge_to_element_connectivity()
# #Global Smoothing
# # Smooth outer surface of mesh (including CSF)
# print("########## Smoothing global mesh ##########")
# iterations = 8
# coeffs = [0.4,-0.2]
# mesh.smooth_mesh(coeffs, iterations)

# # Smooth mesh (excluded CSF)
# print("########## Smoothing mesh exclusing CSF ##########")
# iterations = 8
# coeffs = [0.4,-0.2]
# mesh.smooth_mesh(coeffs, iterations, elementsNotIncluded=[24])

# # Optional Boundary Smoothing
# # Smooth white matter boundary

# print("########## Smoothing white matter ##########")
# non_whitematter_labels_map = material_labels.get_homogenized_labels_map()
# non_whitematter_labels_map.pop('WhiteMatter')
# non_whitematter_labels_map = list(non_whitematter_labels_map.values())
# iterations = 8
# coeffs = [0.4,-0.2]
# mesh.smooth_mesh(coeffs, iterations, elementsNotIncluded=non_whitematter_labels_map)

# # Smooth tumor boundary
# print("########## Smoothing tumour boundary ##########")
# non_lesion_labels_map = material_labels.get_homogenized_labels_map()
# non_lesion_labels_map.pop('Lesion')
# non_lesion_values = list(non_lesion_labels_map.values())
# iterations = 8
# coeffs = [0.4,-0.2]
# mesh.smooth_mesh(coeffs, iterations, elementsNotIncluded=non_lesion_values)

# Write mesh to file
test_labels_map  = Material_Label()
test_labels_map.addLabelToMap('edges_to_check',1000)
filename = "tester_mesh_cleanup"
fileType = "vtk"
print("########## Writing data to " + filename + " as a " + fileType.upper() + " file ##########")
mesh.write_to_file("C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\PythonScripts\\BrainMesher", 
                    filename, labels_map=test_labels_map, filetype=fileType);






