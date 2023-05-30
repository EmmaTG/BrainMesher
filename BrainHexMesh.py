"""
Created on Wed May 10 09:11:56 2023

@author: grife
"""

import numpy as np
# import mne
# import pooch
import nibabel
from BrainModel import BrainModel
from Maze_Solver import Maze_Solver
from PointCloud import PointCloud
from Mesh import Mesh

# Step 1: Using freesurfer and 'recon-all' create mri outputs. Ensure aseg.mgz is created.
t1_file = 'aseg.mgz'
t1 = nibabel.load(t1_file)
# t1.orthoview()
data = np.asarray(t1.dataobj)

data_test = data
# data = data_test
# data = np.ones((6,6,6))
# data[2,0,2] = 0
# data[2,1,2] = 0
# data[2,2,2] = 0
brainModel = BrainModel()

# Step 2: Determine segmentation of brain model via labels map
brainModel.addLabelToMap('BrainStem', 16)
brainModel.addLabelToMap('GreyMatter', [3,42]) # Left, Right
brainModel.addLabelToMap('WhiteMatter' , [2,41,77]); # Left, Right, WM-hypointensities
brainModel.addLabelToMap('Corpuscallosum' , [251,252,253,254,255]); # CC_Posterior, CC_Mid_Posterior, CC_Central, CC_Mid_Anterior, CC_Anterior
brainModel.addLabelToMap('BasalGanglia' , [11,50,12,51,13,52,26,58,62,30]); # Caudate(L&R), Putamen(L&R), Palladium(L&R), Accumbens Area(L&R), vessel(L&R)
brainModel.addLabelToMap('Cerebellum' , [7,46,8,47]); # WM(L&R), GM(L&R)
brainModel.addLabelToMap('Thalamus' , [10,49,28,60]); # Thalamus(L&R), Ventral DC(L&R)
brainModel.addLabelToMap('Hippocampus' , [17,53]); # Left, Right
brainModel.addLabelToMap('Amygdala' , [18,54]); # Left, Right

# brainModel.addLabelToMap('Left-Lateral-Ventricle' , [4]);
# brainModel.addLabelToMap('Right-Lateral-Ventricle' , [43]);
# brainModel.addLabelToMap('Left-Inf-Lat-Vent' , [5]);
# brainModel.addLabelToMap('Right-Inf-Lat-Vent ' , [44]);
# brainModel.addLabelToMap('3rd-Ventricle' , [14]);
# brainModel.addLabelToMap('4th-Ventricle' , [15]);
# brainModel.addLabelToMap('Left-choroid-plexus' , [31]);
# brainModel.addLabelToMap('Right-choroid-plexus' , [63]);
# brainModel.addLabelToMap('Optic-Chiasm' , [85]);

# Homogenize labels
data = brainModel.homogenize_labels(data);

# Coarsen the brain model
voxel_size = 2
data = brainModel.coarsen(voxel_size, data)

# Clean image removing isolated pixels and small holes
brainModel.clean_mesh_data(data);
brainModel.two_d_cleaning(data);
# Remove empty rows/columns and plains from 3D array
brainModel.trim_mesh(data)

# Find and fill voids within model
solver = Maze_Solver(data);
data = solver.find_and_fill_voids();

# # Create point cloud
pointCloud = PointCloud();
pc = pointCloud.create_point_cloud_of_data(data);
# # # Data for visualization
# # pointCloud.view_slice(1, 50); # View slice of point cloud about chosen axis
# pointCloud.view_point_cloud(); # View full 3D point cloud

# Create mesh from point cloud
mesh = Mesh(pc,voxel_size)

Get boundary quads
mesh.locate_boundary_faces()

Smooth mesh
iterations = 6
coeffs = [0.6,-0.2]
mesh.smooth_mesh(coeffs, iterations)

# Write mesh to file
mesh.write_to_file("C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\PythonScripts\\BrainMesher", "tester");






