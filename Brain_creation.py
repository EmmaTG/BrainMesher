# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:34:00 2023

@author: grife
"""

from BrainHexMesh import BrainHexMesh
from Config import ConfigFile;

path = "C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\PythonScripts\\BrainMesher"
fileIn = '\\mri\\aseg_tumor.mgz'
config = ConfigFile(path,fileIn);

brainCreator = BrainHexMesh();
brainCreator.config(config)
config.openConfigFile();

data = brainCreator.import_data(config.fileInPath,config.fileIn);
data = brainCreator.preprocess(data, lesion=True);

pointCloud = brainCreator.make_point_cloud(data)

mesh = brainCreator.make_mesh(pointCloud.pcd);
brainCreator.clean_mesh(mesh);
brainCreator.centerMesh(mesh);

all_labels = config.material_labels.get_homogenized_labels_map()
label_for_unused = all_labels.get("Unused")
mesh.remove_region(label_for_unused)
        
if config.Smooth:
    mesh = brainCreator.smooth_mesh(mesh);       

brainCreator.write_to_file(mesh)  
config.closeConfigFile();    



