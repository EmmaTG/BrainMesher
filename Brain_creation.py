# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:34:00 2023

@author: grife
"""

from BrainHexMesh import BrainHexMesh
from Mesh import Mesh, QuadElement, HexElement
from Material_Label import Material_Label
import numpy as np

class ConfigFile():
    def __init__(self):
        self.readFile = True
        self.fileInPath = "C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\PythonScripts\\BrainMesher"
        self.fileIn = 'aseg_tumor.mgz'
        self.readData = False
        self.data = []        
        self.fileoutPath = "C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\PythonScripts\\BrainMesher"
        self.writeToFile = True
        self.fileout = "tester_brain_tumor"
        self.fileoutTypes = ['vtk','ucd'] # 'ucd' | 'vtk' | 'abaqus'
        self.Coarsen = True
        self.Add_CSF = True
        self.Smooth = True
        self.iterations = 6
        self.coeffs = [0.6,-0.4]
        
        self.Smooth_regions = ['Lesion']
        self.region_iterations = [4]
        self.region_coeffs =[[0.6,-0.4]]
        
        self.material_labels  = Material_Label()
        self.material_labels.addLabelToMap('BrainStem', 16)
        self.material_labels.addLabelToMap('GreyMatter', [3,42]) # Left, Right
        self.material_labels.addLabelToMap('WhiteMatter' , [2,41,77]); # Left, Right, WM-hypointensities
        self.material_labels.addLabelToMap('Corpuscallosum' , [251,252,253,254,255]); # CC_Posterior, CC_Mid_Posterior, CC_Central, CC_Mid_Anterior, CC_Anterior
        self.material_labels.addLabelToMap('BasalGanglia' , [11,50,12,51,13,52,26,58,62,30]); # Caudate(L&R), Putamen(L&R), Palladium(L&R), Accumbens Area(L&R), vessel(L&R)
        self.material_labels.addLabelToMap('Cerebellum' , [7,46,8,47]); # WM(L&R), GM(L&R)
        self.material_labels.addLabelToMap('Thalamus' , [10,49,28,60]); # Thalamus(L&R), Ventral DC(L&R)
        self.material_labels.addLabelToMap('Hippocampus' , [17,53]); # Left, Right
        self.material_labels.addLabelToMap('Amygdala' , [18,54]); # Left, Right
        self.material_labels.addLabelToMap('Lesion' , [25,57]); # Left, Right
        self.material_labels.addLabelToMap('CSF' , [24]); # Left, Right, Lateral(L&R), 3rd, 4th ventricles
        self.material_labels.addLabelToMap('Unused' , [4,43,14,15,31,63,85]); # Lateral(L&R), 3rd, 4th NOTE: These labels will be removed to create empty spaces!!
        # OR
        # material_labels.addLabelToMap('Left-Lateral-Ventricle' , [4]);
        # material_labels.addLabelToMap('Right-Lateral-Ventricle' , [43]);
        # material_labels.addLabelToMap('Left-Inf-Lat-Vent' , [5]);
        # material_labels.addLabelToMap('Right-Inf-Lat-Vent ' , [44]);
        # material_labels.addLabelToMap('3rd-Ventricle' , [14]);
        # material_labels.addLabelToMap('4th-Ventricle' , [15]);
        #
        # material_labels.addLabelToMap('Left-choroid-plexus' , [31]);
        # material_labels.addLabelToMap('Right-choroid-plexus' , [63]);
        # material_labels.addLabelToMap('Optic-Chiasm' , [85]);
        
    def add_data(self, importedData):
        self.readFile = False
        self.readData = True
        self.data = importedData


config = ConfigFile();

brainModel = BrainHexMesh();
brainModel.config(config)

data = brainModel.import_data(config.fileInPath,config.fileIn);
data = brainModel.preprocess(data);
pointCloud = brainModel.make_point_cloud(data)
mesh = brainModel.make_mesh(pointCloud.pcd);

all_labels = config.material_labels.get_homogenized_labels_map()
label_for_unused = all_labels.get("Unused")

element_keys = list(mesh.elements.keys())    
for element_num in element_keys:
    e = mesh.elements[element_num]
    if e.getMaterial().count(label_for_unused):
        mesh.delete_element(element_num)
        
if config.Smooth:
    mesh = brainModel.smooth_mesh(mesh);       
      
brainModel.write_to_file(mesh, material_labels=config.material_labels, boundaryElementMap=boundary_elements_map)


