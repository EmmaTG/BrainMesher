# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:34:00 2023

@author: grife
"""

from BrainHexMesh import BrainHexMesh
from Mesh import Mesh, Element
from Material_Label import Material_Label

class ConfigFile():
    def __init__(self):
        self.readFile = True
        self.fileInPath = "C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\PythonScripts\\BrainMesher"
        self.fileIn = 'aseg_tumor.mgz'
        self.readData = False
        self.data = []        
        self.fileoutPath = "C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\PythonScripts\\BrainMesher"
        self.writeToFile = False
        self.fileout = "tester_ventricle_tumor"
        self.fileoutTypes = ['vtk',]
        self.Coarsen = True
        self.Add_CSF = False
        self.Smooth = False
        self.iterations = 2
        self.coeffs = [0.6,-0.4]
        
        self.Smooth_regions = []
        self.region_iterations = 4
        self.region_coeffs = [0.6,-0.4]
        
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
        self.material_labels.addLabelToMap('Ventricles' , [4,43,14,15]); # Lateral(L&R), 3rd, 4th
        if self.Add_CSF:
            self.material_labels.addLabelToMap('CSF' , [24]); # Left, Right
        
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

config = ConfigFile();
brainModel = BrainHexMesh();
brainModel.config(config)
data,pointCloud,mesh = brainModel.run();

tumor_labels = config.material_labels.get_homogenized_labels_map()
tumor_labels.pop("Lesion")
non_tumor_elements = list(tumor_labels.values())

elementsOnTB = mesh.locate_elements_on_boundary(elementsNotIncluded=non_tumor_elements)
for e in elementsOnTB:
    mesh.elements[e].properties['mat'].append(200)

boundaryElementsOnTB = mesh.locate_boundary_element_map(elementsNotIncluded=non_tumor_elements)
boundary_elements_map = {}
element_number = max(list(mesh.elements.keys()))
for compoundKey,ica in boundaryElementsOnTB.items():
    element_number += 1
    boundary_elements_map[element_number] = Element(element_number, ica, mat=[300])

new_material_labels  = Material_Label()
new_material_labels.addLabelToMap("Tumor elements", 200)
new_material_labels.addLabelToMap("Tumor Boundary elements", 300)
new_material_labels.addLabelToMap("Tumor ", 25)
brainModel.write_to_file(mesh, material_labels=new_material_labels, boundaryElementMap=boundary_elements_map)

