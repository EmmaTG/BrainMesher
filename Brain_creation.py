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
        self.fileIn = 'aseg.mgz'
        self.readData = False
        self.data = []        
        self.fileoutPath = "C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\PythonScripts\\BrainMesher"
        self.writeToFile = True
        self.fileout = "tester_brain_Silvia_H_and_M_4"
        self.fileoutTypes = ['vtk'] # 'ucd' | 'vtk' | 'abaqus'
        self.Coarsen = True
        self.Add_CSF = True
        self.Smooth = True
        self.iterations = 4
        self.coeffs = [0.6,-0.4]
        
        self.Smooth_regions = []
        self.region_iterations = [6]
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
        self.material_labels.addLabelToMap('CSF' , [24,4,43,14,15]); # Left, Right, Lateral(L&R), 3rd, 4th ventricles
        
        # Unused labels (will be set to 0)
        # self.material_labels.addLabelToMap('Ventricles' , [4,43,14,15]); # Lateral(L&R), 3rd, 4th
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
# import numpy as np
# dim = 20
# data = np.ones((dim,dim,dim))*3
# data[:,:,dim-1] = np.ones((dim,dim))*24
# data[10:16,10:16,10:16] = np.ones((6,6,6))*25
# config.add_data(data);


brainModel = BrainHexMesh();
brainModel.config(config)
data = brainModel.import_data();
data = brainModel.preprocess(data);


pointCloud = brainModel.make_point_cloud(data)
mesh = brainModel.make_mesh(pointCloud.pcd);

boundary_elements_map = {}
# element_number = max(list(mesh.elements.keys()))
# mesh.create_node_to_element_connectivity()

# # if config.Add_CSF:
# boundaryElementsOnCSF = mesh.locate_boundary_element_map()
# for compoundKey,ica in boundaryElementsOnCSF.items():
#         [element_num,face] = [int(x) for x in compoundKey.split("-")]
#         if mesh.elements[element_num].properties['mat'].count(24) or mesh.elements[element_num].properties['mat'].count(16):
#             element_number += 1
#             boundary_elements_map[element_number] = Element(element_number, ica, mat=[400])
# config.material_labels.addLabelToMap("CSF elements", 400)


# element_to_delete = []
# tumor_labels = config.material_labels.get_homogenized_labels_map()
# label_for_tumor = tumor_labels.pop("Lesion")
# non_tumor_elements = list(tumor_labels.values())
# boundaryElementsOnTB = mesh.locate_boundary_element_map(elementsNotIncluded=non_tumor_elements)
# for compoundKey,ica in boundaryElementsOnTB.items():
#     [element_num,face] = [int(x) for x in compoundKey.split("-")]
#     node1 = np.array(mesh.nodes[ica[0]])
#     node2 = np.array(mesh.nodes[ica[1]])
#     node3 = np.array(mesh.nodes[ica[2]])
#     v12 = node2-node1
#     v13 = node3-node1
#     normal = np.cross(v12,v13)
#     normal = np.absolute(normal/np.linalg.norm(normal))
#     pcd_data = mesh.elementToPointCloud[element_num]
#     direc, = np.where(normal == 1)
#     if len(direc) != 1:
#         print("Error with normal")
#         print(normal)
#     direc = direc[0]
#     dims = [0,1,2]
#     direcs = np.delete(dims,direc)
#     pc_points = pointCloud.pcd
#     for d in direcs:
#         r,c = np.where(pc_points==pcd_data[d])
#         pc_points = pc_points[r[np.where(c == d)]]
#     pc_points = sorted(pc_points, key=lambda x:x[direc])
#     empty_before = True
#     empty_after = True
#     if(pc_points[0][direc]<(pcd_data[direc]-3)) and (pc_points[-1][direc]>(pcd_data[direc]+3)):
#         if (pcd_data[direc] != pc_points[0][direc]) and (pcd_data[direc] != pc_points[-1][direc]):        
#             for point in pc_points:
#                 if point[direc]<pcd_data[direc] and empty_before:
#                     if point[3] != 24:
#                         empty_before = False
#                 elif point[direc]>pcd_data[direc]:
#                     if not empty_before:
#                         if point[3] != 24:
#                             empty_after = False
#                             break
#                     else:
#                         break
#         if not empty_after and not empty_after:
#             element_number += 1
#             boundary_elements_map[element_number] = Element(element_number, ica, mat=[300])
#             element_to_delete.append(element_num)
# config.material_labels.addLabelToMap("Tumor Boundary elements", 300)


# for element_num, e in mesh.elements.items():
#     if e.properties['mat'] == [label_for_tumor]:
#         element_to_delete.append(element_num)        
# for element_num in element_to_delete:
#     mesh.delete_element(element_num)
    
if config.Smooth:
    mesh = brainModel.smooth_mesh(mesh);
brainModel.write_to_file(mesh, material_labels=config.material_labels, boundaryElementMap=boundary_elements_map)


