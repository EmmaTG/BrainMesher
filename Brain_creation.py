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
        self.fileout = "tester_brain_tumor_1"
        self.fileoutTypes = ['vtk','ucd'] # 'ucd' | 'vtk' | 'abaqus'
        self.Coarsen = True
        self.Add_CSF = False
        self.Smooth = False
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
        # self.material_labels.addLabelToMap('CSF' , [24,4,43,14,15]); # Left, Right, Lateral(L&R), 3rd, 4th ventricles
        
        # Unused labels (will be set to 0)
        self.material_labels.addLabelToMap('Ventricles' , [4,43,14,15,31,63,85]); # Lateral(L&R), 3rd, 4th
        self.material_labels.addLabelToMap('CSF' , [24]); # Left, Right, Lateral(L&R), 3rd, 4th ventricles
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

e_centroids = np.zeros((max(mesh.elements.keys())+1,4))
count = 0
for e_num,element in mesh.elements.items():
    element_centroid = element.calculate_element_centroid(mesh.nodes)
    e_centroids[e_num] = list(element_centroid) + [element.properties['mat'][0]]
    count += 1
e_centroids = np.stack(e_centroids, axis = 0)

grey_matter_label = config.material_labels.get_homogenized_labels_map()['GreyMatter']
white_matter_label = config.material_labels.get_homogenized_labels_map()['WhiteMatter']

mesh.create_node_to_element_connectivity()
boundary_elements_map = {}
boundary_number = max(list(mesh.elements.keys()))
boundaryElementsOnCSF = mesh.locate_boundary_element_map() 
element_to_boundary = {}
count1 = 0
count2 = 0
count3 = 0
print("Locating CSF boundary")
for compoundKey,ica in boundaryElementsOnCSF.items():
        [element_num,face] = [int(x) for x in compoundKey.split("-")]
        [xc,yc,zc,m] = e_centroids[element_num] 
        Boundary = False
        if m == 24:
            elements_inline = e_centroids
            elements_inline = elements_inline[np.where(elements_inline[:,0]==xc)[0],:]
            elements_inline = elements_inline[np.where(elements_inline[:,2]==zc)[0],:]
            elements_inline = elements_inline[elements_inline[:, 1].argsort()]
            current_element_idx, = np.where(elements_inline[:,1]>=yc)
            if len(current_element_idx)!= 0:
                count1 += 1
                current_element_idx = current_element_idx[0]
                elements_inline = elements_inline[:current_element_idx]
            mats = list(elements_inline[:,3])
            if (mats.count(24) + mats.count(grey_matter_label) + mats.count(white_matter_label)) == len(mats):
                Boundary = True;
            else:
                mats = list(filter(lambda x:x!=24,mats))
                if mats[-1] == 3:
                    Boundary = True
        if Boundary:
            boundary_number += 1
            boundary_elements_map[boundary_number] = QuadElement(boundary_number, ica, mat=[400])
            # else:
            #     count += 1
        # if not element_to_boundary.__contains__(element_num):
        #     [xc,yc,zc,m] = e_centroids[element_num]            
        #     Boundary = True
        #     if m == 24 and yc > 160:
        #         elements_inline = np.stack(e_centroids, axis = 0)
        #         elements_inline = elements_inline[np.where(elements_inline[:,0]==xc)[0],:]
        #         elements_inline = elements_inline[np.where(elements_inline[:,2]==zc)[0],:]
        #         elements_inline = elements_inline[elements_inline[:, 1].argsort()]
        #         current_element_idx, = np.where(elements_inline[:,1]>yc)
        #         Boundary = False
        #         if len(current_element_idx)>0:
        #             current_element_idx = current_element_idx[0]
        #             elements_inline = elements_inline[:current_element_idx]
        #             mats = list(elements_inline[:,3])
        #             if (mats.count(24) + mats.count(0)) == len(mats):
        #                 Boundary = True;
        #             else:
        #                 if (mats.count(24) + mats.count(0) + mats.count(grey_matter_label) + mats.count(white_matter_label)) != len(mats):
        #                     for m in mats:
        #                         if m != 24:
        #                             if m == grey_matter_label:
        #                                 Boundary = True
        #                             else:
        #                                 print("not_included")
        #                             break;                                     
        #                 else:
        #                     Boundary = True
        #         else:
        #             Boundary = True 
        #         element_to_boundary[element_num] = Boundary
        # else:
        #     Boundary = element_to_boundary[element_num]
        # if Boundary:
        #     boundary_number += 1
        #     boundary_elements_map[boundary_number] = QuadElement(boundary_number, ica, mat=[400])
config.material_labels.addLabelToMap("CSF elements", 400)

print("Locating Tumor boundary")
Non_lesion_lables = config.material_labels.get_homogenized_labels_map()
lesion_label = Non_lesion_lables.pop("Lesion")
boundaryElementsLesion = mesh.locate_boundary_element_map(elementsNotIncluded=list(Non_lesion_lables.values()))
for compoundKey,ica in boundaryElementsLesion.items():
    boundary_element = QuadElement(boundary_number, ica, mat=[500])
    [element_num,face] = [int(x) for x in compoundKey.split("-")]
    [xc,yc,zc,m] = e_centroids[element_num] 
    element_centroid = [xc,yc,zc]
    face_centroid = boundary_element.calculate_element_centroid(mesh.nodes)
    normal =  np.array(face_centroid) -np.array([xc,yc,zc])
    normal = np.array([int(n) for n in (normal/np.linalg.norm(normal))])
    absnormal = np.array([abs(n) for n in normal])
    direc, = np.where(absnormal == 1)[0]
    dims = [0,1,2]
    dims.remove(direc)
    elements_inline = e_centroids
    for d in dims:
        elements_inline = elements_inline[np.where(elements_inline[:,d]==element_centroid[d])[0],:]
    elements_inline = elements_inline[elements_inline[:, direc].argsort()]
    add_element = True
    if normal[direc] == -1:
        current_element_idx, = np.where(elements_inline[:,direc]<=element_centroid[direc])
        if len(current_element_idx)!= 0:
            current_element_idx = current_element_idx[-1]
            elements_inline = elements_inline[:current_element_idx+1]
        else:
            add_element = False
            print("Error normal -1")
    elif normal[direc] == 1:
        current_element_idx, = np.where(elements_inline[:,direc]>=element_centroid[direc])
        if len(current_element_idx)!= 0:
            current_element_idx = current_element_idx[0]
            elements_inline = elements_inline[current_element_idx:]
        else:
            add_element = False
            print("Error normal +1")
    else:
        add_element = False
        print("Error")
        
    if add_element:
       add_element = False
       mats = list(elements_inline[:,3])
       mats = list(filter(lambda x:x!=24,mats))
       # if normal[direc] == 1:
       #     print(mats)
       if len(mats)>5:
           add_element = True
    if add_element:
        boundary_number += 1
        boundary_elements_map[boundary_number] = boundary_element
config.material_labels.addLabelToMap("TumorBoundary", 500)    

all_labels = config.material_labels.get_homogenized_labels_map()
label_for_ventricles = all_labels.get("Ventricles")

element_keys = list(mesh.elements.keys())    
for element_num in element_keys:
    e = mesh.elements[element_num]
    if e.properties['mat'].count(label_for_ventricles):
        mesh.delete_element(element_num)
        
if config.Smooth:
    mesh = brainModel.smooth_mesh(mesh);
  
# element_keys = list(mesh.elements.keys())    
# for element_num in element_keys:
#     e = mesh.elements[element_num]
#     if e.properties['mat'].count(lesion_label):
#         mesh.delete_element(element_num)       
      
# config.material_labels.removeLabel("Ventricles")
# config.material_labels.removeLabel("Lesion")
brainModel.write_to_file(mesh, material_labels=config.material_labels, boundaryElementMap=boundary_elements_map)

# mesh.write_to_file("C:\\Users\grife\OneDrive\Documents\PostDoc\BrainModels\PythonScripts\BrainMesher", "Tester_brain_tumor",
#                    config.material_labels, filetype='vtk', boundaryElementMap=boundary_elements_map)


