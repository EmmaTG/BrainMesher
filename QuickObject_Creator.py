# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 09:23:11 2023

@author: grife

Script to make easy geometry using brainMesher code

"""

import numpy as np
from scipy import ndimage
from GridBox import GridBox
from PointCloud import PointCloud
from Mesh import Mesh
from Element import Element
from Writer import Writer
from Refine.Refiner import Refiner

from Material_Label import Material_Label


labels = Material_Label()
labels.addLabelToMap("Body1", 3)
dimension = 20;

data = np.ones((dimension,dimension,dimension))*3

pointCloud = PointCloud()
pointCloud.create_point_cloud_from_voxel(data)
mesh = Mesh()
mesh.create_mesh_from_Point_Cloud(pointCloud.pcd, 2)

meshRefiner = Refiner(mesh)
bounds = [17, 23, 37, 41, 21, 27]
meshRefiner.refine_within_region(bounds)
meshRefiner.refine_elements([1,2])
meshRefiner.refine_around_point([0,26,24],5)

# boundary_elements_map = {}
# boundary_number = max(list(mesh.elements.keys()))
# boundaryElements = mesh.locate_boundary_element_map()
# for compoundKey,ica in boundaryElements.items():
#     boundary_number += 1
#     ica_nodes = [mesh.nodes[n] for n in ica]
#     boundary_element = Element(boundary_number, ica_nodes, mat=[200])
#     face_centroid = boundary_element.calculate_element_centroid();
#     [element_num,face] = [int(x) for x in compoundKey.split("-")]
#     centroid = mesh.elements[element_num].calculate_element_centroid();
#     normal = np.array(face_centroid) - np.array(centroid)
#     normal = normal/np.linalg.norm(normal)
#     direc, = np.where(normal == 1)
#     negDirec, = np.where(normal == -1)
#     if not list(direc).count(0):
#         boundary_elements_map[boundary_number] = boundary_element
#     else:
#         boundary_number -= 1
#
# mesh.addBoundaryElements(boundary_elements_map);
#
# labels.addLabelToMap("FixBoundary", 200)
#
# mesh.dataToWrite.append("concentration")
# center = np.array([20,20,20]);
# concentration_radius = 5
# node_touched = []
# for e in mesh.elements.values():
#     centroid = e.calculate_element_centroid();
#     radius = np.linalg.norm(center-centroid);
#     if (radius <= concentration_radius):
#         c = 1 - (radius/concentration_radius);
#         for n in e.ica:
#             if (node_touched.count(n.number)):
#                 if mesh.nodes[n.number].getData("concentration")[0]<c:
#                     mesh.nodes[n.number].addData("concentration",[c]);
#             else:
#                 mesh.nodes[n.number].addData("concentration",[c]);
#             node_touched.append(n.number)
#     else:
#         for n in e.ica:
#             if not (node_touched.count(n.number)):
#                 mesh.nodes[n.number].addData("concentration",[0]);
        

for file_type in ['vtk']:
    writer = Writer()
    writer.openWriter(file_type, "Tester_block_atrophy_small", "C:\\Users\grife\OneDrive\Documents\PostDoc\BrainModels\PythonScripts\BrainMesher\Models\\")
    writer.writeMeshData(mesh)
    writer.closeWriter()

# labels = Material_Label()
# labels.addLabelToMap('Lesion' , [25,57])
# labels.addLabelToMap('Outer' , [1])

# fileInPath = "C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\PythonScripts\\BrainMesher"
# fileIn = 'aseg_tumor.mgz'

# brainModel = BrainModel();
# data = brainModel.import_file(fileInPath, fileIn);
# for x in range(66,113):
#     for y in range(82,125):
#         for z in range(62,110):
#             if data[x,y,z] != 25:
#                 data[x,y,z] = 1
# data = labels.homogenize_material_labels(data)
# data = brainModel.coarsen(2, data)
# brainModel.clean_voxel_data(data)

# newData = brainModel.create_binary_image(data, search=25)
# print("Lesion element size before: {}".format(np.sum(newData)))


# structure1 = ndimage.generate_binary_structure(3,1)
# structure2 = ndimage.generate_binary_structure(3,3)

# newData = ndimage.binary_dilation(newData, structure=structure2, iterations=3).astype(int)
# newData = ndimage.binary_erosion(newData, structure=structure1, iterations=3).astype(int)
# newData = ndimage.binary_erosion(newData, structure=structure2, iterations=1).astype(int)
# newData = ndimage.binary_dilation(newData, structure=structure2, iterations=1).astype(int)
# newData = ndimage.binary_erosion(newData, structure=structure2, iterations=1).astype(int)

# hit_structure1 = np.ones((2,2,2))
# hit_structure1[0,1,0] = 0
# hit_structure2 = np.rot90(hit_structure1)
# hit_structure3 = np.rot90(hit_structure1, k= 2)
# hit_structure4 = np.rot90(hit_structure1, k= 3)
# hit_structure5 = np.rot90(hit_structure1,axes=(1,2))
# hit_structure6 = np.rot90(hit_structure5, k= 1)
# hit_structure7 = np.rot90(hit_structure5, k= 2)
# hit_structure8 = np.rot90(hit_structure5, k= 3)
# total_count = 0;
# count = 1
# iteration = 0
# while (count > 0 and iteration < 10):
#     iteration += 1
#     count = 0
#     count += brainModel.hit_and_miss_3d_2x2x2(newData, hit_structure1, fill=1)
#     count += brainModel.hit_and_miss_3d_2x2x2(newData, hit_structure2, fill=1)
#     count += brainModel.hit_and_miss_3d_2x2x2(newData, hit_structure3, fill=1)
#     count += brainModel.hit_and_miss_3d_2x2x2(newData, hit_structure4, fill=1)
#     count += brainModel.hit_and_miss_3d_2x2x2(newData, hit_structure5, fill=1)
#     count += brainModel.hit_and_miss_3d_2x2x2(newData, hit_structure6, fill=1)
#     count += brainModel.hit_and_miss_3d_2x2x2(newData, hit_structure7, fill=1)
#     count += brainModel.hit_and_miss_3d_2x2x2(newData, hit_structure8, fill=1)
#     total_count += count
# print("Lesion cleaned after {} iterations and {} elements added".format(iteration,total_count))
        
# count_removed = 0
# current_dimensions = data.shape
# for x in range(current_dimensions[0]):
#     if (np.any(newData[x,:,:] == 1) or np.any(data[x,:,:] == 25)):
#         for y in range(current_dimensions[1]):
#             if (np.any(newData[x,y,:] == 1) or np.any(data[x,y,:] == 25)):
#                 for z in range(current_dimensions[2]):
#                     if newData[x,y,z] == 1 and data[x,y,z] != 25 and data[x,y,z] != 0:
#                         data[x,y,z] = 25
#                     if newData[x,y,z] == 0 and (data[x,y,z] == 25):
#                         box = GridBox(data,[x,y,z])
#                         count_removed += 1
#                         idx, = np.where(box.gridBox == 25)
#                         box.gridBox = np.delete(box.gridBox,idx)
#                         replacement_value = box.mode()
#                         data[x,y,z] = replacement_value
                        
# print("previous lesion replaced with non-lesion: {}".format(count_removed))
# print("Lesion element size after: {}".format(np.sum(newData)))

# brainModel.trim_mesh(data)

# x,y,z = np.where(data == 1)
# for xf in range(min(x)-2,max(x)+3):
#     for yf in range(min(y)-2,max(y)+3):
#         for zf in range(min(z)-2,max(z)+3):
#             if data[xf,yf,zf] == 0:
#                 data[xf,yf,zf] = 24

# labels.addLabelToMap('CSF' , [24])


# pointCloud = PointCloud()
# pointCloud.create_point_cloud_from_voxel(data)
# mesh = Mesh()
# mesh.create_mesh_from_Point_Cloud(pointCloud.pcd, 2)

# boundary_elements_map = {}
# boundary_number = max(list(mesh.elements.keys()))
# boundaryElementsOnCSF = mesh.locate_boundary_element_map()
# for compoundKey,ica in boundaryElementsOnCSF.items():
#         [element_num,face] = [int(x) for x in compoundKey.split("-")]
#         face_centroid = calculate_element_centroid(ica,mesh.nodes)
#         element = mesh.elements[element_num]
#         centroid = calculate_element_centroid(element.ica,mesh.nodes)
#         normal = np.array(face_centroid )- np.array(centroid)
#         normal = normal/np.linalg.norm(normal)
#         direc, = np.where(normal == -1)
#         if list(direc).count(0):
#             if mesh.elements[element_num].properties['mat'].count(24):
#                 boundary_number += 1
#                 boundary_elements_map[boundary_number] = Element(boundary_number, ica, mat=[400])

# labels.addLabelToMap("CSF elements", 400)

# Non_lesion_lables = labels.get_homogenized_labels_map()
# lesion_label = Non_lesion_lables.pop("Lesion")
# boundaryElementsLesion = mesh.locate_boundary_element_map(elementsNotIncluded=list(Non_lesion_lables.values()))
# for compoundKey,ica in boundaryElementsLesion.items():
#     boundary_number += 1
#     boundary_elements_map[boundary_number] = Element(boundary_number, ica, mat=[500])

# labels.addLabelToMap("TractionBoundary", 500)

# # mesh.create_node_to_element_connectivity()
# # element_keys = list(mesh.elements.keys())
# # for e in element_keys:
# #     if mesh.elements[e].properties['mat']==[lesion_label]:
# #         mesh.delete_element(e);
        
# # mesh.smooth_mesh([0.6,-0.4], 6)
# for file_type in ['vtk','ucd']:
#     mesh.write_to_file("C:\\Users\grife\OneDrive\Documents\PostDoc\BrainModels\PythonScripts\BrainMesher", "Tester_lesion_original",
#                         labels, filetype=file_type, boundaryElementMap=boundary_elements_map)
        
    
    
    
    


