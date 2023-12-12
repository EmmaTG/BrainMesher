# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 09:23:11 2023

@author: grife

Script to make easy geometry using brainMesher code

"""

import numpy as np
from point_cloud.PointCloud import PointCloud
from mesh.Mesh import Mesh
from writers.HeterogeneityConverter import FourRegionConverter, MaterialsConverterFactory, Heterogeneity
from writers.Writers import Writer
from mesh.refinement.Refiner import Refiner
import mesh.mesh_transformations as mt
from config.Material_Label import Material_Label
from dotenv import load_dotenv
from readers.Reader import Reader
from mesh.PostProcessor import PostProcessor, CreateBoundaryElements
import common.helper_functions as hf
import os

path = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/CoreFoamModels"
filenames = [5,6,8,11]
for curve_number in filenames:
    filename = "sphere_in_cube_4mm_30mm_6mm_curve{}_VTK".format(curve_number)
#     filename_out = "Brain_sphere_3mm_30mm_I6mm_curve{}_Loc1_4R".format(curve_number)
    reader = Reader('vtk')
    reader.openReader(filename, path)
    mesh = reader.getMesh()
    reader.closeReader()
    print("Number of nodes: ", len(mesh.nodes))
    print("Number of elements: ", len(mesh.elements))


# filename = "brain_full_csf_centered_VTK"
# filename_out = "brain_full_csf_centered_rotated"
# path = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/PythonScripts/BrainMesher/IOput/out/"
# path = "C:/Users/grife/OneDrive/Documents/PostDoc/Students/Yashasvi/meshes"
# pathout = "C:/Users/grife/OneDrive/Documents/PostDoc/Students/Yashasvi/meshes"
#
# path = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/CoreFoamModels"
# filename = "Brain_incision_L12mm_W4mm_30mm_Loc1_4R"
# reader = Reader('ucd')
# reader.openReader(filename, path)
# mesh = reader.getMesh()
# reader.closeReader()
# #
# # import mesh.mesh_transformations as mt
# mt.translate_mesh(mesh.nodes,[-50, -18, +30])
# safety_factor = 0.1
# bounds_to_refine = [-1*safety_factor - (25 * 1.1),
#                     safety_factor + 50,  # (xmin,xmax  ...
#                     -1*safety_factor - 3,
#                     safety_factor + 3,  # ymin,ymax ...
#                     -1*safety_factor - 3,
#                     safety_factor + 3]  # zmin,zmax)
# refine_mesh = Refiner(mesh)
# refine_mesh.refine_within_region(bounds_to_refine)
# # mt.rotate_mesh(mesh.nodes,axis=1, degrees=90)
# #
# # # mat_converter = FourRegionConverter()
# # # mat_converter.convert_materials_labels(mesh)
# #
# writer = Writer()
# writer.openWriter('vtk',filename + "_refined",path)
# writer.writeMeshData(mesh)
# writer.closeWriter()
#
# writer = Writer()
# writer.openWriter('vtk',filename_out,pathout)
# writer.writeMeshData(mesh)
# writer.closeWriter()

# reader = Reader('vtk')
# mesh = reader.getMesh()
# path = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/CoreFoamModels"
# for filename_ext in ['curve8']:
# # for filename_ext in ['curve5', 'curve6', 'curve8', 'curve10', 'curve11']:
# #     filename = "sphere_in_cube_4mm_30mm_6mm_{}" .format(filename_ext)
#     filename = "Brain_sphere_4mm_30mm_I6mm_{}_Loc1_4R" .format(filename_ext)
#
#     reader = Reader('vtk')
#     reader.openReader(filename + "_VTK", path)
#     reader.readNodes()
#     reader.readElements()
#     mesh = reader.getMesh()
#     print("Done")
#
#     center_node = [50, 18, -30]
#     mt.translate_mesh(mesh.nodes, [-1 * center_node[0], -1 * center_node[1], -1 * center_node[2]])
#
#     # mat_converter = MaterialsConverterFactory.get_converter(Heterogeneity.FOURR)
#     # mat_converter.convert_materials_labels(mesh)
#
#     # mt.rotate_mesh(mesh.nodes, axis=1, degrees=90)
#
#     writer = Writer()
#     writer.openWriter('vtk', filename, path + "/")
#     writer.writeMeshData(mesh)
#     writer.closeWriter()
#     writer = Writer()
#     writer.openWriter('ucd', filename, path + "/" + filename)
#     writer.writeMeshData(mesh)
#     writer.closeWriter()

# load_dotenv()
#
# home = os.getenv('HOME')
# data_home = home
#
# labels = Material_Label()
# labels.addLabelToMap("Body1", 3)
# # filename = "output_full_brain_VTK"
# filename_out = "regular_grid_centered"
# path = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/Test_Models"
# dimension = 20
#
# data = np.ones((dimension,dimension,dimension))*3
# # data[:,dimension-2:,:2] = np.ones((dimension,1,1))*0
# # data[:8,:8,dimension-2:] = np.ones((8,8,2))*0
#
# pointCloud = PointCloud()
# pointCloud.create_point_cloud_from_voxel(data)
# mesh = Mesh()
# mesh.create_mesh_from_Point_Cloud(pointCloud.pcd, 2)
# mesh.center_mesh_by_region()
#
# post_processor = PostProcessor(None, mesh)
# post_processor = CreateBoundaryElements(post_processor, [])
# post_processor.post_process()
#
# remove_elements = []
# for e in mesh.boundaryElements.values():
#     coords1 = e.ica[0].getCoords()
#     coords2 = e.ica[1].getCoords()
#     coords3 = e.ica[2].getCoords()
#
#     v1v2 = np.array(coords1) - np.array(coords2)
#     v1v3 = np.array(coords1) - np.array(coords3)
#     normal = np.cross(v1v2, v1v3)
#     norm = normal / np.linalg.norm(normal)
#
#     centroid = e.calculate_element_centroid()
#
#     test_point = centroid + norm*(dimension/4)
#     bounds = [-dimension/2., dimension/2., -dimension/2., dimension/2., -dimension/2., dimension/2.]
#     if hf.value_in_square_bounds(test_point, bounds):
#         norm *= -1
#
#     if round(norm[0]) == -1:
#         e.setMaterial(2)
#     elif round(norm[0]) == 1:
#         e.setMaterial(1)
#     else:
#         remove_elements.append(e.num)
#
# for e_num in remove_elements:
#     mesh.boundaryElements.pop(e_num)
#
# writer = Writer()
# writer.openWriter('ucd',filename_out,path)
# writer.writeMeshData(mesh)
# writer.closeWriter()
#
# meshRefiner = Refiner(mesh)
# bounds = [17, 23, 37, 41, 21, 27]
# meshRefiner.refine_within_region(bounds)
# meshRefiner.refine_elements([1,2])
# meshRefiner.refine_around_point([0,26,24],5)

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
        

# for file_type in ['vtk']:
#     writer = Writer()
#     writer.openWriter(file_type, "Tester_block_atrophy_small", "/".join([home, 'Models']))
#     writer.writeMeshData(mesh)
#     writer.closeWriter()
        
    
    
    
    


