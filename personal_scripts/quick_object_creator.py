# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 09:23:11 2023

@author: grife

Script to make easy geometry using brainMesher code

"""

import numpy as np
from point_cloud.PointCloud import PointCloud
from mesh.Mesh import Mesh
from writers.HeterogeneityConverter import FourRegionConverter
from writers.Writer import Writer
from mesh.refinement.Refiner import Refiner
from config.Material_Label import Material_Label
from dotenv import load_dotenv
from readers.Reader import ABQReader, VTKReader
from file_converters.Converter import Converter
import os
# filenames = [4,5,6,9,10]
# for curve_number in filenames:
#     filename = "Brain_sphere_3mm_30mm_I6mm_curve{}_Loc1_9R_VTK".format(curve_number)
    # filename_out = "Brain_sphere_3mm_30mm_I6mm_curve{}_Loc1_4R".format(curve_number)
    # path = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/CoreFoamModels"



filename = "brain_full_csf_centered_VTK"
filename_out = "brain_full_csf_centered_rotated"
path = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/PythonScripts/BrainMesher/IOput/out/"
path = "C:/Users/grife/OneDrive/Documents/PostDoc/Students/Yashasvi/meshes"
pathout = "C:/Users/grife/OneDrive/Documents/PostDoc/Students/Yashasvi/meshes"

reader = VTKReader()
reader.openReader(filename, path)
mesh = reader.getMesh()
reader.closeReader()

import mesh.mesh_transformations as mt
mt.translate_mesh(mesh.nodes,[2,-8,6])
mt.rotate_mesh(mesh.nodes,axis=1, degrees=90)

# mat_converter = FourRegionConverter()
# mat_converter.convert_materials_labels(mesh)

writer = Writer()
writer.openWriter('ucd',filename_out,pathout)
writer.writeMeshData(mesh)
writer.closeWriter()

writer = Writer()
writer.openWriter('vtk',filename_out,pathout)
writer.writeMeshData(mesh)
writer.closeWriter()

# reader = VTKReader()
# mesh = reader.getMesh()
# reader = ABQReader()
# reader.openReader(filename, path)
# reader.readNodes()
# reader.readElements()
# mesh = reader.getMesh()
# print("Done")
#
# writer = Writer()
# writer.openWriter('vtk', "ABQ_reader_test_To_VTK", path + "/")
# writer.writeMeshData(mesh)
# writer.closeWriter()

# load_dotenv()
#
# home = os.getenv('HOME')
# data_home = home
#
# labels = Material_Label()
# labels.addLabelToMap("Body1", 3)
# filename = "output_full_brain_VTK"
# filename_out = "output_non_regular_grid_centered"
# path = "C:/Users/grife/OneDrive/Documents/PostDoc/Students/Yashasvi/meshes"
# dimension = 20
#
# data = np.ones((dimension,dimension,dimension))*3
# data[:,dimension-2:,:2] = np.ones((dimension,1,1))*0
# data[:8,:8,dimension-2:] = np.ones((8,8,2))*0
#
# pointCloud = PointCloud()
# pointCloud.create_point_cloud_from_voxel(data)
# mesh = Mesh()
# mesh.create_mesh_from_Point_Cloud(pointCloud.pcd, 2)
# mesh.center_mesh_by_region()
#
# writer = Writer()
# writer.openWriter('vtk',filename_out,path)
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
        
    
    
    
    


