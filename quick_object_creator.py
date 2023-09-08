# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 09:23:11 2023

@author: grife

Script to make easy geometry using brainMesher code

"""

import numpy as np
from point_cloud.PointCloud import PointCloud
from mesh.Mesh import Mesh
from writers.Writer import Writer
from mesh.refinement.Refiner import Refiner
from config.Material_Label import Material_Label
from dotenv import load_dotenv
import os

load_dotenv()

home = os.getenv('HOME')
data_home = home

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
    writer.openWriter(file_type, "Tester_block_atrophy_small", "/".join([home, 'Models']))
    writer.writeMeshData(mesh)
    writer.closeWriter()
        
    
    
    
    


