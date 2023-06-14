# -*- coding: utf-8 -*-
"""
Created on Thu May 25 15:39:38 2023

@author: grife
"""
import ABQ_UCD_handling as rw
import Smoothing as smooth
from Material_Label import Material_Label

class Element():
    
    def __init__(self, number, ica, **kwargs):        
        self.ica = ica
        self.num = number
        self.properties = {}
        for key,value in kwargs.items():
            self.properties[key] = value;
    
    
class Mesh():
    
    def __init__(self, pointData,voxel_size):
        self.elements = {}
        self.nodes = {}
        self.create_mesh(pointData,voxel_size)
        
    
    def create_mesh(self, pointData, voxel_size):
        [minX,minY,minZ] = pointData[:,:3].min(axis=0)
        [maxX,maxY,maxZ] = pointData[:,:3].max(axis=0)
        
        elementX = maxX-minX+1
        elementY = maxY+1
        elementZ = maxZ+1
        
        elementNo = 0;
        for p in pointData:
            elementNo += 1
            [x,y,z,m] = p
            startNode = (z+1) + (elementZ+1)*y + ((elementZ+1)*(elementY+1))*x
            element_ica = [int(startNode+1), int(startNode), int(startNode+(elementZ+1)) ,int(startNode+(elementZ+1)+1)]
            element_ica_tmp = []
            for i in element_ica:   
                newNode = int(i + (elementZ+1)*(elementY+1))
                element_ica_tmp.append(newNode)
                if not self.nodes.__contains__(i):
                    coords = self.calculate_node_coords(elementX,elementY,elementZ,i,voxel_size)
                    self.nodes[i] = coords
                if not self.nodes.__contains__(newNode):
                    coords = self.calculate_node_coords(elementX,elementY,elementZ,newNode,voxel_size)
                    self.nodes[newNode] = coords 
            element_ica += element_ica_tmp
            element = Element(elementNo, element_ica, mat=[m])
            self.elements[int(elementNo)] = element
        self.clean_mesh();
        self.remove_outer_white_matter(2, 3)
    
    def clean_mesh(self, elementsNotIncluded = [], replace=0):
        self.create_node_to_element_connectivity();
        self.clean_mesh_edges(elementsNotIncluded = elementsNotIncluded, replace=replace);
        self.clean_mesh_nodes(elementsNotIncluded = elementsNotIncluded, replace=replace);
         
    
    def locate_boundary_faces(self, elementsNotIncluded = []):
        print("Identifying boundary faces")
        elementMap = self.create_elements_map(elementsNotIncluded=elementsNotIncluded)
        self.boundary_element_map, self.elements_on_boundary, self.volume_elem_to_boundary = smooth.get_boundary_surfaces(elementMap)
        
    def remove_outer_white_matter(self, white_matter_label, replace_label):
        print("Cleaning brain boundary")
        self.locate_boundary_faces(elementsNotIncluded = [24]);
        for elem in self.elements_on_boundary:
            element = self.elements[elem]
            materials = element.properties['mat']
            if materials.count(white_matter_label):
                materials.remove(white_matter_label)
                materials.insert(0,replace_label)
        
        
    def clean_mesh_nodes(self, elementsNotIncluded = [], replace=0):
        print("Cleaning mesh")
        maxNodeNum = max(self.nodes.keys())
        newNodes = {}
        for key in self.nodes.keys():
            # connectedElements = self.nodeToElements[key]
            All_connectedElements = self.nodeToElements[key]
            connectedElements = []
            for conn_element in All_connectedElements:
                element = self.elements[conn_element]
                add = True
                for el_types in elementsNotIncluded:
                    if element.properties['mat'].count(el_types):
                        add=False
                        break
                if add:
                    connectedElements.append(conn_element)
            if len(connectedElements) == 2:
                element1 = self.elements[connectedElements[0]]
                element2 = self.elements[connectedElements[1]]
                numberSharedNodes = 0
                for node1 in element1.ica:
                    if element2.ica.count(node1)>0:
                        numberSharedNodes += 1
                assert numberSharedNodes>0
                if numberSharedNodes == 1:
                    if (replace != 0) or (len(elementsNotIncluded) != 0):
                        self.replace_element(element2.num,replace=replace)
                    else:
                        self.delete_element(element2.num)
                    # print("Cleaning mesh at node {} between elements {} and {}".format(key,element1.num,element2.num))                    
        #             maxNodeNum += 1
        #             nodeCoords = list(self.nodes[key])
        #             node_idx = element2.ica.index(key)
        #             element2.ica[node_idx] = maxNodeNum
        #             newNodes[maxNodeNum] = nodeCoords
        #             element2.properties['mat'].append(1000)
        #             # Remove from node to element map
        #             self.nodeToElements[key].remove(connectedElements[1])
        #             self.nodeToElements[maxNodeNum] = element2.num
        # for key,node in newNodes.items():
        #     self.nodes[key] = node;
    
    
    def delete_element(self, element_number):
        element = self.elements[element_number]
        ica = element.ica
        for n in ica:
            connectedElements = self.nodeToElements[n]
            connectedElements.remove(element_number)
            if len(connectedElements)==0:
                self.nodeToElements.pop(n)
                self.nodes.pop(n)
        self.elements.pop(element_number)
        
    def replace_element(self, element_number, replace=24):
        element = self.elements[element_number]
        element.properties['mat'] = [replace]
        # ica = element.ica
        # total_connected_elements_material = []
        # for n in ica:
        #     connectedElements = self.nodeToElements[n]
        #     for e in connectedElements:
        #         if self.elements[e].properties['mat'] != element.properties['mat']:
        #             total_connected_elements_material.append(self.elements[e].properties['mat'])
            
    
    def clean_mesh_edges(self, elementsNotIncluded = [], replace=0):
        print("Cleaning mesh edges")
        from Smoothing import get_element_faces
        edgesToElements = self.create_edge_to_element_connectivity(elementsNotIncluded)
        maxNodeNum = max(self.nodes.keys())
        newNodes = {}
        old_node_to_new = {}
        deleted_elements= []
        for edge, edgeConnectedElements in edgesToElements.items():
            if not deleted_elements.count(edgeConnectedElements[0]) and not deleted_elements.count(edgeConnectedElements[1]):
                nodes = [int(n) for n in edge.split("-")]
                element1 = self.elements[edgeConnectedElements[0]]
                element2 = self.elements[edgeConnectedElements[1]]
                faces_e1 = get_element_faces(element1.ica,ordered=True,toString =True)
                faces_e2 = get_element_faces(element2.ica,ordered=True,toString =True) 
                for n in nodes:
                    if not old_node_to_new.__contains__(n):
                        nodeNum = n
                        allConnectedElements = self.nodeToElements[nodeNum] 
                        connectedElements = []
                        for conn_element in allConnectedElements:
                            element = self.elements[conn_element]
                            add = True
                            for el_types in elementsNotIncluded:
                                if element.properties['mat'].count(el_types):
                                    add=False
                                    break
                            if add:
                                connectedElements.append(conn_element)
                        elementsConnected_To_new_node_map = {}
                        if len(connectedElements) <= 4:
                            if (replace != 0) or (len(elementsNotIncluded) != 0):
                                self.replace_element(element2.num,replace=replace)
                            else:
                                if not deleted_elements.count(element2.num):
                                    self.delete_element(element2.num)
                                    deleted_elements.append(element2.num)
                        # maxNodeNum += 1
                        # newNodeNum = maxNodeNum
                        # nodeCoords = list(self.nodes[nodeNum])
                        # node_idx = element2.ica.index(nodeNum)
                        # element2.ica[node_idx] = newNodeNum
                        # newNodes[newNodeNum] = nodeCoords
                        # old_node_to_new[nodeNum]=newNodeNum
                        # if elementsConnected_To_new_node_map.__contains__(newNodeNum):
                        #     elementsConnected_To_new_node_map[newNodeNum].append(element2.num)
                        # else:
                        #     elementsConnected_To_new_node_map[newNodeNum] = [element2.num]
                        # if len(connectedElements) > 2:
                        #     for e in connectedElements:
                        #         if e != element2.num and e != element1.num:
                        #             connectedEle = self.elements[e]
                        #             faces_e = get_element_faces(connectedEle.ica,ordered=True,toString =True)
                        #             shared_Faces = False
                        #             for f in faces_e:
                        #                 if faces_e2.count(f) or faces_e1.count(f):
                        #                     shared_Faces = True
                        #                     if faces_e2.count(f):
                        #                         node_idx = connectedEle.ica.index(nodeNum)
                        #                         connectedEle.ica[node_idx] = newNodeNum
                        #                         elementsConnected_To_new_node_map[newNodeNum].append(e)
                        #             if not shared_Faces:                                    
                        #                 ele = self.elements[e]
                        #                 maxNodeNum += 1
                        #                 newNodeNum2 = maxNodeNum
                        #                 nodeCoords = list(self.nodes[nodeNum])
                        #                 node_idx = ele.ica.index(nodeNum)
                        #                 ele.ica[node_idx] = newNodeNum2
                        #                 newNodes[newNodeNum2] = nodeCoords
                        #                 old_node_to_new[nodeNum] = newNodeNum2
                        #                 elementsConnected_To_new_node_map[newNodeNum2] = [ele.num]
                                    
                        # # Remove from node to element map
                        # for newNodeNum,elementsConnected_To_new_node in elementsConnected_To_new_node_map.items():
                        #     for e in elementsConnected_To_new_node:
                        #         self.nodeToElements[nodeNum].remove(e)
                        #     self.nodeToElements[newNodeNum] = elementsConnected_To_new_node
                    
            # for key,node in newNodes.items():
            #     self.nodes[key] = node;

            
    def smooth_mesh(self, coeffs, iterations, elementsNotIncluded = []):
        print("Starting mesh smoothing")
        elementMap = self.create_elements_map(elementsNotIncluded=elementsNotIncluded)
        if len(elementMap)>0:
            from element_functions import create_node_to_elem_map
            self.locate_boundary_faces(elementsNotIncluded=elementsNotIncluded);
            # self.clean_mesh(elementsNotIncluded=elementsNotIncluded)
            surfaceNodeConnectivity = smooth.create_surface_connectivity(self.boundary_element_map)  
            elementMapFull = self.create_elements_map()
            nodeToElemMap = create_node_to_elem_map(elementMapFull)
            for iteration in range(iterations):
                smooth.perform_smoothing(iteration, coeffs, surfaceNodeConnectivity, self.nodes, elementMapFull, nodeToElemMap=nodeToElemMap)
        else:
            print("No elements selected to smooth")
    
    def create_elements_map(self, elementsNotIncluded = []):
        elementMap = {}
        for elementNo, element in self.elements.items():
            add = True
            for el_types in elementsNotIncluded:
                if element.properties['mat'].count(el_types):
                    add=False
                    break;
            if add:
                elementMap[elementNo] = element.ica
        return elementMap
    
    def create_node_to_element_connectivity(self):
        self.nodeToElements = {}
        for element in self.elements.values():
            for node in element.ica:
                connectedElements = []
                if self.nodeToElements.__contains__(node):
                    connectedElements = self.nodeToElements[node]
                else:
                    self.nodeToElements[node] = connectedElements
                connectedElements.append(element.num)
                
    def create_edge_to_element_connectivity(self, elementsNotIncluded= []):
        from element_functions import get_edges
        from Smoothing import get_element_faces
        edgesToElements_tmp = {}
        joined_edges_and_no_faces_count = 0
        joined_edges_and_no_faces = []
        for element in self.elements.values():
            add = True
            for el_types in elementsNotIncluded:
                if element.properties['mat'].count(el_types):
                    add=False
                    break
            if add:
                edges = get_edges(element.ica)
                for edge in edges:
                    connectedElements = []
                    if edgesToElements_tmp.__contains__(edge):
                        connectedElements = edgesToElements_tmp[edge]
                    else:
                        edgesToElements_tmp[edge] = connectedElements
                    connectedElements.append(element.num)
                    
        edgesToElements = {}
        for edge, elements in edgesToElements_tmp.items():
            if len(elements) == 2:
                element1 = elements[0]
                element2 = elements[1]
                faces1 = get_element_faces(self.elements[element1].ica, ordered = True, toString = True)
                faces2 = get_element_faces(self.elements[element2].ica, ordered = True, toString = True)
                shared_face = False
                for face in faces1:                    
                    if faces2.count(face):
                        shared_face = True
                        break
                if not shared_face:
                    joined_edges_and_no_faces_count += 1
                    edgesToElements[edge]= list(elements)
        return edgesToElements
                

    
    def calculate_node_coords(self,elementX,elementY,elementZ,i,size):
        coordx = int((i-1)/((elementZ+1)*(elementY+1)))
        tmp = i - (coordx*((elementZ+1)*(elementY+1)))
        coordy = int((tmp-1)/(elementZ+1))
        coordz = (tmp - (coordy*(elementZ+1)))-1
        return [float(d) for d in [coordx*size, coordy*size, coordz*size]]
    
    def write_to_file(self, path, filename, labels_map, filetype="abaqus", boundary=True):
        print("Writing mesh data to file in "+ filetype + " format")
        if (path[-1] != "\\"):
            path += "\\" 
        elementMap = self.create_elements_map()
        material_mapping = labels_map.create_material_sets(self.elements,file_format=filetype)
        if (filetype.lower() == "ucd"):
            if boundary:
                if not hasattr(self, "boundary_element_map"):
                    self.locate_boundary_faces();  
            
            homogenized_labels_map = labels_map.get_homogenized_labels_map();
            rw.writeUCD(path, filename, self.nodes, elementMap, boundaryElementMap=self.boundary_element_map, 
                        elementToElsetMap=material_mapping, elset_number_Mappings=homogenized_labels_map)
            
        elif (filetype.lower() == "vtk"):
            elementToMaterial = {}
            rw.writeVTK(path, filename, self.nodes, elementMap, elementToMaterial=material_mapping)              
              
        else:
            rw.writeABQ(path, filename, self.nodes, elementMap, elsetsMap=material_mapping)
            

    
                   
        
         
        