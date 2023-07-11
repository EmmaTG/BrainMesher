# -*- coding: utf-8 -*-
"""
Created on Thu May 25 15:39:38 2023

@author: grife
"""
import ABQ_UCD_handling as rw
import Smoothing as smooth
from abc import ABC, abstractmethod

class IElement(ABC): 
    @abstractmethod
    def get_faces(self, stringyfy, order):
       pass
    
    @abstractmethod
    def get_edges(self, stringyfy, order):
        pass
    
class Element():
    def __init__(self, number, ica, **kwargs):        
        self.ica = ica
        self.num = number
        self.properties = {}
        for key,value in kwargs.items():
            self.properties[key] = value;
    
    def setMaterial(self,mat):
        try:
            list(mat)
        except TypeError:
            mat = [mat]
        self.properties['mat'] = list(mat)    
    
    
    def getMaterial(self):
        return self.properties['mat']
     
            
    def get_nodes_involved(self, faces, stringyfy=True, order=True):
       
        node_faces = []
        for f in faces:
            new_face = []
            for n in f:
                new_face.append(self.ica[n-1])                
            if order:
                new_face = sorted(new_face)
            if stringyfy:
                node_faces.append("-".join(str(x) for x in new_face))
            else:
                node_faces.append(new_face)
        return node_faces  
    
    def calculate_element_centroid(self,nodeMap):
        """
        Calculates element centroid

        Parameters
        ----------
        nodeMap : Map(int, array)
            Map of node numbers (keys) to coordinates (values).

        Returns
        -------
        centroid : float
            Element centroid.

        """
        centroid = [0,0,0]
        for n in self.ica:
            coords = nodeMap[n]
            for i in range(3):
                centroid[i] += coords[i]
                
        for i in range(3):
            centroid[i] /= len(self.ica)
        return centroid 
    
class HexElement(Element, IElement):
    
    def __init__(self, number, ica, **kwargs): 
        Element.__init__(self, number, ica, **kwargs);
    
    def get_faces(self, stringyfy=True, order=True):
        face_ABQ = [[1,2,3,4],
                [5,8,7,6],
                [1,5,6,2],
                [2,6,7,3],
                [3,7,8,4],
                [4,8,5,1]]
        return super().get_nodes_involved(face_ABQ, stringyfy=stringyfy, order=order)
    
    def get_edges(self, stringyfy=True, order=True):
        edge_classification = [[0,1],[1,2],[2,3],[3,0],
                               [4,5],[5,6],[6,7],[7,4],
                               [0,4],[1,5],[2,6],[3,7]]
        return super().get_nodes_involved(edge_classification, stringyfy=stringyfy, order=order)


class QuadElement(Element, IElement):
    
    def __init__(self, number, ica, **kwargs): 
        Element.__init__(self, number, ica, **kwargs);

    def get_faces(self, stringyfy=True, order=True):
        face_ABQ = [[1,2,3,4]] 
        return super().get_nodes_involved(face_ABQ, stringyfy=stringyfy, order=order)
    
    def get_edges(self, stringyfy=True, order=True):
        edge_classification = [[0,1],[1,2],[2,3],[3,0]]
        return super().get_nodes_involved(edge_classification, stringyfy=stringyfy, order=order)
    
    
class Mesh():
    
    def __init__(self, pointData,voxel_size):
        self.elements = {}
        self.nodes = {}
        self.elementToPointCloud = {}
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
                if not self.nodes.get(i,False):
                    coords = self.calculate_node_coords(elementX,elementY,elementZ,i,voxel_size)
                    self.nodes[i] = coords
                if not self.nodes.get(newNode,False):
                    coords = self.calculate_node_coords(elementX,elementY,elementZ,newNode,voxel_size)
                    self.nodes[newNode] = coords 
            element_ica += element_ica_tmp
            element = HexElement(elementNo, element_ica, mat=[m])
            self.elements[int(elementNo)] = element
            self.elementToPointCloud[int(elementNo)] = [x,y,z,m]
    
    
    def clean_mesh(self, elementsNotIncluded = [], replace=0):
        iteration = 0
        count = 1
        total_count = 0
        self.create_node_to_element_connectivity();
        while ((count > 0) and (iteration<10)):
            count = 0;
            iteration += 1
            count += self.clean_mesh_edges(elementsNotIncluded = elementsNotIncluded, replace=replace);
            count += self.clean_mesh_nodes(elementsNotIncluded = elementsNotIncluded, replace=replace);
            if ((replace != 0) and (iteration == 1)):
                elementsNotIncluded.append(replace);
            total_count += count;
        print(str(total_count) + " elements deleted/replaced due to poor node/edge connectivity in " + str(iteration) + " iterations")
         
    
    def locate_boundary_element_map(self, elementsNotIncluded = []):
        print("Identifying boundary faces")
        elementMap = self.create_elements_map(elementsNotIncluded=elementsNotIncluded)
        return smooth.get_boundary_element_map(elementMap)
    
    def locate_elements_on_boundary(self, elementsNotIncluded = []):
        print("Identifying boundary elements")
        elementMap = self.create_elements_map(elementsNotIncluded=elementsNotIncluded)
        return smooth.get_elements_on_boundary(elementMap)
        
    def replace_outer_region(self, white_matter_label, replace_label, elementsNotIncluded=[]):
        print("Cleaning brain boundary")
        elements_on_boundary = self.locate_elements_on_boundary(elementsNotIncluded=elementsNotIncluded)
        for elem in elements_on_boundary:
            element = self.elements[elem]
            materials = element.properties['mat']
            if materials.count(white_matter_label):
                materials.remove(white_matter_label)
                materials.insert(0,replace_label)
        
        
    def clean_mesh_nodes(self, elementsNotIncluded = [], replace=0):
        print("Cleaning mesh node")
        cleaned_elements = []
        cleaned_nodes= []
        node_keys = list(self.nodes.keys())
        for key in node_keys:
            if not cleaned_nodes.count(key):
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
                    assert numberSharedNodes>0;
                    if numberSharedNodes == 1:
                        if not cleaned_elements.count(element2.num):
                            cleaned_elements.append(element2.num)
                            cleaned_nodes = cleaned_nodes + element2.ica
                            if (replace != 0) or (len(elementsNotIncluded) != 0):
                                self.replace_element(element2.num,replace=replace)
                            else:
                                self.delete_element(element2.num) 
        # print(str(len(cleaned_elements)) + " element deleted due to poor node connectivity")
        return len(cleaned_elements)
    
    def delete_element(self, element_number):
        if self.elements.get(element_number,False):    
            element = self.elements[element_number]
            ica = element.ica
            for n in ica:
                connectedElements = self.nodeToElements[n]
                connectedElements.remove(element_number)
                if len(connectedElements) == 0:
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
        edgesToElementsMap = self.create_edge_to_element_connectivity(elementsNotIncluded)
        edgesToElements = self.get_edge_without_shared_face(edgesToElementsMap)
        old_node_to_new = {}
        cleaned_elements = []
        for edge, edgeConnectedElements in edgesToElements.items():
            if not cleaned_elements.count(edgeConnectedElements[0]) and not cleaned_elements.count(edgeConnectedElements[1]):
                nodes = [int(n) for n in edge.split("-")]
                element1 = self.elements[edgeConnectedElements[0]]
                element2 = self.elements[edgeConnectedElements[1]] 
                for n in nodes:
                    if not old_node_to_new.get(n,False):
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
                        if len(connectedElements) <= 4:
                            if not cleaned_elements.count(element2.num):
                                cleaned_elements.append(element2.num)
                                if (replace != 0) or (len(elementsNotIncluded) != 0):
                                    self.replace_element(element2.num,replace=replace)
                                else:                                
                                    self.delete_element(element2.num)
        
        # print(str(len(cleaned_elements)) + " element deleted due to poor edge connectivity")
        return len(cleaned_elements)                 
            
    def smooth_mesh(self, coeffs, iterations, elementsNotIncluded = []):
        print("Starting mesh smoothing")
        boundary_element_map = self.locate_boundary_element_map(elementsNotIncluded)
        if len(boundary_element_map)>0:
            node_to_boundary_element_map = self.create_node_to_elem_map(boundary_element_map)
            surfaceNodeConnectivity = smooth.create_surface_connectivity(boundary_element_map,node_to_boundary_element_map)  
            elementMapFull = self.create_elements_map()
            nodeToElemMap = self.create_node_to_elem_map(elementMapFull)
            for iteration in range(iterations):
                smooth.perform_smoothing(iteration, coeffs, surfaceNodeConnectivity, self.nodes, elementMapFull, nodeToElemMap=nodeToElemMap)
        else:
            print("No elements selected to smooth")
            
    def create_node_to_elem_map(elementMap):
        "Creating node to element connectivity"
        nodeToElemMap = {}
        for e,ica in elementMap.items():        
            for node in ica:
                if nodeToElemMap.get(node,False):
                    elements = nodeToElemMap[node]
                else:
                    elements = []
                elements.append(e)
                nodeToElemMap[node] = elements  
        return nodeToElemMap
    
    def create_elements_map(self, elements = [], elementsNotIncluded = [], elementsIncluded = []):
        elementMap = {}
        if len(elements) == 0:
            elements = self.elements
        for elementNo, element in elements.items():
            add = True
            for el_types in elementsNotIncluded:
                if element.properties['mat'].count(el_types):
                    add=False
                    break;
            if add:                
                if len(elementsIncluded)>0:
                    for el_types in elementsIncluded:
                        if element.properties['mat'].count(el_types):
                            elementMap[elementNo] = element.ica
                else:
                    elementMap[elementNo] = element.ica
        return elementMap
    
    def create_node_to_element_connectivity(self):
        self.nodeToElements = {}
        for element in self.elements.values():
            for node in element.ica:
                connectedElements = []
                if self.nodeToElements.get(node,False):
                    connectedElements = self.nodeToElements[node]
                else:
                    self.nodeToElements[node] = connectedElements
                connectedElements.append(element.num)
                
    def create_edge_to_element_connectivity(self, elementsNotIncluded= [], elementsMap={}):
        edgesToElements_tmp = {}
        if len(elementsMap)>0:
            elementValues = elementsMap.values()
        else:
            elementValues = self.elements.values()
        for element in elementValues:
            add = True
            for el_types in elementsNotIncluded:
                if element.properties['mat'].count(el_types):
                    add=False
                    break
            if add:
                edges = element.get_edges(stringyfy=True, order=True)
                for edge in edges:
                    connectedElements = []
                    if edgesToElements_tmp.get(edge,False):
                        connectedElements = edgesToElements_tmp[edge]
                    else:
                        edgesToElements_tmp[edge] = connectedElements
                    connectedElements.append(element.num)
        return edgesToElements_tmp
                           
                
    def get_edge_without_shared_face(self,edgesToElements_map):
        edgesToElements = {}
        for edge, elements in edgesToElements_map.items():
            if len(elements) == 2:
                element1 = elements[0]
                element2 = elements[1]
                faces1 = self.elements[element1].get_faces(order = True, stringyfy = True)
                faces2 = self.elements[element2].get_faces(order = True, stringyfy = True)
                shared_face = False
                for face in faces1:                    
                    if faces2.count(face):
                        shared_face = True
                        break
                if not shared_face:
                    edgesToElements[edge]= list(elements)
        return edgesToElements
    
    def calculate_node_coords(self,elementX,elementY,elementZ,i,size):
        coordx = int((i-1)/((elementZ+1)*(elementY+1)))
        tmp = i - (coordx*((elementZ+1)*(elementY+1)))
        coordy = int((tmp-1)/(elementZ+1))
        coordz = (tmp - (coordy*(elementZ+1)))-1
        return [float(d) for d in [coordx*size, coordy*size, coordz*size]]
    
    def write_to_file(self, path, filename, labels_map, filetype="abaqus", boundaryElementMap={}):
        print("Writing mesh data to file in "+ filetype + " format")
        if (path[-1] != "\\"):
            path += "\\" 
        # Map of element number (key) to element ica (value)
        element_ica_Map = self.create_elements_map()
        # for VTK and UCD files: Map of element number (key) to list of material numbers (value)
        # for Abaqus files: Map of material names (key) to elements numbers
        element_to_material_mapping = labels_map.create_material_sets(self.elements,file_format=filetype) 
        boundary_element_ica_map = {}
        boundary_material_mapping = {}
        if len(boundaryElementMap)>0:
            # Map of boundary element number [quad] (key) to element [quad] ica (value)
            boundary_element_ica_map = self.create_elements_map(boundaryElementMap)
            # for VTK and UCD files: Map of boundary element number (key) to list of boundary numbers (value)
            # for Abaqus files: Map of boundary names (key) to boundary elements numbers
            boundary_material_mapping = labels_map.create_material_sets(boundaryElementMap,file_format=filetype)
        if (filetype.lower() == "ucd"):
            rw.writeUCD(path, filename, self.nodes, 
                        element_ica_Map, elementToMaterialMap=element_to_material_mapping, 
                        boundaryElementMap=boundary_element_ica_map, boundaryElementToMaterial=boundary_material_mapping)
            
        elif (filetype.lower() == "vtk"):
            rw.writeVTK(path, filename, self.nodes, element_ica_Map, elementToMaterial=element_to_material_mapping, 
                        boundaryElementMap=boundary_element_ica_map, boundaryElementToMaterial=boundary_material_mapping)              
              
        else:
            boundary_nodes = {}
            boundary_nodes_list = []
            for ica in boundary_element_ica_map.values():
                for n in ica:
                    if not boundary_nodes_list.count(n):
                        boundary_nodes_list.append(n)            
            boundary_nodes['Boundary nodes'] = boundary_nodes_list
            rw.writeABQ(path, filename, self.nodes, element_ica_Map, elsetsMap=element_to_material_mapping, nsetMap=boundary_nodes)
            

    
                   
        
         
        