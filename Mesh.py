# -*- coding: utf-8 -*-
"""
Created on Thu May 25 15:39:38 2023

@author: grife
"""
# import ABQ_UCD_handling as rw
import Smoothing as smooth
from abc import ABC, abstractmethod

class INode():
    def __init__(self, number, coords):        
        self.number = number
        self.coords = coords
        self.data = {}
    
    def addData(self,name,value):
        self.data[name] = value;
    
    
    def setCoords(self, coords):
        self.coords = coords;
    
    def getCoords(self):
        return self.coords;
        
class Node(INode):    
    pass;

class IElement(ABC): 
    @abstractmethod
    def get_faces(self, stringyfy, order):
       pass
    
    @abstractmethod
    def get_edges(self, stringyfy, order):
        pass
    
class ElementCalculations():
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
            coords = nodeMap[n].getCoords()
            for i in range(3):
                centroid[i] += coords[i]
                
        for i in range(3):
            centroid[i] /= len(self.ica)
        return centroid ;
    
class Element(ElementCalculations):
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
    
class MeshTransformations():
    @staticmethod
    def rotate_mesh(nodeMap, axis=0, degrees=90):
        from math import pi,cos,sin
        import numpy as np
        deg = degrees*pi/180
        if (axis == 0):
            R = np.array([[1,0,0],[0, cos(deg), -1*sin(deg)],[0, sin(deg), cos(deg)]])      
        if (axis == 1):
            R = np.array([[cos(deg),0,sin(deg)],[0, 1, 0],[-1*sin(deg), 0, cos(deg)]])
        if (axis == 2):
            R = np.array([[cos(deg), -1*sin(deg), 0],[sin(deg), cos(deg), 0],[0, 0, 1]])
        for n in nodeMap.values():
            coords = n.getCoords();
            newCoords = np.matmul(np.array(coords),R)
            coords[0] = round(newCoords[0],3)
            coords[1] = round(newCoords[1],3)
            coords[2] = round(newCoords[2],3)
    
    @staticmethod            
    def scale_mesh(nodeMap,scale=[1,1,1],reduce=True):
        if reduce:
            scale[0] = 1/scale[0]
            scale[1] = 1/scale[1]
            scale[2] = 1/scale[2]        
        for n in nodeMap.values():
            coords = n.getCoords();
            coords[0] = coords[0]*scale[0]
            coords[1] = coords[1]*scale[1]
            coords[2] = coords[2]*scale[2]

    @staticmethod
    def translate_mesh(nodeMap,distance=[1,1,1]):        
        for n in nodeMap.values():
            coords = n.getCoords();
            coords[0] = coords[0]+distance[0]
            coords[1] = coords[1]+distance[1]
            coords[2] = coords[2]+distance[2]
    
class MeshUtils():

    def create_elements_map(self, elements, elementsNotIncluded, elementsIncluded = []):
        elementMap = {}
        for elementNo, element in elements.items():
            add = True
            for el_types in elementsNotIncluded:
                if element.getMaterial().count(el_types):
                    add=False
                    break;
            if add:                
                if len(elementsIncluded)>0:
                    for el_types in elementsIncluded:
                        if element.getMaterial().count(el_types):
                            elementMap[elementNo] = element
                else:
                    elementMap[elementNo] = element
        return elementMap
    
    def create_elements_ica_map(self, elements):
        elementMap = {}
        for elementNo, element in elements.items():
            elementMap[elementNo] = element.ica
        return elementMap

    def locate_boundary_element_map(self, elements, elementsNotIncluded = []):
        print("Locating boundary elements")
        elementMap = self.create_elements_map(elements, elementsNotIncluded)
        face_to_elems_map = {}
        surface_face_to_elems_map = {}
        for e,element in elementMap.items():
            list_of_faces = elementMap[e].get_faces(True,True)
            for face_key in list_of_faces:                                             # Create map key 
                if face_to_elems_map.get(face_key,False):                            # Check if face key already in map
                   face_to_elems_map[face_key].append(e)                    # key already in face so append element to array (NOT surface face)
                   if surface_face_to_elems_map.get(face_key,False):                   # If previously classified as a free surface; remove from this map
                       del surface_face_to_elems_map[face_key]
                else:
                    face_to_elems_map[face_key] = [e]                                   # If not in map, add to map
                    surface_face_to_elems_map[face_key] = e
            
        boundary_element_map = {}
        for face_key,e in surface_face_to_elems_map.items():    
            faces = elementMap[e].get_faces(True,True)
            for face_num,f in enumerate(faces):
                if f == face_key:
                    compound_key = "-".join([str(e),str(face_num)])
                    boundary_element_map[compound_key] = elementMap[e].get_faces(False,False)[face_num]
                    break    
        
        return boundary_element_map
    
        
    def locate_elements_on_boundary(self, elements, elementsNotIncluded = []):
        print("Locating elements on the boundary")
        elementMap = self.create_elements_map(elements, elementsNotIncluded)
        face_to_elems_map = {}
        surface_face_to_elems_map = {}
        for e,element in elementMap.items():       
            list_of_faces = elementMap[e].get_faces(True,True)
            for face_key in list_of_faces:                                             # Create map key 
                if face_to_elems_map.get(face_key,False):                            # Check if face key already in map
                   connected_elements =  face_to_elems_map[face_key]                    # key already in face so append element to array (NOT surface face)
                   connected_elements.append(e)
                   if surface_face_to_elems_map.get(face_key,False):                   # If previously classified as a free surface; remove from this map
                       del surface_face_to_elems_map[face_key]
                else:
                    face_to_elems_map[face_key] = [e]                                   # If not in map, add to map
                    surface_face_to_elems_map[face_key] = e
            
        elements_on_boundary = []
        for face_key,e in surface_face_to_elems_map.items():    
            if not elements_on_boundary.count(e):
                    elements_on_boundary.append(e) 
        
        return elements_on_boundary
    
    def create_node_to_elem_map(self,elementICAMap):
        "Creating node to element connectivity"
        nodeToElemMap = {}
        for e,ica in elementICAMap.items():        
            for node in ica:
                if nodeToElemMap.get(node,False):
                    elements = nodeToElemMap[node]
                else:
                    elements = []
                elements.append(e)
                nodeToElemMap[node] = elements  
        return nodeToElemMap
    
    def create_surface_connectivity(self, boundary_element_map, nodeToBoundaryElementMap):
        ## Create surface connectivity map
        print("Creating node surface connectivty map")
        surfaceNodeConnectivity = {}
        for node,compoundKeys in nodeToBoundaryElementMap.items():
            connectedNodes = []
            for f in compoundKeys:
                faceICA = boundary_element_map[f] 
                idx = faceICA.index(node)
                idx1 = idx + 1 if idx < 3 else 0
                idx2 = idx -1 if idx > 0 else 3
                connectedNodes.append(faceICA[idx1])
                connectedNodes.append(faceICA[idx2])
            surfaceNodeConnectivity[node] = list(set(connectedNodes))
        return surfaceNodeConnectivity
    
    def create_edge_to_element_connectivity(self, elementsMap, elementsNotIncluded= []):
        edgesToElements_tmp = {}
        for element in elementsMap.values():
            add = True
            for el_types in elementsNotIncluded:
                if element.getMaterial().count(el_types):
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
    
    
  
class Mesh():
    
    def __init__(self):
        self.elements = {}
        self.nodes = {}
        self.boundaryElements = {}
        self.elementToPointCloud = {}
        self.dataToWrite = {}
        self._meshUtils = MeshUtils();
        
    def addBoundaryElements(self,boundaryElementsMap):
        self.boundaryElements = self.boundaryElements | boundaryElementsMap;
    
    def getBoundingBox(self):
        maxV = [-1000,-1000,-1000]
        minV = [1000,1000,1000]        
        for node in self.nodeMap.values():
            n = node.getCoords();
            for d in range(3):
                if maxV[d]<n[d]:
                    maxV[d] = n[d]
                if minV[d]>n[d]:
                    minV[d] = n[d]
        return maxV + minV
    
    def locate_boundary_element_map(self,elementsNotIncluded=[]):
        return self._meshUtils.locate_boundary_element_map(self.elements,elementsNotIncluded = elementsNotIncluded)
        
    def remove_region(self,region_value):
            element_keys = list(self.elements.keys())    
            for element_num in element_keys:
                e = self.elements[element_num]
                if e.getMaterial().count(region_value):
                    self.delete_element(element_num)
                    
    def create_mesh_from_Point_Cloud(self, pointData, voxel_size):
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
                    self.nodes[i] = Node(i,coords)
                if not self.nodes.get(newNode,False):
                    coords = self.calculate_node_coords(elementX,elementY,elementZ,newNode,voxel_size)
                    self.nodes[newNode] = Node(newNode,coords) 
            element_ica += element_ica_tmp
            element = HexElement(elementNo, element_ica, mat=[m])
            self.elements[int(elementNo)] = element
            self.elementToPointCloud[int(elementNo)] = [x,y,z,m]
        self.create_node_to_element_connectivity();
    
    
    def clean_mesh(self, elementsNotIncluded = [], replace=0):
        iteration = 0
        count = 1
        total_count = 0
        while ((count > 0) and (iteration<10)):
            count = 0;
            iteration += 1
            count += self.clean_mesh_edges(elementsNotIncluded = elementsNotIncluded, replace=replace);
            count += self.clean_mesh_nodes(elementsNotIncluded = elementsNotIncluded, replace=replace);
            if ((replace != 0) and (iteration == 1)):
                elementsNotIncluded.append(replace);
            total_count += count;
        print(str(total_count) + " elements deleted/replaced due to poor node/edge connectivity in " + str(iteration) + " iterations")
        
    def replace_outer_region(self, white_matter_label, replace_label, elements_on_boundary):
        print("Cleaning brain boundary")
        for elem in elements_on_boundary:
            element = self.elements[elem]
            materials = element.getMaterial()
            if materials.count(white_matter_label):
                materials.remove(white_matter_label)
                materials.insert(0,replace_label)
        
        
    def clean_mesh_nodes(self, elementsNotIncluded = [], replace=0):
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
                        if element.getMaterial().count(el_types):
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
        element.setMaterial(replace)            
    
    def clean_mesh_edges(self, elementsNotIncluded = [], replace=0):
        edgesToElementsMap = self._meshUtils.create_edge_to_element_connectivity(self.elements, elementsNotIncluded)
        edgesToElements = self.get_edge_without_shared_face(edgesToElementsMap)
        old_node_to_new = {}
        cleaned_elements = []
        for edge, edgeConnectedElements in edgesToElements.items():
            if not cleaned_elements.count(edgeConnectedElements[0]) and not cleaned_elements.count(edgeConnectedElements[1]):
                nodes = [int(n) for n in edge.split("-")]
                element1 = self.elements[edgeConnectedElements[0]]
                element2 = self.elements[edgeConnectedElements[1]] 
                for n in nodes:
                    nodeNum = n
                    if not old_node_to_new.get(nodeNum,False):                        
                        allConnectedElements = self.nodeToElements[nodeNum] 
                        connectedElements = []
                        for conn_element in allConnectedElements:
                            element = self.elements[conn_element]
                            add = True
                            for el_types in elementsNotIncluded:
                                if element.getMaterial().count(el_types):
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
        return len(cleaned_elements)                 
            
    def smooth_mesh(self, coeffs, iterations, boundary_element_map):
        print("Starting mesh smoothing")
        if len(boundary_element_map)>0:
            node_to_boundary_element_map = self._meshUtils.create_node_to_elem_map(boundary_element_map)
            surfaceNodeConnectivity = self._meshUtils.create_surface_connectivity(boundary_element_map,node_to_boundary_element_map)
            elementICAMap = self._meshUtils.create_elements_ica_map(self.elements)
            nodeToElemMap = self._meshUtils.create_node_to_elem_map(elementICAMap)
            for iteration in range(iterations):
                smooth.perform_smoothing(iteration, coeffs, surfaceNodeConnectivity, self.nodes, elementICAMap, nodeToElemMap=nodeToElemMap)
        else:
            print("No elements selected to smooth")
    
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
    

        
         
        