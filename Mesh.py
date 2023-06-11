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
        self.locate_boundary_faces(elementsNotIncluded = [24.0])
        self.remove_outer_white_matter(2, 3)
        
    def locate_boundary_faces(self, elementsNotIncluded = []):
        print("Identifying boundary faces")
        elementMap = self.create_elements_map(elementsNotIncluded=elementsNotIncluded)
        self.boundary_element_map, self.elements_on_boundary, self.volume_elem_to_boundary = smooth.get_boundary_surfaces(elementMap)
        
    def remove_outer_white_matter(self, white_matter_label, replace_label):
        print("Cleaning brain boundary")
        for elem in self.elements_on_boundary:
            element = self.elements[elem]
            materials = element.properties['mat']
            if materials.count(white_matter_label):
                materials.remove(white_matter_label)
                materials.insert(0,replace_label)
        
        
    def clean_mesh(self, elementsNotIncluded = []):
        print("Cleaning mesh")
        self.create_node_to_element_connectivity();
        maxNodeNum = max(self.nodes.keys())
        newNodes = {}
        for key in self.nodes.keys():
            if not hasattr(self,"nodeToElements"):
                self.create_node_to_element_connectivity()
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
                    # print("Cleaning mesh at node {} between elements {} and {}".format(key,element1.num,element2.num))                    
                    maxNodeNum += 1
                    nodeCoords = list(self.nodes[key])
                    node_idx = element2.ica.index(key)
                    element2.ica[node_idx] = maxNodeNum
                    newNodes[maxNodeNum] = nodeCoords
                    # Remove from node to element map
                    self.nodeToElements[key].remove(connectedElements[1])
                    self.nodeToElements[maxNodeNum] = element2.num
        for key,node in newNodes.items():
            self.nodes[key] = node;
                
    def smooth_mesh(self, coeffs, iterations, elementsNotIncluded = []):
        print("Starting mesh smoothing")
        elementMap = self.create_elements_map(elementsNotIncluded=elementsNotIncluded)
        if len(elementMap)>0:
            self.locate_boundary_faces(elementsNotIncluded=elementsNotIncluded);
            self.clean_mesh(elementsNotIncluded=elementsNotIncluded)
            surfaceNodeConnectivity = smooth.create_surface_connectivity(self.boundary_element_map)
            for iteration in range(iterations):
                smooth.perform_smoothing(iteration, coeffs, surfaceNodeConnectivity, self.nodes, elementMap)
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
                
    def create_edge_to_element_connectivity(self):
        from element_functions import get_edges
        self.edgesToElements = {}
        joined_edges = 0
        for element in self.elements.values():
            if not element.properties['mat'].count(24):
                edges = get_edges(element.ica)
                for edge in edges:
                    connectedElements = []
                    if self.edgesToElements.__contains__(edge):
                        connectedElements = self.edgesToElements[edge]
                    else:
                        self.edgesToElements[edge] = connectedElements
                    connectedElements.append(element.num)
        for edge, elements in self.edgesToElements.items():
            if len(elements) == 2:
                joined_edges += 1
            
        print (joined_edges)
                

    
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
            for e_num,e in self.elements.items():
                materials = e.properties['mat'] 
                elementToMaterial[e_num] = int(materials[0])
            rw.writeVTK(path, filename, self.nodes, elementMap, elementToMaterial=elementToMaterial)              
              
        else:
            rw.writeABQ(path, filename, self.nodes, elementMap, elsetsMap=material_mapping)
            

    
                   
        
         
        