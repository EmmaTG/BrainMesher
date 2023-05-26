# -*- coding: utf-8 -*-
"""
Created on Thu May 25 15:39:38 2023

@author: grife
"""
import ABQ_UCD_handling as rw
import Smoothing as smooth

class Element():
    
    def __init__(self, number, ica, **kwargs):        
        self.ica = ica
        self.num = number
        for key,value in kwargs.items():
            self.key = value;
    
    
class Mesh():
    
    def __init__(self, pointData,voxel_size):
        self.elements = {}
        self.nodes = {}
        self.mat_sets = {}
        for m in set(pointData[:,3]):
            self.mat_sets[m] = []
        self.create_mesh(pointData,voxel_size)
        
    
    def create_mesh(self, pointData, voxel_size):
        print("Creating mesh from point cloud data")
        [minX,minY,minZ] = pointData[:,:3].min(axis=0)
        [maxX,maxY,maxZ] = pointData[:,:3].max(axis=0)
        
        elementX = maxX-minX+1
        elementY = maxY-minY+1
        elementZ = maxZ-minZ+1
        
        elementNo = 0;
        for p in pointData:
            elementNo += 1
            [x,y,z,m] = p
            self.mat_sets[m].append(elementNo)
            startNode = (z+1) + (elementX+1)*y + ((elementZ+1)*(elementY+1))*x
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
            element = Element(elementNo, element_ica, mat=m)
            self.elements[int(elementNo)] = element
        self.clean_mesh();
        
    def locate_boundary_faces(self):
        print("Identifying boundary faces")
        elementMap = {}
        for elementNo, element in self.elements.items():
            elementMap[elementNo] = element.ica
        self.boundary_element_map, self.elements_on_boundary, self.volume_elem_to_boundary = smooth.get_boundary_surfaces(elementMap)  
        
        
    def clean_mesh(self):
        print("Cleaning mesh")
        self.create_node_to_element_connectivity();
        maxNodeNum = max(self.nodes.keys())
        newNodes = {}
        for key in self.nodes.keys():
            if not hasattr(self,"nodeToElements"):
                self.create_node_to_element_connectivity()
            connectedElements = self.nodeToElements[key]
            if len(connectedElements) == 2:
                element1 = self.elements[connectedElements[0]]
                element2 = self.elements[connectedElements[1]]
                numberSharedNodes = 0
                for node1 in element1.ica:
                    if element2.ica.count(node1)>0:
                        numberSharedNodes += 1
                assert numberSharedNodes>0
                if numberSharedNodes == 1:
                    print("Cleaning mesh at node {} between elements {} and {}".format(key,element1.num,element2.num))                    
                    maxNodeNum += 1
                    nodeCoords = list(self.nodes[key])
                    node_idx = element2.ica.index(key)
                    element2.ica[node_idx] = maxNodeNum
                    newNodes[maxNodeNum] = nodeCoords
        for key,node in newNodes.items():
            self.nodes[key] = node;
                
    def smooth_mesh(self, coeffs, iterations):
        print("Starting mesh smoothing")
        if not hasattr(self, "boundary_element_map"):
            self.locate_boundary_faces();
        elementMap= {}
        for elementNo, element in self.elements.items():
            elementMap[elementNo] = element.ica
        surfaceNodeConnectivity = smooth.create_surface_connectivity(self.boundary_element_map)
        for iteration in range(iterations):
            smooth.perform_smoothing(iteration, coeffs, surfaceNodeConnectivity, self.nodes, elementMap)
        
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
        
    def calculate_node_coords(self,elementX,elementY,elementZ,i,size):
        coordx = int((i-1)/((elementZ+1)*(elementY+1)))
        tmp = i - (coordx*((elementZ+1)*(elementY+1)))
        coordy = int((tmp-1)/(elementZ+1))
        coordz = (tmp - (coordy*(elementZ+1)))-1
        return [float(d) for d in [coordx*size, coordy*size, coordz*size]]
    
    def write_to_file(self, path, filename, filetype="abaqus", boundary=True):
        print("Writing mesh data to file in "+ filetype + " format")
        if (path[-1] != "\\"):
            path += "\\" 
        elementMap = {}
        for elementNo, element in self.elements.items():
            elementMap[elementNo] = element.ica
        if (filetype.lower() == "ucd"):
            if boundary:
                if not hasattr(self, "boundary_element_map"):
                    self.locate_boundary_faces();                
            rw.writeUCD(path, filename, self.nodes, elementMap, boundaryElementMap=self.boundary_element_map)
        else:
            rw.writeABQ(path, filename, self.nodes, elementMap)
            

    
                   
        
         
        