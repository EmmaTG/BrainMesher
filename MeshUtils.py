# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 12:22:25 2023

@author: grife

Module of mesh utility functions
"""

def create_elements_map(elements, elementsNotIncluded, elementsIncluded = []):
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
    
def create_elements_ica_map(elements):
    elementMap = {}
    for elementNo, element in elements.items():
        elementMap[elementNo] = element.ica
    return elementMap

def locate_boundary_element_map(elements, elementsNotIncluded = []):
    print("Locating boundary elements")
    elementMap = create_elements_map(elements, elementsNotIncluded)
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

    
def locate_elements_on_boundary(elements, elementsNotIncluded = []):
    print("Locating elements on the boundary")
    elementMap = create_elements_map(elements, elementsNotIncluded)
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

def create_node_to_elem_map(elementICAMap):
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

def create_surface_connectivity(boundary_element_map, nodeToBoundaryElementMap):
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

def create_edge_to_element_connectivity(elementsMap, elementsNotIncluded= []):
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
    
    
  

