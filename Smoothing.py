# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 09:41:58 2023

@author: grife
"""
def get_element_faces(node_ica, ordered = False, toString = False):
    face_ABQ = [[1,2,3,4],
                [5,8,7,6],
                [1,5,6,2],
                [2,6,7,3],
                [3,7,8,4],
                [4,8,5,1]]
    faces = []
    for f in face_ABQ:
        face =[]
        for node_num in f:
            face.append(node_ica[node_num-1])
        if ordered:
            face.sort()
        if toString:
            face = "-".join([str(i) for i in face])
        faces.append(face)        
    return faces
def calculateQualityMetric(coords, UCD_Values=False):
    import numpy as np
    if UCD_Values:
        coords = np.concatenate((coords[-4:],coords[:4]))
    J = 0
    neighbours = [[1,3,4],
            [2,0,5],
            [3,1,6],
            [0,2,7],
            [7,5,0],
            [4,6,1],
            [5,7,2],
            [6,4,3]]     
    for i in range(8):
        neighbourNodes = neighbours[i];
        vertexNodeCoords = coords[i];
        A = np.zeros([3,3]);
        for count,n in enumerate(neighbourNodes):
            neighbourNodeCoords = coords[n];
            A[count] = neighbourNodeCoords - vertexNodeCoords
        if (np.linalg.det(A) > 0):
            inv_A = np.linalg.inv(A)
            kTk = np.linalg.norm(A)*np.linalg.norm(inv_A)
            J += (kTk/3)**2
        else:
            # print("Error: determinant is negative! Element Tangled" ) 
            return -100000
    J = J/8
    J = 1/J
    return J

def calculateAngles(coords):
    import numpy as np
    from math import acos,pi
    coords =  np.concatenate((coords[-4:],coords[:4]))
    neighbours = [[1,3,4],
            [2,0,5],
            [3,1,6],
            [0,2,7],
            [7,5,0],
            [4,6,1],
            [5,7,2],
            [6,4,3]]
    angles = np.zeros([8,3])
    for i in range(8):
        neighbourNodes = neighbours[i];
        vertexNodeCoords = coords[i];
        points = [0,1,2,1,2]
        for count in range(3):
            A = neighbourNodes[points[count]]
            vertexA = coords[A]
            vertexB = vertexNodeCoords
            C = neighbourNodes[points[count + 1]]
            vertexC = coords[C]
            ab = vertexA-vertexB
            cb = vertexC - vertexB
            norm_ab = np.linalg.norm(ab)
            norm_cb = np.linalg.norm(cb)            
            quotient = np.dot(ab,cb)/(norm_ab*norm_cb)
            angle = acos(quotient)*180/pi
            angles[i,count]= angle
    return angles
        
def calculateCurvature(coords,currentNodeCoords):
    dim = 3
    curvature = [0]*dim
    for coord in coords:
       for i in range(dim):
          curvature[i] += (coord[i]/len(coords))
    for i in range(dim):
        curvature[i] = curvature[i] - currentNodeCoords[i]
    return curvature


def get_elements_on_boundary(elementMap):
    #Get elementMapment
    #Create faces
    # Make map key: face(str: node numbers ordered numerically and joined by "-"); value: element number
    # If face strign not in map, check in list of free faces
        # if in free faces, remove from free faces and add to fa
        # else add in free faced
    # elementMap = {1: [1,2,3,4,5,6,7,8,9], 2: [ 10,11,12,13,1,2,3,4], 3: [3,15,16,17,7,18,19,20]} ##TEST INPUT
    print("Locating elements the boundary")
    face_to_elems_map = {}
    surface_face_to_elems_map = {}
    # free_faces = []
    # nodeToElem = create_node_to_elem_map(elementMap)
    for e,ica in elementMap.items():       
        list_of_faces = get_element_faces(ica,ordered = True, toString = True)
        for face_key in list_of_faces:                                             # Create map key 
            if face_to_elems_map.__contains__(face_key):                            # Check if face key already in map
               connected_elements =  face_to_elems_map[face_key]                    # key already in face so append element to array (NOT surface face)
               connected_elements.append(e)
               if surface_face_to_elems_map.__contains__(face_key):                   # If previously classified as a free surface; remove from this map
                   del surface_face_to_elems_map[face_key]
            else:
                face_to_elems_map[face_key] = [e]                                   # If not in map, add to map
                surface_face_to_elems_map[face_key] = e
        
    elements_on_boundary = []
    for face_key,e in surface_face_to_elems_map.items():    
        if not elements_on_boundary.count(e):
                elements_on_boundary.append(e) 
    
    return elements_on_boundary

def get_boundary_element_map(elementMap):
    #Get elementMapment
    #Create faces
    # Make map key: face(str: node numbers ordered numerically and joined by "-"); value: element number
    # If face strign not in map, check in list of free faces
        # if in free faces, remove from free faces and add to fa
        # else add in free faced
    # elementMap = {1: [1,2,3,4,5,6,7,8,9], 2: [ 10,11,12,13,1,2,3,4], 3: [3,15,16,17,7,18,19,20]} ##TEST INPUT
    print("Locating boundary elements")
    face_to_elems_map = {}
    surface_face_to_elems_map = {}
    # free_faces = []
    # nodeToElem = create_node_to_elem_map(elementMap)
    for e,ica in elementMap.items():       
        list_of_faces = get_element_faces(ica,ordered = True, toString = True)
        for face_key in list_of_faces:                                             # Create map key 
            if face_to_elems_map.__contains__(face_key):                            # Check if face key already in map
               connected_elements =  face_to_elems_map[face_key]                    # key already in face so append element to array (NOT surface face)
               connected_elements.append(e)
               if surface_face_to_elems_map.__contains__(face_key):                   # If previously classified as a free surface; remove from this map
                   del surface_face_to_elems_map[face_key]
            else:
                face_to_elems_map[face_key] = [e]                                   # If not in map, add to map
                surface_face_to_elems_map[face_key] = e
        
    boundary_element_map = {}
    for face_key,e in surface_face_to_elems_map.items():    
        faces = get_element_faces(elementMap[e],True,True)
        for face_num,f in enumerate(faces):
            if f == face_key:
                compund_key = "-".join([str(e),str(face_num)])
                boundary_element_map[compund_key] = get_element_faces(elementMap[e])[face_num]
                break    
    
    return boundary_element_map

# def get_volume_elem_to_boundary(elementMap):
#     #Get elementMapment
#     #Create faces
#     # Make map key: face(str: node numbers ordered numerically and joined by "-"); value: element number
#     # If face strign not in map, check in list of free faces
#         # if in free faces, remove from free faces and add to fa
#         # else add in free faced
#     # elementMap = {1: [1,2,3,4,5,6,7,8,9], 2: [ 10,11,12,13,1,2,3,4], 3: [3,15,16,17,7,18,19,20]} ##TEST INPUT
#     from element_functions import calculate_max_number
#     print("Locating boundary elements and surfaes")
#     face_to_elems_map = {}
#     surface_face_to_elems_map = {}
#     # free_faces = []
#     # nodeToElem = create_node_to_elem_map(elementMap)
#     for e,ica in elementMap.items():       
#         list_of_faces = get_element_faces(ica,ordered = True, toString = True)
#         for face_key in list_of_faces:                                             # Create map key 
#             if face_to_elems_map.__contains__(face_key):                            # Check if face key already in map
#                connected_elements =  face_to_elems_map[face_key]                    # key already in face so append element to array (NOT surface face)
#                connected_elements.append(e)
#                if surface_face_to_elems_map.__contains__(face_key):                   # If previously classified as a free surface; remove from this map
#                    del surface_face_to_elems_map[face_key]
#             else:
#                 face_to_elems_map[face_key] = [e]                                   # If not in map, add to map
#                 surface_face_to_elems_map[face_key] = e
        

#     print("Creating volume_elem_to_boundary map")
#     volume_elem_to_boundary = {}
#     b_elem_num = calculate_max_number(elementMap)+1
#     for face_key,e in surface_face_to_elems_map.items():    
#         faces = get_element_faces(elementMap[e],True,True)
#         for face_num,f in enumerate(faces):
#             if f == face_key:
#                 volume_elem_to_boundary[b_elem_num] = e 
#                 b_elem_num += 1
#                 break    
    
#     return volume_elem_to_boundary

def create_surface_connectivity(boundary_element_map):
    ## Create surface connectivity map
    from element_functions import create_node_to_elem_map
    nodeToBoundaryElementMap = create_node_to_elem_map(boundary_element_map)
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

def perform_smoothing(iteration, coeffs, surfaceNodeConnectivity, nodeMap, elementMap, nodeToElemMap,
                      bounds = [-100000,100000,-100000,100000,-100000,100000], inBounds=True):

    from element_functions import value_in_square_bounds
    print("Iteration: " + str(iteration+1))    
    import numpy as np 
    newNodePositions = {}
    badElements = []
    tangled_elements = []
    nodesUnsmoothed = []
    coeff = coeffs[iteration%2]
    for node,connected in surfaceNodeConnectivity.items():
        currentNodeCoords = list(nodeMap[node])
        if (value_in_square_bounds(currentNodeCoords, bounds, inside=inBounds) ):
            coords = []
            for n in connected:
                coords.append(nodeMap[n])
            curvature = calculateCurvature(coords,currentNodeCoords)
            newCoords = [0,0,0]
            for i in range(3):
                newCoords[i] = currentNodeCoords[i] + coeff*curvature[i]
            newNodePositions[node] = newCoords
            for e in nodeToElemMap[node]:
                elemCoords = np.zeros([8,3])
                for count,n in enumerate(elementMap[e]):
                    if newNodePositions.__contains__(n):
                        elemCoords[count] = newNodePositions[n]
                    else:
                        elemCoords[count] = nodeMap[n]
                metric = calculateQualityMetric(elemCoords);
                if metric < 0.2:
                    newNodePositions.pop(node)
                    if metric < -10000:
                        tangled_elements.append(e)
                    else:
                        badElements.append(e)
                        nodesUnsmoothed.append(node)
                    break
                    
    badElements = list(set(badElements))
    nodesUnsmoothed= list(set(nodesUnsmoothed))
    print("Number of unsmoothed nodes: " + str(len(nodesUnsmoothed)))
    print("Number of elements affected: " + str(len(badElements)))
    print("Number of tangled elements: " + str(len(tangled_elements)))
    for node, newcoords in newNodePositions.items():
        nodeMap[node] = newcoords
