# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 13:39:52 2023

@author: grife
"""

def calculate_max_number(map_of_interest):
    """
    Calculates the maximum number in a map of keys with numbers

    Parameters
    ----------
    map_of_interest : map(int,array)
        Map of array with nubes (keys) and associated data (values.

    Returns
    -------
    max_node_num : int
        maximum number in keys rounded up.

    """
    from math import ceil 
    max_node_num = max(list(map_of_interest.keys()))
    len_of_num = (len(str(max_node_num))-1)
    max_node_num = ceil(max_node_num/(10**len_of_num))*(10**len_of_num)
    return max_node_num

def increment_numbers(map_of_old, current_map):
    """
    Increments current numbers based on existsing map numbers

    Parameters
    ----------
    map_of_old : Map(int,array)
        Map of existsign numbers.
    current_map : Map(int,array)
        Map of new numbers.

    Returns
    -------
    replacement_Map : Map(int,array)
        Map of new incremented numbers.

    """
    startnum = calculate_max_number(current_map)
    replacement_Map = {}
    old_nums = list(map_of_old.keys())
    for n in old_nums:
        value = map_of_old.pop(n)
        replacement_Map[n+startnum] = value
    return replacement_Map

def calculate_element_centroid(node_ica, nodeMap):
    """
    Calculates element centroid

    Parameters
    ----------
    node_ica : array(int)
        Array of node numbers associated with element.
    nodeMap : Map(int, array)
        Map of node numbers (keys) to coordinates (values).

    Returns
    -------
    centroid : float
        Element centroid.

    """
    centroid = [0,0,0]
    for n in node_ica:
        coords = nodeMap[n]
        for i in range(3):
            centroid[i] += coords[i]
            
    for i in range(3):
        centroid[i] /= len(node_ica)
    return centroid

def value_in_square_bounds(n_coords, bounds, inside=True):
    """
    Determines in 3D coordinates are within the given square bounds

    Parameters
    ----------
    n_coords : array(float)
        3D coordinates.
    bounds : array(float)
        3D square bounds.
    inbounds : bool, optional
        determine is check is for inside bounds given or outside bounds given.

    Returns
    -------
    bool
        Confirmation of in bounded square or not.

    """
    if inside:
        if ( (n_coords[0] > bounds[0]) and (n_coords[0] < bounds[1]) and
                (n_coords[1] > bounds[2]) and (n_coords[1] < bounds[3]) and
                (n_coords[2] > bounds[4]) and (n_coords[2] < bounds[5]) ):
                return True
    else:
        if ( ((n_coords[0] < bounds[0]) or (n_coords[0] > bounds[1])) or
                ((n_coords[1] < bounds[2]) or (n_coords[1] > bounds[3])) or
                ((n_coords[2] < bounds[4]) or (n_coords[2] > bounds[5])) ):
                return True
    return False
        
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
    for coords in nodeMap.values():
        newCoords = np.matmul(np.array(coords),R)
        coords[0] = round(newCoords[0],3)
        coords[1] = round(newCoords[1],3)
        coords[2] = round(newCoords[2],3)

        
def location_rotation(loc,mapOfValues,reverse=False):
    if (loc==4) or (loc==6):
        if not reverse:
            rotate_mesh(mapOfValues,1,90)
        else:
            rotate_mesh(mapOfValues,1,-90)            
    if loc==5:
        if not reverse:
            rotate_mesh(mapOfValues,1,-90)
        else:
            rotate_mesh(mapOfValues,1,90)
            
def scale_mesh(nodeMap,scale=[1,1,1],reduce=True):
    if reduce:
        scale[0] = 1/scale[0]
        scale[1] = 1/scale[1]
        scale[2] = 1/scale[2]
    for coords in nodeMap.values():
        coords[0] = coords[0]*scale[0]
        coords[1] = coords[1]*scale[1]
        coords[2] = coords[2]*scale[2]


def translate_mesh(nodeMap,distance=[1,1,1]):
    for coords in nodeMap.values():
        coords[0] = coords[0]+distance[0]
        coords[1] = coords[1]+distance[1]
        coords[2] = coords[2]+distance[2]
def round_to_interval(value,interval):
    return interval*round(value/interval)    

def get_edges(ica):
    edges = []
    edge_classification = [[0,1],[1,2],[2,3],[3,0],
                           [4,5],[5,6],[6,7],[7,4],
                           [0,4],[1,5],[2,6],[3,7]]
    for edge in edge_classification:
        edges.append("-".join(str(x) for x in sorted([ica[edge[0]],ica[edge[1]]])))
    return edges
        
def get_face_normal(nodeMap,ica):
    import numpy as np
    v1v2 = np.array(nodeMap[ica[0]]) -np.array(nodeMap[ica[1]])
    v1v3 = np.array(nodeMap[ica[0]]) -np.array(nodeMap[ica[2]])
    normal = np.cross(v1v2,v1v3)
    return normal/np.linalg.norm(normal)
            
def remove_elements_in_square_bounds(elementMap,nodeMap,bounds):
    extraction_Elements = []
    extraction_Nodes = []
    for e_num,ica in elementMap.items():
        centroid = calculate_element_centroid(ica, nodeMap);
        if value_in_square_bounds(centroid, bounds):
            extraction_Elements.append(e_num)
            for n in ica:
                # if not extraction_Nodes.count(n):
                extraction_Nodes.append(n)
    
    return extraction_Elements, extraction_Nodes

def create_node_to_elem_map(elementMap):
    "Creating node to element connectivity"
    nodeToElemMap = {}
    for e,ica in elementMap.items():        
        for node in ica:
            if nodeToElemMap.__contains__(node):
                elements = nodeToElemMap[node]
            else:
                elements = []
            elements.append(e)
            nodeToElemMap[node] = elements  
    return nodeToElemMap
    
def replace_duplicate_nodes(nodes_changed_map, elems_changed_maps, nodeMap, decimal_places=6, bounds=[-1000]):
    
    print("Replacing duplicated nodes")
    max_node_num = calculate_max_number(nodeMap)
    nodeLocationMapX = {}    
    for n,coords in nodeMap.items():
        [xcoord, ycoord, zcoord] = [round(i,decimal_places) for i in coords]
        if nodeLocationMapX.__contains__(xcoord):
            nodeLocationMapY = nodeLocationMapX[xcoord]
        else:                    
            nodeLocationMapY = {}                
        if nodeLocationMapY.__contains__(ycoord):
            nodeLocationMapZ = nodeLocationMapY[ycoord]
        else:                    
            nodeLocationMapZ = {}                
        if nodeLocationMapZ.__contains__(zcoord):
            nodeLocationMapZ[zcoord] = n
        else:
            nodeLocationMapZ[zcoord] = n
            nodeLocationMapY[ycoord] = nodeLocationMapZ
            nodeLocationMapX[xcoord] = nodeLocationMapY
            
    swapMap = {}
    for n,orig_coords in nodes_changed_map.items():
        # if not nodeMap.__contains__(n):
        [xcoord, ycoord, zcoord] = [round(i,decimal_places) for i in orig_coords]
        if nodeLocationMapX.__contains__(xcoord):
            nodeLocationMapY = nodeLocationMapX[xcoord]
        else:                    
            nodeLocationMapY = {}                
        if nodeLocationMapY.__contains__(ycoord):
            nodeLocationMapZ = nodeLocationMapY[ycoord]
        else:                    
            nodeLocationMapZ = {}                
        if nodeLocationMapZ.__contains__(zcoord):
            swapMap[n] = nodeLocationMapZ[zcoord]
        else:
            nodeLocationMapZ[zcoord] = n
            nodeLocationMapY[ycoord] = nodeLocationMapZ
            nodeLocationMapX[xcoord] = nodeLocationMapY
            nodeMap[n] = orig_coords
    
    for elem_num,ica in elems_changed_maps.items():
        for pos,n in enumerate(ica):
            ica[pos] = n+max_node_num
            if swapMap.__contains__(n+max_node_num):
                ica[pos] = swapMap[n+max_node_num]
    return  

# import ABQ_UCD_handling as functions

# path = 'C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\MatLab\\INPFiles\\'
# filename = "ellipse_in_cube_2mm_coarse"
# nodeMap, elementMap, elementSetsMap = functions.readABQ(path, filename)

# rotate_mesh(nodeMap, axis=1, degrees=90)
# scale_mesh(nodeMap,scale=[6,6,6])
# functions.writeABQ(path, filename+"_rotated_scaled", nodeMap, elementMap)














