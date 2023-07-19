# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 09:41:58 2023

@author: grife
"""

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


def perform_smoothing(iteration, coeffs, surfaceNodeConnectivity, nodeMap, elementICAMap, nodeToElemMap,
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
                for count,n in enumerate(elementICAMap[e]):
                    if newNodePositions.get(n,False):
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
