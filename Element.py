# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 14:51:17 2023

@author: grife
"""

from abc import ABC, abstractmethod

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