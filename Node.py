# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 14:50:29 2023

@author: grife
"""


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