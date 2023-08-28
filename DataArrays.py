# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 09:21:52 2023

@author: grife
"""

class LookupTable():
    def __init__(self,name):
        self.tableName = name;
    
    def addLookupArrays(self, array):
        self.lookupArray = array;
    


class IData():
    
    def setName(self,name):
        self.name = name;
     
    def setDataType(self, dataType):
        self.dataType = dataType;
        
class ScalarData(IData):
    
    def setComponentNumber(self, compNumber):
        self.numComp = compNumber;
    
    
    def setLookUpTable(self, lookUpTable):
        self.LookupTable = lookUpTable;
        
    def setArrayData(self, )
        