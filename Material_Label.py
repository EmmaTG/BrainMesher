# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:08:53 2023

@author: grife
"""

import numpy as np
import collections.abc

class Material_Label:
    
    def __init__(self):
        self.labelsMap = {}
        self.inverseLabelsMap = {}
    
    def clear_labels_mapself(self):
        self.labelsMap = {}
        self.inverseLabelsMap = {}
        
    def addLabelToMap(self, name, numbers):
        if not (hasattr(self, "labelsMap")):
            self.create_labels_map();
        
        storedLabelValues = []
        for x in list(self.labelsMap.values()):
            storedLabelValues += list(x)
        labelsArr = []
        if (self.labelsMap.__contains__(name)):
            labelsArr = self.labelsMap[name]
        else:
            self.labelsMap[name] = labelsArr;          
        try:    
            if (isinstance(numbers, (collections.abc.Sequence, np.ndarray))):
                for num in numbers:                    
                    num = int(num)
                    if storedLabelValues.count(num):
                        self.labelsMap.pop(name)
                        print("Error. The value '{}' already exists in the map".format(num))
                        return -1
            else:
                numbers = [int(numbers)]
        except ValueError or TypeError:
            print("Error: input number '{}' not compatible".format(numbers))
            self.labelsMap.pop(name)
            return -1
            
        
        labelsArr += list(numbers);
        for n in numbers:
            self.inverseLabelsMap[n] = name
        return 0;
    def get_homogenized_labels_map(self):
        h_label_map = {}
        for key,values in self.labelsMap.items():
            h_label_map[key]=values[0]
        return h_label_map
        
    
    def homogenize_material_labels(self, data):
        print("Homogenizing data according to chosen labels")
        current_dimensions = data.shape
        newData = np.zeros(current_dimensions, int)
        for x in range(current_dimensions[0]):
            if (np.sum(data[x,:,:]) > 0):
                for y in range(current_dimensions[1]):
                    if (np.sum(data[x,y,:]) > 0):
                        for z in range(current_dimensions[2]):
                            data_value = data[x,y,z]
                            if self.inverseLabelsMap.__contains__(data_value):
                                label_name = self.inverseLabelsMap[data_value]
                                label_number = self.labelsMap[label_name][0]
                                newData[x,y,z] = label_number
        return newData;
    
    def create_material_sets(self, elements, file_format="abaqus"):
        print("Creating material sets")        
        if (file_format.lower() == "ucd"):    
            elementToMat = {}
            for num,element in elements.items():
                materials = element.properties['mat']
                material_list = []
                for material in materials:
                    mat_name = self.inverseLabelsMap[material]
                    material_list.append(mat_name)
                elementToMat[num] = material_list
            return elementToMat                
        else:            
            materialToElements = {}
            for materialName in self.inverseLabelsMap.values():
                materialToElements[materialName] = []
                
            for num,element in elements.items():
                materials = element.properties['mat']
                for material in materials:
                    mat_name = self.inverseLabelsMap[material]
                    materialToElements[mat_name].append(num)
            return materialToElements
            