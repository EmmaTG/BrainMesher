# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 15:33:29 2023

@author: grife
"""
import numpy as np

class BrainVoxelData():  
    
    @staticmethod      
    def create_binary_image(current_data, search=-1):
        current_dimensions = current_data.shape
        newData = np.zeros(current_dimensions, int)
        for x in range(current_dimensions[0]):
            if (np.sum(current_data[x,:,:]) > 0):
                for y in range(current_dimensions[1]):
                    if (np.sum(current_data[x,y,:]) > 0):
                        for z in range(current_dimensions[2]):
                            if (search== -1) and (current_data[x,y,z] != 0):
                                newData[x,y,z] = 1
                            elif (search != -1) and (current_data[x,y,z] == search):
                                newData[x,y,z] = 1
                                
        return newData
    
    @staticmethod      
    def get_bounding_box(data):
        maxValues = [-100000,-100000,-100000]
        minValues = [100000,100000,100000]
        current_dimensions = data.shape
        for x in range(current_dimensions[0]):
            if (np.sum(data[x,:,:]) > 0):
                for y in range(current_dimensions[1]):
                    if (np.sum(data[x,y,:]) > 0):
                        for z in range(current_dimensions[2]):
                            if (data[x,y,z] == 1):
                                coords = [x,y,z]
                                for d in range(3):
                                    if (maxValues[d] < coords[d]):
                                        maxValues[d] = coords[d]
                                    if (minValues[d] > coords[d]):
                                        minValues[d] = coords[d]
        return maxValues + minValues