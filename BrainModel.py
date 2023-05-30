# -*- coding: utf-8 -*-
"""
Created on Fri May 12 11:02:36 2023

@author: grife
"""

import numpy as np
import collections.abc
from scipy import stats
from scipy import ndimage
from Maze_Solver import Maze_Solver
from Vertex import Vertex 
from GridBox import GridBox

class BrainModel():
    
    def __init__(self):        
        self.labelsMap = {}
        self.inverseLabelsMap = {}
        
    def check_data_dims(self,data):        
        if len(data.shape) == 3:
            return True
        else:
            print('Error, data shape not 3 dimensional!')
            return False       
        
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

    def homogenize_labels(self, data):
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
    
    def coarsen(self, new, original_data):        
        print("Coarsening mesh by a factor of " + str(new))
        current_dimensions = original_data.shape
        new_dimensions = [int(p) for p in np.floor(np.array(current_dimensions)/new)];
        newData = np.zeros(new_dimensions, int)
        for x in np.arange(0,(current_dimensions[0]-1),new):
            top_x = x+new
            if (np.sum(original_data[x:top_x,:,:]) > 0):
                for y in np.arange(0,(current_dimensions[1]-1),new):
                    top_y = y+ new
                    if (np.sum(original_data[x:top_x,y:top_y,:]) > 0):
                        for z in np.arange(0,(current_dimensions[2]-1),new):
                            top_z = z+new
                            gridBox = original_data[x:top_x, y:top_y, z:top_z].reshape(-1)
                            if (np.sum(gridBox) > 0):
                                [modes,count] = stats.find_repeats(gridBox)
                                modeIndices, = np.where(count == max(count))
                                modeIndex = modeIndices[0]
                                replacedValue = modes[modeIndex]
                                if (len(modes)>1) and (count[0]==count[1]):
                                    xtopL = top_x+1 if top_x<current_dimensions[0] else current_dimensions[0]-1
                                    ytopL = top_y+1 if top_y<current_dimensions[1] else current_dimensions[1]-1
                                    ztopL = top_z+1 if top_z<current_dimensions[2] else current_dimensions[2]-1
                                    gridBox = original_data[x:xtopL, y:ytopL, z:ztopL].reshape(-1)
                                    [modes,count] = stats.find_repeats(gridBox)
                                    modeIndices, = np.where(count == max(count)) 
                                    modeIndex = modeIndices[0]                                
                                    replacedValue = modes[modeIndex] 
                                newData[int(x/new),int(y/new),int(z/new)] = replacedValue
        return newData        
    
    def create_binary_image(self, current_data):
        current_dimensions = current_data.shape
        newData = np.zeros(current_dimensions, int)
        for x in range(current_dimensions[0]):
            if (np.sum(current_data[x,:,:]) > 0):
                for y in range(current_dimensions[1]):
                    if (np.sum(current_data[x,y,:]) > 0):
                        for z in range(current_dimensions[2]):
                            if (current_data[x,y,z] != 0):
                                newData[x,y,z] = 1
        return newData
    
    def clean_mesh_data(self, start_data):  
        print("Performing cleaning operations on data")
        cleaned = self.create_binary_image(start_data)
        print(np.sum(cleaned))
        structure = ndimage.generate_binary_structure(3,3) 
        print("Filling in holes (1)")
        cleaned = ndimage.binary_fill_holes(cleaned, structure=structure).astype(int)
        # print(np.sum(cleaned))
        print("Perfroming binary erosion")
        cleaned = ndimage.binary_erosion(cleaned, structure=structure).astype(int)
        # print(np.sum(cleaned))
        print("Perfroming binary dilation")
        cleaned = ndimage.binary_dilation(cleaned, structure=structure).astype(int)
        # print(np.sum(cleaned))
        print("Filling in holes (2)")
        cleaned = ndimage.binary_fill_holes(cleaned, structure=structure).astype(int)
        print(np.sum(cleaned))        
        print("Complete")        
        self.assign_materials_labels(start_data,cleaned)
        return
    
    def two_d_cleaning(self, start_data):
        print("Filling in 2D holes")
        cleaned = self.create_binary_image(start_data)
        self.two_d_fill(cleaned)
        self.assign_materials_labels(start_data,cleaned)

    def two_d_fill(self, cleaned):
        hit_structure = np.ones((3,3))
        hit_structure[1,1] = 0
        
        for x in np.arange(0,(cleaned.shape[0])):
            if (np.sum(cleaned[x,:,:]) > 0): 
                old = cleaned[x,:,:]
                sum_old = np.sum(old)
                g = ndimage.binary_hit_or_miss(old, structure1=hit_structure).astype(int)
                loc = np.where(g == 1)
                for r in range(len(loc[0])):
                    px = loc[0][r]
                    py = loc[1][r]
                    cleaned[x,px,py] = 1
                    
        for y in np.arange(0,(cleaned.shape[1])):
            if (np.sum(cleaned[:,y,:]) > 0): 
                old = cleaned[:,y,:]
                sum_old = np.sum(old)
                g = ndimage.binary_hit_or_miss(old, structure1=hit_structure).astype(int)
                loc = np.where(g == 1)
                for r in range(len(loc[0])):
                    px = loc[0][r]
                    py = loc[1][r]
                    cleaned[px,y,py] = 1
                    
        for z in np.arange(0,(cleaned.shape[2])):
            if (np.sum(cleaned[:,:,z]) > 0): 
                old = cleaned[:,:,z]
                sum_old = np.sum(old)
                g = ndimage.binary_hit_or_miss(old, structure1=hit_structure).astype(int)
                loc = np.where(g == 1)
                for r in range(len(loc[0])):
                    px = loc[0][r]
                    py = loc[1][r]
                    cleaned[px,py,z] = 1
        
   
    def assign_materials_labels(self, labelled_data, end):
        assert labelled_data.shape == end.shape
        dimensions = labelled_data.shape;
        problem_areas = []
        for x in range(dimensions[0]):
            for y in range(dimensions[1]):
                for z in range(dimensions[2]):
                    if (labelled_data[x,y,z] == 0) and (end[x,y,z] != 0):
                        box = GridBox(labelled_data,[x,y,z])
                        replacement_value = box.mode() 
                        labelled_data[x,y,z] = replacement_value                              
                    elif (labelled_data[x,y,z] != 0) and (end[x,y,z] == 0):
                        labelled_data[x,y,z] = 0; 
                
    
    def trim_mesh(self, data):
        current_dimensions = data.shape
        start = 0
        end = current_dimensions[0]
        for x in range(current_dimensions[0]):
            if np.sum(data[x,:,:]) != 0:
                start = x
                break
        for rev_x in range(1,current_dimensions[0]):
            x = current_dimensions[0]-rev_x
            if np.sum(data[x,:,:]) != 0:
                end = x
                break
        # print("x trim")
        # print(start,end)
        data = data[start-1:end+1,:,:]
        
        start = 0
        end = current_dimensions[1]
        for y in range(current_dimensions[1]):
            if np.sum(data[:,y,:]) != 0:
                start = y
                break
        for rev_y in range(1,current_dimensions[1]):
            y = current_dimensions[1]-rev_y
            if np.sum(data[:,y,:]) != 0:
                end = y
                break
        # print("y trim")
        # print(start,end)
        data = data[:,start-1:end+1,:]
        
        start = 0
        end = current_dimensions[2]
        for z in range(current_dimensions[2]):
            if np.sum(data[:,:,z]) != 0:
                start = z
                break
        for rev_z in range(1,current_dimensions[2]):
            z = current_dimensions[2]-rev_z
            if np.sum(data[:,:,z]) != 0:
                end = z
                break
        # print("y trim")
        # print(start,end)
        data = data[:,:,start-1:end+1]
        
        return data
        
    
    
        
                                        


            
                                
                            
            
        
    
    
    