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

def in_hull(p, hull):
    """
    Test if points in `p` are in `hull`

    `p` should be a `NxK` coordinates of `N` points in `K` dimensions
    `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
    coordinates of `M` points in `K`dimensions for which Delaunay triangulation
    will be computed
    """
    from scipy.spatial import Delaunay
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)

    return hull.find_simplex(p)>=0

class BrainModel():
        
    def check_data_dims(self,data):        
        if len(data.shape) == 3:
            return True
        else:
            print('Error, data shape not 3 dimensional!')
            return False       
    
    def coarsen(self, new, original_data):        
        from GridBox import GridBox
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
                                unique, counts = np.unique(gridBox, return_counts=True)
                                num_values = dict(zip(unique, counts))
                                if num_values.__contains__(4):
                                    replacedValue = 4 
                                elif num_values.__contains__(251):
                                    replacedValue = 251                                       
                                elif (len(modes)>1) and (len(modeIndices)>1):
                                    if (modeIndices[0]==0) or (modeIndices[1]==0):
                                        if (modeIndices[0]==0):
                                            replacedValue = modes[modeIndices[1]]
                                        elif (modeIndices[1]==0):
                                            replacedValue = modes[modeIndices[0]]
                                    else:
                                        # xtopL = top_x+1 if top_x<current_dimensions[0] else current_dimensions[0]-1
                                        # ytopL = top_y+1 if top_y<current_dimensions[1] else current_dimensions[1]-1
                                        # ztopL = top_z+1 if top_z<current_dimensions[2] else current_dimensions[2]-1
                                        xbot = x-1 if x>0 else 0
                                        ybot = y-1 if y>0 else 0
                                        zbot = z-1 if z>0 else 0
                                        gridBox = original_data[xbot:top_x, ybot:top_y, zbot:top_z].reshape(-1)
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
        cleaned = self.fill_in_holes(cleaned, structure)
        # print(np.sum(cleaned))
        print("Perfroming binary erosion")
        cleaned = self.binary_erosion(cleaned, structure)
        # print(np.sum(cleaned))
        print("Perfroming binary dilation")
        cleaned = self.binary_dilation(cleaned, structure)
        # print(np.sum(cleaned))
        print("Filling in holes (2)")
        cleaned = self.fill_in_holes(cleaned, structure)
        print(np.sum(cleaned))        
        print("Complete")        
        self.assign_materials_labels(start_data,cleaned)
        return
    
    def fill_in_holes(self, data, structure = None):
        if np.any(structure == None):
            structure = self.create_structure()
        return ndimage.binary_fill_holes(data, structure=structure).astype(int)
    
    def binary_dilation(self, data, structure=None):
        if np.any(structure == None):
            structure = self.create_structure()
        return ndimage.binary_dilation(data, structure=structure).astype(int)
    
    def binary_erosion(self, data, structure=None):
        if np.any(structure == None):
            structure = self.create_structure()
        return ndimage.binary_erosion(data, structure=structure).astype(int)
    
    def create_structure(self):       
        structure = ndimage.generate_binary_structure(3,3)
        return structure
    
    def two_d_cleaning(self, start_data):
        print("Performing 2D cleaning operations on data")
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
                        if replacement_value == None:
                            labelled_data[x,y,z] = 0
                        else:
                            labelled_data[x,y,z] = replacement_value                              
                    elif (labelled_data[x,y,z] != 0) and (end[x,y,z] == 0):
                        labelled_data[x,y,z] = 0; 
                
    def add_CSF(self,data,layers=2):
        from PointCloud import PointCloud
        from scipy.spatial import Delaunay
        
        current_dimensions = data.shape
        newData = np.zeros(current_dimensions, int)
        
        xs,ys,zs = np.where(data == 3)
        for [x,y,z] in np.column_stack((xs,ys,zs)):
            newData[x,y,z] = 3
        
        xs,ys,zs = np.where(data == 2)
        for [x,y,z] in np.column_stack((xs,ys,zs)):
            newData[x,y,z] = 3
            
        xs,ys,zs = np.where(data == 25)
        for [x,y,z] in np.column_stack((xs,ys,zs)):
            newData[x,y,z] = 3
            
        xs,ys,zs = np.where(data == 57)
        for [x,y,z] in np.column_stack((xs,ys,zs)):
            newData[x,y,z] = 3

        inflated_CSF = self.binary_dilation(newData)
        for i in range(layers-1):
            inflated_CSF = self.binary_dilation(inflated_CSF)

        xs,ys,zs = np.where(inflated_CSF == 1)
        current_dimensions = inflated_CSF.shape
        for [x,y,z] in np.column_stack((xs,ys,zs)):
            if newData[x,y,z] == 0:
                newData[x,y,z] = 3 
            
        # Create point cloud
        pointCloud = PointCloud();
        pc = pointCloud.create_point_cloud_of_data(newData);

        xmin_tot,ymin_tot,zmin_tot = [int(p) for p in np.min(pc[:,:3],axis=0)]
        xmax_tot,ymax_tot,zmax_tot = [int(p) for p in np.max(pc[:,:3],axis=0)]

        # ymax_tot = 70

        print("Filling in CSF z-dim")
        for z in range(zmin_tot,zmax_tot+1):
            points = pointCloud.get_slice(2,z);
            points = points[:,:2]    
            hull = Delaunay(points)
            
            xmin,ymin = np.min(points, axis=0)
            xmax,ymax = np.max(points, axis=0)
            
            for x in range(int(xmin),int(xmax+1)):
                for y in range(int(ymin),int(ymax+1)):
                    if (data[x,y,z] == 0) and (y<ymax_tot):
                        if in_hull([x,y], hull):
                            pointCloud.add_point_to_cloud([x,y,z,24])
                            data[x,y,z] = 24
                            newData[x,y,z] = 24

        print("Filling in CSF x-dim")
        for x in range(xmin_tot,xmax_tot+1):
            points = pointCloud.get_slice(0,x);
            points = points[:,1:3]    
            hull = Delaunay(points)
            
            min1d,min2d = np.min(points, axis=0)
            max1d,max2d = np.max(points, axis=0)
            
            for y in range(int(min1d),int(max1d+1)):
                for z in range(int(min2d),int(max2d+1)):
                    if (data[x,y,z] == 0) and (y<ymax_tot):
                        if in_hull([y,z], hull):
                            pointCloud.add_point_to_cloud([x,y,z,24])
                            data[x,y,z] = 24  
                            newData[x,y,z] = 24  
            
        print("Filling in CSF y-dim")
        for y in range(ymin_tot,ymax_tot+1):
            points = pointCloud.get_slice(1,y);
            points = points[:,[0,2]]    
            hull = Delaunay(points)
            
            min1d,min2d = np.min(points, axis=0)
            max1d,max2d = np.max(points, axis=0)
            
            for x in range(int(min1d),int(max1d+1)):
                for z in range(int(min2d),int(max2d+1)):
                    if (data[x,y,z] == 0) and (y<ymax_tot):
                        if in_hull([x,z], hull):     
                            pointCloud.add_point_to_cloud([x,y,z,24])
                            data[x,y,z] = 24 
                            newData[x,y,z] = 24     
    
    
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
        
    
    
        
                                        


            
                                
                            
            
        
    
    
    