"""
Created on Wed May 10 09:11:56 2023

@author: grife
"""

import numpy as np
# import mne
# import pooch
import nibabel
from BrainModel import BrainModel
from Maze_Solver import Maze_Solver
from PointCloud import PointCloud
from Mesh import Mesh

class BrainHexMesh():
    
    def __init__(self):
        self.configured = False
    
    def config(self, configFile):
        # Config
        if configFile.readFile:
            self.readFile = True
            self.fileInPath = configFile.fileInPath
            self.fileIn = configFile.fileIn        
        elif configFile.readData:
            self.readData = True
            self.data = self.data 
        
        self.writeToFile = configFile.writeToFile
        self.fileoutPath = configFile.fileoutPath
        self.fileout = configFile.fileout
        self.fileoutTypes = configFile.fileoutTypes
        self.Coarsen = configFile.Coarsen
        self.Add_CSF = configFile.Add_CSF
        self.Smooth = configFile.Smooth
        self.iterations = configFile.iterations
        self.coeffs = configFile.coeffs
        
        self.Smooth_regions = configFile.Smooth_regions
        self.region_iterations = configFile.region_iterations
        self.region_coeffs = configFile.region_coeffs
        self.configured = True
        self.material_labels = configFile.material_labels
    
    def run(self):
        if self.configured:
            self.configured = False
            data = None
            if self.readFile:
                # Step 1: Using freesurfer and 'recon-all' create mri outputs. Ensure aseg.mgz is created.
                t1_file = "\\".join([self.fileInPath,self.fileIn])
                t1 = nibabel.load(t1_file)
                # t1.orthoview()
                data = np.asarray(t1.dataobj)
            elif self.readData:
                data = self.data
            else:
                print("Error")
                return -1
            # Step 2: Determine segmentation of brain model via labels map
            
            # Homogenize labels
            data = self.material_labels.homogenize_material_labels(data);
            
            brainModel = BrainModel()
            
            # Coarsen the brain model
            voxel_size= 1;
            if self.Coarsen:
                print("########## Coarsening data ##########")
                voxel_size = 2
                data = brainModel.coarsen(voxel_size, data)
            
            # Clean image removing isolated pixels and small holes
            print("########## Performing cleaning operations on the data ##########")
            brainModel.clean_mesh_data(data);
            brainModel.two_d_cleaning(data);
            
            # Remove empty rows/columns and plains from 3D array
            brainModel.trim_mesh(data)
            
            # Find and fill voids within model
            print("########## Removing voids from data ##########")
            solver = Maze_Solver(data);
            data = solver.find_and_fill_voids();
            
            # Create CSF layer around GM
            if self.Add_CSF:
                print("########## Adding layers of CSF ##########")
                brainModel.add_CSF(data,layers=1);
                      
            # Create point cloud
            print("########## Creating point cloud from data ##########")
            pointCloud = PointCloud();
            pc = pointCloud.create_point_cloud_of_data(data);
            # # Data for visualization
            # pointCloud.view_slice(0, 50); # View slice of point cloud about chosen axis
            # pointCloud.view_point_cloud(); # View full 3D point cloud
            
            # Create mesh from point cloud
            print("########## Creating mesh from point cloud ##########")
            mesh = Mesh(pointCloud.pcd,voxel_size)
            if self.Add_CSF:
                mesh.clean_mesh(elementsNotIncluded=[24], replace=24)
            mesh.remove_outer_white_matter(2, 3)
            
            # # Smoothing
            # Smooth outer surface of mesh (including CSF)
            if self.Smooth:
                print("########## Smoothing global mesh ##########")
                mesh.smooth_mesh(self.coeffs, self.iterations)
            
                # Smooth mesh (excluded CSF)
                if self.Add_CSF:
                    print("########## Smoothing mesh excluding CSF ##########")
                    mesh.smooth_mesh(self.coeffs, self.iterations, elementsNotIncluded=[24])
                
                # Optional Boundary Smoothing
                    for region in self.Smooth_regions:
                        # Smooth regional boundary
                        print("########## Smoothing white matter ##########")
                        non_whitematter_labels_map = self.material_labels.get_homogenized_labels_map()
                        non_whitematter_labels_map.pop(region)
                        non_whitematter_labels_map = list(non_whitematter_labels_map.values())
                        mesh.smooth_mesh(self.region_coeffs, self.region_iterations, elementsNotIncluded=non_whitematter_labels_map)
            
            if self.writeToFile:
                self.write_to_file(mesh)
            return data, pointCloud, mesh

    def write_to_file(self, mesh, material_labels=None, boundaryElementMap={}):
        
        if material_labels == None:
            material_labels = self.material_labels
        # Write mesh to file
        for fileType in self.fileoutTypes:
            print("########## Writing data to " + self.fileout + " as a " + fileType.upper() + " file ##########")
            mesh.write_to_file(self.fileoutPath, self.fileout, labels_map=material_labels, 
                               filetype=fileType, boundaryElementMap=boundaryElementMap);




