"""
Created on Wed May 10 09:11:56 2023

@author: grife
"""

# import numpy as np
# import mne
# import pooch
# import nibabel
from BrainModel import BrainModel
from Maze_Solver import Maze_Solver
from PointCloud import PointCloud
from Mesh import Mesh, MeshUtils

class BrainHexMesh():
    
    def __init__(self):
        self.configured = False
    
    def config(self, configFile):
        # Config
        if configFile.readFile:
            # self.readFile = True
            # self.fileInPath = configFile.fileInPath
            # self.fileIn = configFile.fileIn 
            self.readData = False
            self.data = []       
        elif configFile.readData:
            self.readData = True
            self.data = configFile.data 
            self.readFile = False
            self.fileInPath = ""
            self.fileIn = "" 
        
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
    
    def import_data(self,fileInPath, fileIn): 
        self.brainModel = BrainModel();
        return self.brainModel.import_file(fileInPath,fileIn)
    
    def preprocess(self,data, lesion=True, edemicTissue=True):        
        # Step 2: Determine segmentation of brain model via labels map
        
        # Homogenize labels
        label_number = 0
        if self.material_labels.labelsMap.get('Ventricles',False):
            label_number = self.material_labels.labelsMap['Ventricles'][0]
        
        # if self.material_labels.labelsMap.get('Corpuscallosum',False):
        #     if not (list(np.unique(data)).count(self.material_labels.labelsMap['Corpuscallosum'][0])):
        #         self.add_corpus_callosum(path_cc, filename_cc, current_data)
            
        data = self.material_labels.homogenize_material_labels(data, replace = label_number);      
         
        
        # Coarsen the brain model
        self.VOXEL_SIZE= 1;
        if self.Coarsen:
            print("########## Coarsening data ##########")
            self.VOXEL_SIZE = 2
            data = self.brainModel.coarsen(self.VOXEL_SIZE, data)        
        
        ######### DEBUG COMMENT
        # Clean image removing isolated pixels and small holes
        print("########## Performing cleaning operations on the data ##########")
        
        self.brainModel.clean_voxel_data(data); 
        if lesion:
            self.brainModel.clean_lesion(data)
        if edemicTissue:
            self.brainModel.add_edemic_tissue(data)        
        
        # Remove empty rows/columns and plains from 3D array
        self.brainModel.trim_mesh(data)
        
        # Find and fill voids within model
        print("########## Removing voids from data ##########")
        solver = Maze_Solver(data);
        data = solver.find_and_fill_voids();
        ######### DEBUG COMMENT
        
        # Create CSF layer around GM
        if self.Add_CSF:
            print("########## Adding layers of CSF ##########")
            self.brainModel.add_CSF(data,layers=0);
        return data       
        
    def make_point_cloud(self, data):
            # Create point cloud
        print("########## Creating point cloud from data ##########")
        pointCloud = PointCloud();
        pointCloud.create_point_cloud_from_voxel(data);
        # # Data for visualization
        # pointCloud.view_slice(0, 50); # View slice of point cloud about chosen axis
        # pointCloud.view_point_cloud(); # View full 3D point cloud
        return pointCloud
    
    def get_point_cloud_from_mesh(self, mesh):
        pointCloudFromMesh = PointCloud()
        pointCloudFromMesh.create_point_cloud_from_mesh(mesh.elements, mesh.nodes)
        return pointCloudFromMesh
        
    def make_mesh(self, pc_data):
        print("########## Creating mesh from point cloud ##########")
        mesh = Mesh()        
        mesh.create_mesh_from_Point_Cloud(pc_data,self.VOXEL_SIZE)
        self.clean_mesh(mesh);
        return mesh;
    
    def clean_mesh(self,mesh):
        meshUtils = MeshUtils();
        ######### DEBUG COMMENT
        if self.Add_CSF:
            # Clean grey matter boundary
            mesh.clean_mesh(elementsNotIncluded=[24], replace=24)
            elementsOnBoundary = meshUtils.locate_elements_on_boundary(mesh.elements)
            mesh.replace_outer_region(3, 24, elementsOnBoundary)
        # Clean CSF boundary    
        mesh.clean_mesh()
        # Clean white matter boundary
        mesh.clean_mesh(elementsNotIncluded=[24,3], replace=2)
        ######### 
        # Replace any white matter on boundary with grey matter
        elementsOnBoundary = meshUtils.locate_elements_on_boundary(mesh.elements, elementsNotIncluded = [24])
        mesh.replace_outer_region(2, 3, elementsOnBoundary)        
        
    def add_corpus_callosum(self,path_cc,filename_cc,current_data):
        self.brainModel.add_corpus_callosum(path_cc,filename_cc,current_data)
    
    def smooth_mesh(self, mesh):
        meshUtils = MeshUtils();
        # Optional Boundary Smoothing
        count = 0
        for region in self.Smooth_regions:
            # Smooth regional boundary
            print("########## Smoothing Regions ##########")
            non_whitematter_labels_map = self.material_labels.get_homogenized_labels_map()
            non_whitematter_labels_map.pop(region)
            non_whitematter_labels_map = list(non_whitematter_labels_map.values())
            boundary_element_map = meshUtils.locate_boundary_element_map(mesh.elements, elementsNotIncluded = non_whitematter_labels_map)
            coeffs = self.region_coeffs[count]
            iterations = self.region_iterations[count]
            mesh.smooth_mesh(coeffs, iterations, boundary_element_map)
            count += 1
        # Smooth mesh (excluded CSF)
        if self.Add_CSF:
            print("########## Smoothing mesh excluding CSF ##########")
            # label = self.material_labels.get_homogenized_labels_map()
            # ventricles_label = label.pop("Ventricles")
            boundary_element_map = meshUtils.locate_boundary_element_map(mesh.elements, elementsNotIncluded=[24])
            mesh.smooth_mesh(self.coeffs, self.iterations, boundary_element_map)
            
        # # Smoothing
        # Smooth outer surface of mesh (including CSF)
        print("########## Smoothing global mesh ##########")
        boundary_element_map = meshUtils.locate_boundary_element_map(mesh.elements)
        mesh.smooth_mesh(self.coeffs, self.iterations, boundary_element_map)        
        return mesh
    
        

    def write_to_file(self, mesh):
        
        # Write mesh to file
        for fileType in self.fileoutTypes:
            print("########## Writing data to " + self.fileout + " as a " + fileType.upper() + " file ##########")
            mesh.write_to_file(self.fileoutPath, self.fileout);




