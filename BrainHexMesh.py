"""
Created on Wed May 10 09:11:56 2023

@author: grife
"""

import numpy as np
from Maze import Maze
from InverseMaze import InverseMaze
import VoxelDataUtils as bm
from Maze_Solver import Maze_Solver
from Mesh import Mesh, QuadElement
from Writer import Writer    
from abc import ABC, abstractmethod;
from scipy import ndimage, stats
import warnings

class IpreProcessAction(ABC):
    """
    Command interface for additional preprocessing step outisde of cleaning
    """
    
    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'performAction') and 
                callable(subclass.performAction) or 
                NotImplemented)
    
    @abstractmethod 
    def performAction(self, data, **kwargs):
        raise NotImplementedError;
        
class CleanLesion(IpreProcessAction):
    """
    Cleans added lesion to ensure smooth lesion boundary
    """
    
    def __init__(self,lesion_label):
        self.label = lesion_label
    
    def performAction(self, current_data):
        print("########## Cleaning Lesion ##########")
        ## TODO: Better creation of featureless lesion (possibly the same way featurless csf is created)
        newData = bm.create_binary_image(current_data, search=self.label)
        lesionOGSize = np.sum(newData)
        if lesionOGSize > 0:
            print("Lesion element size before: {}".format(lesionOGSize))
        
            structure1 = ndimage.generate_binary_structure(3,1)
            structure2 = ndimage.generate_binary_structure(3,3)
            
            newData = ndimage.binary_dilation(newData, structure=structure2, iterations=1).astype(int)
            newData = ndimage.binary_erosion(newData, structure=structure1, iterations=1).astype(int)
            newData = ndimage.binary_dilation(newData, structure=structure2, iterations=1).astype(int)
            newData = ndimage.binary_erosion(newData, structure=structure2, iterations=2).astype(int)
        
            hit_structure1 = np.ones((2,2,2))
            hit_structure1[0,1,0] = 0
            hit_structure2 = np.rot90(hit_structure1)
            hit_structure3 = np.rot90(hit_structure1, k= 2)
            hit_structure4 = np.rot90(hit_structure1, k= 3)
            hit_structure5 = np.rot90(hit_structure1,axes=(1,2))
            hit_structure6 = np.rot90(hit_structure5, k= 1)
            hit_structure7 = np.rot90(hit_structure5, k= 2)
            hit_structure8 = np.rot90(hit_structure5, k= 3)
            total_count = 0;
            count = 1
            iteration = 0
            while (count > 0 and iteration < 10):
                iteration += 1
                count = 0
                count += bm.hit_and_miss_3d_2x2x2(newData, hit_structure1, fill=1)
                count += bm.hit_and_miss_3d_2x2x2(newData, hit_structure2, fill=1)
                count += bm.hit_and_miss_3d_2x2x2(newData, hit_structure3, fill=1)
                count += bm.hit_and_miss_3d_2x2x2(newData, hit_structure4, fill=1)
                count += bm.hit_and_miss_3d_2x2x2(newData, hit_structure5, fill=1)
                count += bm.hit_and_miss_3d_2x2x2(newData, hit_structure6, fill=1)
                count += bm.hit_and_miss_3d_2x2x2(newData, hit_structure7, fill=1)
                count += bm.hit_and_miss_3d_2x2x2(newData, hit_structure8, fill=1)
                total_count += count
            print("Lesion cleaned after {} iterations and {} elements added".format(iteration,total_count))
            
            count_removed = 0
            current_dimensions = current_data.shape
            for x in range(current_dimensions[0]):
                if (np.any(newData[x,:,:] == 1) or np.any(current_data[x,:,:] == self.label)):
                    for y in range(current_dimensions[1]):
                        if (np.any(newData[x,y,:] == 1) or np.any(current_data[x,y,:] == self.label)):
                            for z in range(current_dimensions[2]):
                                if newData[x,y,z] == 1:
                                    current_data[x,y,z] = self.label
                                if newData[x,y,z] == 0 and (current_data[x,y,z] == self.label):
                                    current_data[x,y,z] = 0;
            print("previous lesion replaced with non-lesion: {}".format(count_removed))
            finalSize = np.sum(newData)
            print("Lesion element size after: {}".format(finalSize))
            if (finalSize == 0):
                warnings.warn("No lesion elements found in data after cleaning")
        else:
            warnings.warn("No lesion elements found in data")
        
class AddEdemicTissue(IpreProcessAction):
    """
    Add layers of edemic tissue to outside of lesion
    """
    
    def __init__(self, lesion_label = 25, edemic_tissue_label = 29, layers = 1 ):
        self.label = lesion_label
        self.edemic_tissue_label = edemic_tissue_label
        self.layers = layers;
    
    def performAction(self, current_data):
        print("########## Creating edemic Tissue ##########")
        newData = bm.create_binary_image(current_data, search = self.label)   
        if np.sum(newData)>0:            
            structure2 = ndimage.generate_binary_structure(3,3)
            newData = ndimage.binary_dilation(newData, structure=structure2, iterations=self.layers).astype(int)
            current_dimensions = current_data.shape
            for x in range(current_dimensions[0]):
                if (np.any(newData[x,:,:] == 1) or np.any(current_data[x,:,:] == self.label)):
                    for y in range(current_dimensions[1]):
                        if (np.any(newData[x,y,:] == 1) or np.any(current_data[x,y,:] == self.label)):
                            for z in range(current_dimensions[2]):
                                if newData[x,y,z] == 1 and current_data[x,y,z] != self.label:
                                    current_data[x,y,z] = self.edemic_tissue_label
        else:
            warnings.warn("No edemic tissue added as no lesion elements were found in data")
            
class CoarsenData(IpreProcessAction):
    """
    Coarsen voxels in data
    """
    
    def __init__(self, voxel_size = 2):
        self.voxel_size = voxel_size
    
    def performAction(self, data):
        print("Coarsening mesh by a factor of " + str(self.voxel_size))
        original_data = np.copy(data);
        current_dimensions = original_data.shape
        new_dimensions = [int(p) for p in np.floor(np.array(current_dimensions)/self.voxel_size)];
        data *= 0;
        for x in np.arange(0,(current_dimensions[0]-1),self.voxel_size):
            top_x = x + self.voxel_size
            if (np.sum(original_data[x:top_x,:,:]) > 0):
                for y in np.arange(0,(current_dimensions[1]-1),self.voxel_size):
                    top_y = y + self.voxel_size
                    if (np.sum(original_data[x:top_x,y:top_y,:]) > 0):
                        for z in np.arange(0,(current_dimensions[2]-1),self.voxel_size):
                            top_z = z + self.voxel_size
                            gridBox = original_data[x:top_x, y:top_y, z:top_z].reshape(-1)
                            if (np.sum(gridBox) > 0):
                                [modes,count] = stats.find_repeats(gridBox)
                                modeIndices, = np.where(count == max(count))
                                modeIndex = modeIndices[0]
                                replacedValue = modes[modeIndex]
                                unique, counts = np.unique(gridBox, return_counts=True)
                                num_values = dict(zip(unique, counts))
                                if num_values.get(4,False):
                                    replacedValue = 4 
                                elif num_values.get(251,False):
                                    replacedValue = 251                                       
                                elif (len(modes)>1) and (len(modeIndices)>1):
                                    if (modeIndices[0]==0) or (modeIndices[1]==0):
                                        if (modeIndices[0]==0):
                                            replacedValue = modes[modeIndices[1]]
                                        elif (modeIndices[1]==0):
                                            replacedValue = modes[modeIndices[0]]
                                    else:
                                        xbot = x-1 if x>0 else 0
                                        ybot = y-1 if y>0 else 0
                                        zbot = z-1 if z>0 else 0
                                        gridBox = original_data[xbot:top_x, ybot:top_y, zbot:top_z].reshape(-1)
                                        [modes,count] = stats.find_repeats(gridBox)
                                        modeIndices, = np.where(count == max(count)) 
                                        modeIndex = modeIndices[0]                                
                                        replacedValue = modes[modeIndex]                                 
                                data[int(x/self.voxel_size),int(y/self.voxel_size),int(z/self.voxel_size)] = replacedValue
        


class IBoundaryTest(ABC):
    """
    Interface for tests on CSF boundaries. If certain criteria have to be met for a boundary element to be added.
    """
    
    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'validElement') and 
                callable(subclass.validElement) or 
                NotImplemented)
    
    @abstractmethod
    def validElement(self, element_num):
        raise NotImplementedError 
        
class OnlyOnLabel(IBoundaryTest):
    """
    Boundary test to add elements that are only attached to CSF elements 
    """
    
    def __init__(self, mesh, label):
        self.mesh = mesh;
        self.label = label;
    
    def validElement(self, element_num):
        mat = self.mesh.elements[element_num].getMaterial();
        if mat.count(self.label):            
            return True;
        return False; 
    
class OpenBottomCSF(IBoundaryTest):
    """
    Boundary test to add elements that are attached to CSF elements and also are not directly below the subcortical structures
    """
    
    def __init__(self, mesh):
        e_centroids = np.zeros((max(mesh.elements.keys())+1,4))
        count = 0
        for e_num,element in mesh.elements.items():
            element_centroid = element.calculate_element_centroid()
            e_centroids[e_num] = list(element_centroid) + [element.getMaterial()[0]]
            count += 1
        self.e_centroids = np.stack(e_centroids, axis = 0)
    
    def validElement(self, element_num):
        [xc,yc,zc,m] = self.e_centroids[element_num]
        if m == 24:            
            elements_inline = self.e_centroids
            elements_inline = elements_inline[np.where(elements_inline[:,0]==xc)[0],:]
            elements_inline = elements_inline[np.where(elements_inline[:,2]==zc)[0],:]
            elements_inline = elements_inline[elements_inline[:, 1].argsort()]
            current_element_idx, = np.where(elements_inline[:,1]>=yc)
            if len(current_element_idx)!= 0:
                current_element_idx = current_element_idx[0]
                elements_inline = elements_inline[:current_element_idx]
            mats = list(elements_inline[:,3])
            
            grey_matter_label = 3
            white_matter_label = 2
            csf_label = 24
            
            valid_element = ((mats.count(csf_label) + mats.count(grey_matter_label) + mats.count(white_matter_label)) == len(mats))
            return valid_element;
        return False;



class BrainHexMesh():
    """
    A class used to facilitate 3D brain model creation.

    Attributes
    ----------
    configured : boolean
        Indicator whether config file has been given
    material_labels : Material_Label
         Object describing material labels to be used in model creation
    config: Config
        Configuration file for model creation

    Methods
    -------
    config(configFile)
        configures class according to prefernces in configFile
    import_data(path="", file="")
        imports data according to configuration file
    preprocess(data, lesion=False, edemicTissue=1, unusedLabel="Unused")
        performs all preprocessing steps:\n
        1. Homogenize data labels\n
        2. Corasen model, if requested\n
        3. Clean voxel data\n
        4. Clean lesion, if requested\n
        4. Add layers of edemic tissue, if edemicTissue specified\n
        5. Trim data\n
        6. Identify and fill erroneous voids in the data\n
        7. Add layers of CSF, if requested
    make_mesh(pc_data)
        create mesh from pointcloud data. Raises error if VOXEL_SIZE has not been set
    clean_mesh(mesh,wm=True)
        cleans inputed mesh to ensure better smoothing and remove erroneously connected edges.
        cleans outer boundary, grey matter and, if wm=True, white matter boundary
    createCSFBoundary(mesh,elementNUmber)
        creates boundary elements on CSF with element material= elementNUmber
    add_region(cc_data,current_data, region_value)
        overwrite current_data with data given in cc_data and replaces label with region_value
    smooth_mesh(mesh)
        performs Laplacian smoothing of mesh based on config file preferences.
        Smoothing options include: specific regional smoothing, smoothing boundar exclusind CSF, global outer mesh smoothing
    write_to_file(mesh)
        write mesh data to file accordign to filetypes specifed in config file
    __validCSFBoundary(mats), private
        private method to determine if row of elements is valid for use as a boundary csf element
    """
    
    def __init__(self):
        self.configured = False
    
    def config(self, configFile):
        """
        Imports configuration file defining preferences w.r.t. model creation.

        Parameters
        ----------
        configFile : Config
            The configuration settings for the model

        """
        self.config = configFile;
        self.material_labels = configFile.material_labels 
        self.configured = True  
    
    def import_data(self, path="", file=""):
        """
        Return voxel data object by either importing it from a file or returning it from the config file

        Parameters
        ----------
        path : str, options
            path to input file if differnt from config file
        file : str, options
            input file if differnt from config file
            
        Outputs
        ----------
        3D array:
            3D data array of voxels with label numbers
            
        Errors
        ----------
        Error raised if config file not initialized beforehand 
        """
        assert self.configured, "Config file has not been set for this. Please run config(cf -> ConfigFile) before importing data" 
        if self.config.readData:
            return self.config.data;
        if (path == "") and (file == ""):
            path = self.config.fileInPath
            file = self.config.fileIn
        return bm.import_file(path,file)
    
    def homogenize_data(self, data, unusedLabel):
        """
        Homogenizes data based on materials labels

        Parameters
        ----------
        data : 3D array
            voxel data
        unusedLabel : string
            label name of voxel label numbers to be removed at the end of the model creation
            
        Outputs
        ----------
        3D array:
            3D data array of voxels with label numbers
        """
        # Homogenize labels
        label_number = 0
        if self.material_labels.labelsMap.get(unusedLabel,False):
            label_number = self.material_labels.labelsMap[unusedLabel][0]
                    
        # Replace regions with multiple labels with only one label, if label is not required replace with unused/0
        data = self.material_labels.homogenize_material_labels(data, replace = label_number); 
        return data;
        
    
    def preprocess(self, data, *args, unusedLabel="unusedLabel", add_CSF_Function=None):
        """
        Performs all preprocessing steps on voxel data:\n
        2. Corasen model, if requested\n
        3. Clean voxel data\n
        4. Clean lesion, if requested\n
        4. Add layers of edemic tissue, if edemicTissue specified\n
        5. Trim data\n
        6. Identify and fill erroneous voids in the data\n
        7. Add layers of CSF, if requested\n

        Parameters
        ----------
        data : 3D array
            voxel data
        lesion : boolean, optional
            parameter to decide if lesion should be made featureless
            Default is False
        edemicTissue : number, optional
            parameter toif and then number of layers of edemic tissue to be added
            Default is 1
        unusedLabel : string, optional
            label name of voxel label numbers to be removed at the end of the model creation
            Default is "unusedLabel"
            
        Outputs
        ----------
        3D array:
            3D data array of voxels with label numbers
            
        Errors
        ----------
        Error raised lesion smoothing requested but no lesion material exists in mateirls labels object
        """                     
        
        # Coarsen the brain model
        self.VOXEL_SIZE= 1;
        if self.config.Coarsen:
            print("########## Coarsening data ##########")
            self.VOXEL_SIZE = 2
            data = bm.coarsen(self.VOXEL_SIZE, data)  
            
        print("########## Performing cleaning operations on the data ##########")
            # Clean image removing isolated pixels and small holes
        bm.clean_voxel_data(data); 
        
        for a in args:
            a.performAction(data);
        
        change = True;
        iterationCount = 0;        
        while(change):
            iterationCount += 1
            print("### Iteration number " + str(iterationCount))
            # Find and fill erroneous voids within model
            print("########## Removing voids from data ##########")
            maze = Maze(data);
            solver = Maze_Solver(maze);
            voidsToFill = solver.find_voids();
            data = solver.fill_voids(voidsToFill);
            
            print("########## Removing disconnected regions from data ##########")
            cont_data = bm.create_binary_image(data);
            cont_data = cont_data-1;
            cont_data = cont_data*(-1);
            
            maze2 = InverseMaze(cont_data)
            solver2 = Maze_Solver(maze2);
            voidsToFill = solver2.find_voids(); 
            
            for key in voidsToFill:
                [x,y,z]  = [int(x) for x in key.split("-")]
                data[x,y,z] = 0
                
            change = bm.clean_voxel_data(data);
        
        # Create CSF layer around GM
        if (self.config.Add_CSF) and (not add_CSF_Function is None):
            print("########## Adding layers of CSF ##########")
            add_CSF_Function(data, layers=self.config.layers)
            
            print("########## Checking for voids in csf data ##########")
            csfMaze = Maze(data)
            solver3 = Maze_Solver(csfMaze);
            voidsToFill = solver3.find_voids();
            data = solver3.fill_voids(voidsToFill);
        return data       
        
    def make_mesh(self, pc_data):
        """
        Makes mesh from voxel point Cloud data

        Parameters
        ----------
        pc_data : nx4 array
            point data of n points with columns 0:3 specifying coordinates 
            and column 3 giving the material label
            
        Outputs
        ----------
        mesh: Mesh
            Mesh object
            
        Errors
        ----------
        Error raised if voxel size has not been specified, i.e. pre-processing has not been performed
        """
        
        assert hasattr(self, "VOXEL_SIZE"), "Voxel size has not been specified, ensure" + \
            " you have run preprocess() method on voxel data before progressing"
            
        print("########## Creating mesh from point cloud ##########")
        mesh = Mesh()        
        mesh.create_mesh_from_Point_Cloud(pc_data,self.VOXEL_SIZE)
        return mesh;
    
    def clean_mesh(self,mesh, wm=True):
        """
        Cleans mesh by removing/replacing poorly connected elements and/or 
        nodes on the followign boundaries:\n
        1. Grey matter, if CSF added\n
        2. Outer boundary\n
        3. White matter, if specified

        Parameters
        ----------
        mesh: Mesh
            Mesh object to be cleaned
        wm: boolean, optional
            parameter to indiciate white matter boundary cleaning needed
            Default is True
            
        """
        
        if self.config.Add_CSF:
            # Clean grey matter boundary
            print("####### Cleaning grey matter boundary #######")
            mesh.clean_mesh(elementsNotIncluded=[24], replace=24)
            elementsOnBoundary = mesh.locate_elements_on_boundary()
            mesh.replace_outer_region(3, 24, elementsOnBoundary)
            
        # Clean outer boundary
        print("####### Cleaning outer boundary #######")    
        mesh.clean_mesh()
        
        # Clean white matter boundary
        if wm:
            print("####### Cleaning white matter boundary #######")
            mesh.clean_mesh(elementsNotIncluded=[24,3], replace=2)
        
        # Replace any white matter on boundary with grey matter
        elementsOnBoundary = mesh.locate_elements_on_boundary(elementsNotIncluded = [24])
        mesh.replace_outer_region(2, 3, elementsOnBoundary)        


    def createBoundary(self, mesh, elementNUmber, elementsNotIncluded = [], boundaryTest=None):
        """
        Creates CSF boundary elements. 
        Does not create CSF Boundary elements below subcortical structures at base of brain
    
        Parameters
        ----------
        mesh: Mesh
            Mesh object
        elementNUmber: int
            materials number for boundary elements
            
        Outputs
        ----------
        boundary_elements_map: Map(int,QuadElement)
            Map of boundary elements to thier element numbers        
            
        """            
        # mesh.create_node_to_element_connectivity()
        boundary_elements_map = {}
        boundary_number = max(mesh.elements.keys()) if len(mesh.boundaryElements) == 0 else max(mesh.boundaryElements.keys())
        boundaryElements = mesh.locate_boundary_element_map(elementsNotIncluded = elementsNotIncluded) 
        print("Locating CSF boundary")
        for compoundKey,ica in boundaryElements.items():
            boundary_number += 1;
            ica_nodes = [mesh.nodes[n] for n in ica]
            boundary_element = QuadElement(boundary_number, ica_nodes, mat=[elementNUmber])
            [element_num,face] = [int(x) for x in compoundKey.split("-")]
            # [xc,yc,zc,m] = e_centroids[element_num]
            Boundary = True;
            if not boundaryTest is None:
                Boundary = boundaryTest.validElement(element_num)
            if Boundary:
                boundary_elements_map[boundary_number] = boundary_element
            else:
                boundary_number -= 1;
        return boundary_elements_map;
        
    def add_region(self,cc_data,current_data, region_value):
        """
        ADd region to voxel data by overwriting current data

        Parameters
        ----------
        cc_data : 3D array
            voxel data of new region to be added
        current_data : 3D array
            voxel data to eb overwritten
        region_value : 3D array
            label to be assigned to overwritten data
            
        Errors
        ----------
        Error raised cc_data and current_data are not the same size
        """
        
        assert np.all(cc_data.shape == current_data.shape),"Region associated with added data is not the same size as the current data set"
        bm.override_voxel_data(cc_data,current_data, region_value)
    
    def smooth_mesh(self, mesh):
        """
        Smooth various regions of the mesh as specified in config file

        Parameters
        ----------
        mesh: Mesh
            Mesh object to be smoothed
        """
        # Optional Boundary Smoothing
        count = 0
        for region in self.config.Smooth_regions:
            # Smooth regional boundary
            print("########## Smoothing Regions ##########")
            non_whitematter_labels_map = self.material_labels.get_homogenized_labels_map()
            non_whitematter_labels_map.pop(region)
            non_whitematter_labels_map = list(non_whitematter_labels_map.values())
            coeffs = self.config.region_coeffs[count]
            iterations = self.config.region_iterations[count]
            mesh.smooth_mesh(coeffs, iterations, elementsNotIncluded = non_whitematter_labels_map)
            count += 1
        # Smooth mesh (excluded CSF)
        if self.config.Add_CSF:
            print("########## Smoothing mesh excluding CSF ##########")
            # label = self.material_labels.get_homogenized_labels_map()
            # ventricles_label = label.pop("Ventricles")
            mesh.smooth_mesh(self.config.coeffs, self.config.iterations, elementsNotIncluded=[24])
            
        # # Smoothing
        # Smooth outer surface of mesh (including CSF)
        print("########## Smoothing global mesh ##########")
        mesh.smooth_mesh(self.config.coeffs, self.config.iterations)        
        # return mesh    
        

    def write_to_file(self, mesh):
        """
        Wriet mesh data to various fiel types as specified in config file

        Parameters
        ----------
        mesh: Mesh
            Mesh object to be written
        """        
        # Write mesh to file
        for fileType in self.config.fileoutTypes:
            print("########## Writing data as a " + fileType.upper() + " file ##########")             
            writer = Writer()
            writer.openWriter(fileType, self.config.fileout, self.config.fileoutPath)
            writer.writeMeshData(mesh)
            writer.closeWriter();




