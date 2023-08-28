"""
Created on Wed May 10 09:11:56 2023

@author: grife
"""

import numpy as np
from Maze import Maze
from InverseMaze import InverseMaze
import MeshUtils as mu
import VoxelDataUtils as bm
from Maze_Solver import Maze_Solver
from Mesh import Mesh, QuadElement
from Writer import Writer

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
        
    
    def preprocess(self,data, lesion=False, edemicTissue=1, unusedLabel="unusedLabel"):
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
        if lesion:
            # Smooth and remove features of lesions
            lesion_label = self.material_labels.labelsMap.get("Lesion",[-1000])[0]
            assert lesion_label != -1000, "Lesion label not defined in materials label"
            bm.clean_lesion(data,lesion_label)
        if edemicTissue>0:
            # Add edemicTissue number of layers of tissue around lesion
            self.material_labels.addLabelToMap('EdemicTissue' , [29]);
            bm.add_edemic_tissue(data, edemicTissue, lesion_label, 29)        
        
        # Remove empty rows/columns and plains from 3D array
        # data = bm.trim_data(data)
        
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
        if self.config.Add_CSF:
            print("########## Adding layers of CSF ##########")
            bm.add_CSF(data,layers=self.config.layers);
            
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

        
    def createCSFBoundary(self,mesh,elementNUmber, fullyClosed=False):
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
        if (self.config.Add_CSF):
            e_centroids = np.zeros((max(mesh.elements.keys())+1,4))
            count = 0
            for e_num,element in mesh.elements.items():
                element_centroid = element.calculate_element_centroid()
                e_centroids[e_num] = list(element_centroid) + [element.getMaterial()[0]]
                count += 1
            e_centroids = np.stack(e_centroids, axis = 0)
            
            # mesh.create_node_to_element_connectivity()
            boundary_elements_map = {}
            boundary_number = max(list(mesh.elements.keys()))
            boundaryElementsOnCSF = mesh.locate_boundary_element_map() 
            print("Locating CSF boundary")
            for compoundKey,ica in boundaryElementsOnCSF.items():
                boundary_number += 1;
                ica_nodes = [mesh.nodes[n] for n in ica]
                boundary_element = QuadElement(boundary_number, ica_nodes, mat=[elementNUmber])
                [element_num,face] = [int(x) for x in compoundKey.split("-")]
                [xc,yc,zc,m] = e_centroids[element_num]
                Boundary = fullyClosed
                if m == 24:
                    if not fullyClosed:
                        if m == 24:
                            boundary_number += 1             
                            elements_inline = e_centroids
                            elements_inline = elements_inline[np.where(elements_inline[:,0]==xc)[0],:]
                            elements_inline = elements_inline[np.where(elements_inline[:,2]==zc)[0],:]
                            elements_inline = elements_inline[elements_inline[:, 1].argsort()]
                            current_element_idx, = np.where(elements_inline[:,1]>=yc)
                            if len(current_element_idx)!= 0:
                                current_element_idx = current_element_idx[0]
                                elements_inline = elements_inline[:current_element_idx]
                            mats = list(elements_inline[:,3])
                            Boundary = self.__validCSFBoundary(mats);
                    if Boundary:
                        boundary_elements_map[boundary_number] = boundary_element
                else:
                    boundary_number -= 1;
            return boundary_elements_map;
        return {};
                
    def __validCSFBoundary(self, mats): 
        """
        Private method to determine if elements in selected row are valid to be added to boundary elements
    
        Parameters
        ----------
        mats : array(int)
            array of ints defining element material in row normal to boundary element
            
        Outputs
        ----------
        valid_element: boolean
            
        """   
        grey_matter_label = self.material_labels.get_homogenized_labels_map()['GreyMatter']
        white_matter_label = self.material_labels.get_homogenized_labels_map()['WhiteMatter']
        csf_label = self.material_labels.get_homogenized_labels_map()['CSF']
        
        valid_element = ((mats.count(csf_label) + mats.count(grey_matter_label) + mats.count(white_matter_label)) == len(mats))
        return valid_element;
        
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




