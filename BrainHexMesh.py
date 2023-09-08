"""
Created on Wed May 10 09:11:56 2023

@author: grife
"""

import numpy as np

from writers.HeterogeneityConverter import MaterialsConverterFactory
from voxel_data.void_filler import Maze, InverseMaze, Maze_Solver
import voxel_data.voxel_data_utils as bm
from voxel_data.PreProcessingActions import IpreProcessAction
from mesh.Mesh import Mesh
from mesh.Element import QuadElement
from writers.Writer import Writer
from voxel_data import PreProcessingActions as pp
from readers import Importer as inp

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
        smoothing options include: specific regional smoothing, smoothing boundar exclusind CSF, global outer mesh smoothing
    write_to_file(mesh)
        write mesh data to file accordign to filetypes specifed in config file
    __validCSFBoundary(mats), private
        private method to determine if row of elements is valid for use as a boundary csf element
    """
    
    def __init__(self):
        self.configured = False
    
    def setConfig(self, configFile):
        """
        Imports configuration file defining preferences w.r.t. model creation.

        Parameters
        ----------
        configFile : ConfigFile
            The configuration settings for the model

        """
        self.config = configFile
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
        assert self.configured, "config file has not been set for this. Please run config(cf -> ConfigFile) before importing data"
        if self.config.readData:
            return self.config.data
        if (path == "") and (file == ""):
            path = self.config.fileInPath
            file = self.config.fileIn
        importer = inp.ImportFromFile(path, file)
        return importer.getData()
    
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
        data = self.material_labels.homogenize_material_labels(data, replace = label_number) 
        return data

    def preprocessFactory(self, data, config, csfConfig=""):
        csf_function = None
        if csfConfig == "full":
            csf_function = bm.add_full_CSF
        elif csfConfig == "partial":
            csf_function = bm.add_partial_CSF

        if config == "lesion":
            lesion_label = self.material_labels.labelsMap.get("Lesion", [-1000])[0]
            assert lesion_label != -1000, "Lesion label not defined in materials label"
            # Add edemicTissue number of layers of tissue around lesion
            self.material_labels.addLabelToMap('EdemicTissue', [29])
            return self.preprocess(data, pp.CleanLesion(lesion_label), pp.AddEdemicTissue(),
                                   add_CSF_Function=csf_function)
        elif config == "atrophy":
            return self.preprocess(data, add_CSF_Function=csf_function)
        elif config == "basic":
            if self.config.lesion:
                lesion_label = self.material_labels.labelsMap.get("Lesion", [-1000])[0]
                assert lesion_label != -1000, "Lesion label not defined in materials label"
                return self.preprocess(data, pp.CleanLesion(lesion_label), add_CSF_Function=csf_function)
            return self.preprocess(data, add_CSF_Function=csf_function)
        else:
            raise KeyError("{} is not a valid pre processing key".format(config))
                
    def preprocess(self, data, *args, add_CSF_Function=None):
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
        self.VOXEL_SIZE= 1
        if self.config.Coarsen:
            print("########## Coarsening data ##########")
            self.VOXEL_SIZE = 2
            data = bm.coarsen(self.VOXEL_SIZE, data)

        print("########## Performing cleaning operations on the data ##########")
            # Clean image removing isolated pixels and small holes
        bm.clean_voxel_data(data)
        
        for a in args:
            if isinstance(a, IpreProcessAction):
                a.performAction(data)

        change = True
        iterationCount = 0        
        while(change):
            iterationCount += 1
            print("### Iteration number " + str(iterationCount))
            # Find and fill erroneous voids within model
            print("########## Removing voids from data ##########")
            maze = Maze.Maze(data)
            solver = Maze_Solver.Maze_Solver(maze)
            voidsToFill = solver.find_voids()
            data = solver.fill_voids(voidsToFill)
            
            print("########## Removing disconnected regions from data ##########")
            cont_data = bm.create_binary_image(data)
            cont_data = cont_data-1
            cont_data = cont_data*(-1)
            
            maze2 = InverseMaze.InverseMaze(cont_data)
            solver2 = Maze_Solver.Maze_Solver(maze2)
            voidsToFill = solver2.find_voids() 
            
            for key in voidsToFill:
                [x,y,z]  = [int(x) for x in key.split("-")]
                data[x,y,z] = 0
                
            change = bm.clean_voxel_data(data)
        
        # Create CSF layer around GM
        if (self.config.Add_CSF) and (not add_CSF_Function is None):
            print("########## Adding layers of CSF ##########")
            add_CSF_Function(data, layers=self.config.layers)
            
            print("########## Checking for voids in csf data ##########")
            csfMaze = Maze.Maze(data)
            solver3 = Maze_Solver.Maze_Solver(csfMaze)
            voidsToFill = solver3.find_voids()
            data = solver3.fill_voids(voidsToFill)
        else:
            self.config.Add_CSF = False
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
            mesh object
            
        Errors
        ----------
        Error raised if voxel size has not been specified, i.e. pre-processing has not been performed
        """
        
        assert hasattr(self, "VOXEL_SIZE"), "Voxel size has not been specified, ensure" + \
            " you have run preprocess() method on voxel data before progressing"
            
        print("########## Creating mesh from point cloud ##########")
        mesh = Mesh()        
        mesh.create_mesh_from_Point_Cloud(pc_data,self.VOXEL_SIZE)
        return mesh
    
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
            mesh object to be cleaned
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
            mesh object
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
            boundary_number += 1
            ica_nodes = [mesh.nodes[n] for n in ica]
            boundary_element = QuadElement(boundary_number, ica_nodes, mat=[elementNUmber])
            [element_num,face] = [int(x) for x in compoundKey.split("-")]
            # [xc,yc,zc,m] = e_centroids[element_num]
            Boundary = True
            if not boundaryTest is None:
                Boundary = boundaryTest.validElement(element_num)
            if Boundary:
                boundary_elements_map[boundary_number] = boundary_element
            else:
                boundary_number -= 1
        return boundary_elements_map
        
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
            mesh object to be smoothed
        """
        # Optional Boundary smoothing
        count = 0
        for region in self.config.Smooth_regions:
            # Smooth regional boundary
            print("########## smoothing Regions ##########")
            non_whitematter_labels_map = self.material_labels.get_homogenized_labels_map()
            non_whitematter_labels_map.pop(region)
            non_whitematter_labels_map = list(non_whitematter_labels_map.values())
            coeffs = self.config.region_coeffs[count]
            iterations = self.config.region_iterations[count]
            mesh.smooth_mesh(coeffs, iterations, elementsNotIncluded = non_whitematter_labels_map)
            count += 1
        # Smooth mesh (excluded CSF)
        if self.config.Add_CSF:
            print("########## smoothing mesh excluding CSF ##########")
            # label = self.material_labels.get_homogenized_labels_map()
            # ventricles_label = label.pop("Ventricles")
            mesh.smooth_mesh(self.config.coeffs, self.config.iterations, elementsNotIncluded=[24])
            
        # # smoothing
        # Smooth outer surface of mesh (including CSF)
        print("########## smoothing global mesh ##########")
        mesh.smooth_mesh(self.config.coeffs, self.config.iterations)        
        # return mesh    

    def __convert_heterogeneity__(self, mesh):
        converter = MaterialsConverterFactory.get_converter(self.config.converter_type)
        converter.convert_materials_labels(mesh)

        self.config.writeToConfig("Heterogeneity level", self.config.converter_type)
        print("Heterogeneity level: {}".format(self.config.converter_type))

    def write_to_file(self, mesh):
        """
        Wriet mesh data to various fiel types as specified in config file

        Parameters
        ----------
        mesh: Mesh
            mesh object to be written
        """        
        # Write mesh to file
        self.__convert_heterogeneity__(mesh)
        for fileType in self.config.fileoutTypes:
            print("########## Writing data as a " + fileType.upper() + " file ##########")             
            writer = Writer()
            writer.openWriter(fileType, self.config.fileout, self.config.fileoutPath)
            writer.writeMeshData(mesh)
            writer.closeWriter()




