# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 13:55:28 2023

@author: grife
"""
from config.Material_Label import Material_Label
from numpy import array
from writers.HeterogeneityConverter import Heterogeneity


class ConfigFile:
    """
    A class used to set up rpeferences for brain creation.

    Attributes
    ----------
    readFile : boolean
        Parameter to determine if data shoudl eb read in from file
    fileInPath : string
        Path to input file
    fileIn : string
        Input file name (.mgz)
    readData : boolean
        Parameter to determine if data should be given explicitly
    data: 3D array 
        Orginal data imported in 
    writeToFile : boolean
        Paremeter to determine if mesh should be written to file      
    fileoutPath : string
        Path to output file, deafault = 'fileInPath'
    fileout : string
        Output file name, deafault = 'fileIn'
    fileoutTypes : array(string)
        File types to which data must be written. Options: 'ucd' | 'vtk' | 'abaqus'
    Coarsen : boolean
        Parameter to determine if data should be corasened to element size = 2
    Add_CSF : boolean
        Parameter to determine if cerebrospinal fluid (CSF) should be added to brain
    csf_layers : int
        Number of layers of CSF to be added
    Smooth : boolean
        Parameter to determine if mesh should be smoothed
    iterations : int
        Number of smoothing iteration to be performed
    coeffs : array(int), len = 2
        Coeeficients to be used in smoothing operation
    Smooth_regions : array(string)
        List of specific regions to smoothed
    region_iterations : array(int)
        List of the number of smoothing iteration to be performed on each region specified
    region_coeffs : array(array(int)) 
        List of the coefficients to be used in smoothing operation on each region specified
    f: file
        Output file of configuration settings
    material_labels : Material_Label
         Object describing material labels to be used in model creation

    Methods
    -------
    add_data(importedData)
        Method to add data explicitly
    openConfigFile()
        Opens configuarion log file and writes configuration setting
    writeToConfig(name, value)
        Method to write additionsl key-value pair data to config file
    closeConfigFile()
        Closes log file
    """
    def __init__(self, inputPath, inputFilename, outputPath, outputFilename):
        """
        Parameters
        ----------
        inputPath : string
            Path to input file.
        inputFilename : TYPE
            Input file name (.mgz).
        outputPath : string
            Path to where output should be saved.
        outputFilename : TYPE
            Output file name.

        """
        self.readFile = True
        self.fileInPath = inputPath
        self.fileIn = inputFilename
        self.readData = False
        self.data = []        
        self.fileoutPath = outputPath
        self.writeToFile = True
        self.fileout = outputFilename
        self.fileoutTypes = ['vtk']  # 'ucd' | 'vtk' | 'abaqus'
        self.preprocess = 'basic'  # 'basic' | 'debug' | 'atrophy' | 'lesion'
        self.Coarsen = True
        self.Add_CSF = 'none'
        self.csf_layers = 1
        self.Smooth = True
        self.iterations = 4
        self.coeffs = [0.6, -0.4]
        self.lesion = True
        self.atrophy = False
        self.ventricles = False
        self.external_cc = True
        self.Refine = False
        
        self.Smooth_regions = ['Lesion']
        self.region_iterations = [4]
        self.region_coeffs =[[0.6, -0.4]]

        self.converter_type = Heterogeneity.NINER
        
        self.material_labels  = Material_Label()
        self.material_labels.addLabelToMap('BrainStem', 16)
        self.material_labels.addLabelToMap('GreyMatter', [3,42]) # Left, Right
        self.material_labels.addLabelToMap('WhiteMatter', [2,41,77]) # Left, Right, WM-hypointensities
        self.material_labels.addLabelToMap('Corpuscallosum', [251,252,253,254,255]) # CC_Posterior, CC_Mid_Posterior, CC_Central, CC_Mid_Anterior, CC_Anterior
        self.material_labels.addLabelToMap('BasalGanglia', [11,50,12,51,13,52,26,58,62,30]) # Caudate(L&R), Putamen(L&R), Palladium(L&R), Accumbens Area(L&R), vessel(L&R)
        self.material_labels.addLabelToMap('Cerebellum', [7,46,8,47]) # WM(L&R), GM(L&R)
        self.material_labels.addLabelToMap('Thalamus' , [10,49,28,60]) # Thalamus(L&R), Ventral DC(L&R)
        self.material_labels.addLabelToMap('Hippocampus' , [17,53]) # Left, Right
        self.material_labels.addLabelToMap('Amygdala' , [18,54]) # Left, Right
        self.material_labels.addLabelToMap('Lesion' , [25,57]) # Left, Right
        self.material_labels.addLabelToMap('CSF' , [24]) # Left, Right, Lateral(L&R), 3rd, 4th ventricles
        self.material_labels.addLabelToMap('Unused' , [4,43,14,15,31,63,85]) # Lateral(L&R), 3rd, 4th NOTE: These labels will be removed to create empty spaces!!
        # OR
        # material_labels.addLabelToMap('Left-Lateral-Ventricle' , [4])
        # material_labels.addLabelToMap('Right-Lateral-Ventricle' , [43])
        # material_labels.addLabelToMap('Left-Inf-Lat-Vent' , [5])
        # material_labels.addLabelToMap('Right-Inf-Lat-Vent ' , [44])
        # material_labels.addLabelToMap('3rd-Ventricle' , [14])
        # material_labels.addLabelToMap('4th-Ventricle' , [15])
        #
        # material_labels.addLabelToMap('Left-choroid-plexus' , [31])
        # material_labels.addLabelToMap('Right-choroid-plexus' , [63])
        # material_labels.addLabelToMap('Optic-Chiasm' , [85])

    def set_materials_label(self, newLabels):
        self.material_labels = newLabels

    def add_data(self, importedData):
        """
        Adds external data to config file

        Parameters
        ----------
        importedData : array
            Array of data to be used

        Raises
        -------
        Error raised if data is not a 3Dimensional array

        """
        self.readFile = False
        self.readData = True
        assert len(array(importedData).shape) == 3
        self.data = importedData
        
    def openConfigFile(self):
        """
        Opens and write data to config log file

        """
        self.f = open("/".join([self.fileoutPath, self.fileout]) + ".txt", 'w')
        if self.readFile:
            self.f.write("Input file: " + self.fileInPath + self.fileIn + "\n")
        else:
            self.f.write("Data input from user")
        self.f.write("Write output to file: " + str(self.writeToFile) + "\n")
        if self.writeToFile:            
            self.f.write("Output written to: " + self.fileoutPath + self.fileout + "\n")
            self.f.write("Output file types: " + ", ".join(self.fileoutTypes) + "\n")
        
        self.f.write("Coarsen: " + str(self.Coarsen) + "\n")
        self.f.write("Add CSF: " + str(self.Add_CSF) + "\n")
        if (self.Add_CSF):
            self.f.write("Layers of CSF: " + str(self.csf_layers) + "\n")
        self.f.write("Smooth global mesh: " + str(self.Smooth) + "\n")
        if (self.Smooth):            
            self.f.write("Iterations: " + str(self.iterations) + ", coeffs: " + ", ".join([str(x) for x in self.coeffs]) + "\n")
        for count,r in enumerate(self.Smooth_regions):
            self.f.write("Smooth region: " + r + "\n")
            self.f.write("Iterations: " + str(self.region_iterations[count]) + ", coeffs: " + \
                         ", ".join([str(x) for x in self.region_coeffs[count]]) + "\n")
        
    def writeToConfig(self, name, value):
        """
        Write additionsl key-value pair data to config file.

        Parameters
        ----------
        name : string
            key string.
        value : string
            value string.

        """
        self.f.write("{}: {}\n".format(name, str(value)))
        
    def closeConfigFile(self):
        """
        Closes config log file

        """
        self.writeToConfig("Regions included", "")
        for r in self.material_labels.get_homogenized_labels_map().keys():
            self.writeToConfig("\t", r)
        self.f.write("COMPLETE\n")
        self.f.close()
