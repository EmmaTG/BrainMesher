# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 13:55:28 2023

@author: grife
"""
from config.Material_Label import Material_Label
from numpy import array
from writers.HeterogeneityConverter import Heterogeneity
import configparser


class ConfigFile:
    """
    A class used to set up rpeferences for brain creation.

    Attributes
    ----------
    read_file : boolean
        Parameter to determine if data shoudl eb read in from file
    file_in_path : string
        Path to input file
    file_in : string
        Input file name (.mgz)
    read_data : boolean
        Parameter to determine if data should be given explicitly
    data: 3D array 
        Orginal data imported in 
    write_to_file : boolean
        Paremeter to determine if mesh should be written to file      
    file_out_path : string
        Path to output file, deafault = 'fileInPath'
    fileout : string
        Output file name, deafault = 'fileIn'
    fileout_types : array(string)
        File types to which data must be written. Options: 'ucd' | 'vtk' | 'abaqus'
    coarsen : boolean
        Parameter to determine if data should be corasened to element size = 2
    add_csf : boolean
        Parameter to determine if cerebrospinal fluid (CSF) should be added to brain
    csf_layers : int
        Number of layers of CSF to be added
    smooth : boolean
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
        config = configparser.ConfigParser()
        config.read('model_config.ini')
        sections = config.sections()
        print(sections)
        default_config = config['DEFAULT']
        for key in default_config:
            print(key)
        # file read and write settings
        self.read_file = default_config['read_file']
        self.file_in_path = inputPath
        self.file_in = inputFilename
        self.read_data = False
        self.data = []        
        self.file_out_path = outputPath
        self.write_to_file = True
        self.fileout = outputFilename
        self.fileout_types = ['vtk']  # 'ucd' | 'vtk' | 'abaqus'
        self.external_cc = False

        # preprocessign option
        self.coarsen = True

        # model features to include
        self.add_csf = 'partial'  # 'none' | 'full' | 'partial'
        self.lesion = True
        self.lesion_layers = 1
        self.atrophy = True
        self.csf_layers = 1

        # Smoothing features
        self.smooth = False
        self.iterations = 4
        self.coeffs = [0.6, -0.4]
        self.Smooth_regions = ['Lesion']
        self.region_iterations = [4]
        self.region_coeffs = [[0.6, -0.4]]

        # Boundary elements
        self.boundary_element_numbers = []  # [400, 500, 600]
        self.excluded_regions = [[]]  # [[], [], []]
        self.boundary_tests = []  # ['OpenBottomCSF', 'OnlyOnLabel-Ventricles', 'OnlyOnLabel-Lesion']
        self.boundary_element_numbers.reverse()
        self.excluded_regions.reverse()
        self.boundary_tests.reverse()

        # Refining features
        self.refine = True

        # Material converter preference
        self.converter_type = Heterogeneity.NINER

        # Material labels (PRESET)
        self.material_labels = Material_Label()
        self.material_labels.addLabelToMap('BrainStem', 16)
        self.material_labels.addLabelToMap('GreyMatter', [3, 42])  # Left, Right
        self.material_labels.addLabelToMap('WhiteMatter', [2, 41, 77])  # Left, Right, WM-hypointensities
        self.material_labels.addLabelToMap('CorpusCallosum', [251, 252, 253, 254, 255])  # CC_Posterior, CC_Mid_Posterior, CC_Central, CC_Mid_Anterior, CC_Anterior
        self.material_labels.addLabelToMap('BasalGanglia', [11, 50, 12, 51, 13, 52, 26, 58, 62, 30])  # Caudate(L&R), Putamen(L&R), Palladium(L&R), Accumbens Area(L&R), vessel(L&R)
        self.material_labels.addLabelToMap('Cerebellum', [7, 46, 8, 47])  # WM(L&R), GM(L&R)
        self.material_labels.addLabelToMap('Thalamus', [10, 49, 28, 60])  # Thalamus(L&R), Ventral DC(L&R)
        self.material_labels.addLabelToMap('Hippocampus', [17, 53])  # Left, Right
        self.material_labels.addLabelToMap('Amygdala', [18, 54])  # Left, Right
        self.material_labels.addLabelToMap('Ventricles',
                                           [4, 43, 5, 44, 14, 15])  # Lateral(L&R), 3rd, 4th, Inf-Lat-Vent(L&R)

        unused_labels = [0, 31, 63, 85]  # choroid-plexus (L&R), Optic-Chiasm
        if self.add_csf:
            self.material_labels.addLabelToMap('CSF', [24])  # CSF
        else:
            unused_labels.append(24)

        if self.lesion:
            self.material_labels.addLabelToMap('Lesion', [25, 57])  # Left, Right
        else:
            unused_labels += [25, 57]

        self.material_labels.addLabelToMap('Unused', unused_labels)

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
        self.read_file = False
        self.read_data = True
        assert len(array(importedData).shape) == 3
        self.data = importedData
        
    def openConfigFile(self):
        """
        Opens and write data to config log file

        """
        self.f = open("/".join([self.file_out_path, self.fileout]) + ".txt", 'w')
        if self.read_file:
            self.f.write("Input file: " + self.file_in_path + self.file_in + "\n")
        else:
            self.f.write("Data input from user")
        self.f.write("Write output to file: " + str(self.write_to_file) + "\n")
        if self.write_to_file:
            self.f.write("Output written to: " + self.file_out_path + self.fileout + "\n")
            self.f.write("Output file types: " + ", ".join(self.fileout_types) + "\n")
        
        self.f.write("Coarsen: " + str(self.coarsen) + "\n")
        self.f.write("Add CSF: " + str(self.add_csf) + "\n")
        if (self.add_csf):
            self.f.write("Layers of CSF: " + str(self.csf_layers) + "\n")
        self.f.write("Smooth global mesh: " + str(self.smooth) + "\n")
        if (self.smooth):
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
