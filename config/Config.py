# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 13:55:28 2023

@author: grife
"""
import warnings

from config.Material_Label import Material_Label
from numpy import array
from writers.HeterogeneityConverter import Heterogeneity
import configparser
import os


class ConfigFile:
    """
    A class used to set up preferences for brain creation.
    """

    def __init__(self, configFilePath):
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
        self.f = None
        self.data = []

        config = configparser.ConfigParser()
        config.read(configFilePath)
        default_config = config['DEFAULT']

        self.config_dict = {}

        # file read settings
        self.config_dict['read_file'] = default_config.getboolean('read_file', False)
        self.config_dict['file_in_path'] = default_config.get('inputPath')
        if default_config.get('file_in_path') == '':
            self.config_dict['file_in_path'] = os.getcwd().replace("\\", "/")
        self.config_dict['file_in'] = default_config.get('file_in')

        self.config_dict['read_data'] = default_config.getboolean('read_data', False)
        self.config_dict['external_cc'] = default_config.getboolean('external_cc', False)

        # write settings
        self.config_dict['file_out_path'] = default_config.get('file_out_Path')
        if default_config.get('file_out_path') == '':
            self.config_dict['file_out_path'] = self.config_dict.get('file_in_path')
        self.config_dict['file_out_path'] = "/".join([self.config_dict.get('file_in_path'), "Models"])
        if not os.path.exists(self.config_dict['file_out_path']):
            os.mkdir(self.get('file_out_path'))

        self.config_dict['write_to_file'] = default_config.getboolean('write_to_file', False)
        self.config_dict['fileout'] = default_config.get('fileout', 'unnamed_test')

        types_string = default_config.get('fileout_types')
        types = [x.strip() for x in types_string.split(",")]
        self.config_dict['fileout_types'] = types  # 'ucd' | 'vtk' | 'abaqus'

        # preprocessing option
        self.config_dict['coarsen'] = default_config.getboolean('coarsen', False)

        # model features to include
        lesion_config = config['lesion']
        self.config_dict['lesion'] = lesion_config.getboolean('lesion', False)
        self.config_dict['lesion_layers'] = lesion_config.getint('lesion_layers', 1)

        csf_config = config['csf']
        self.config_dict['add_csf'] = csf_config.get('add_csf', 'full').lower()  # 'none' | 'full' | 'partial'
        self.config_dict['csf_layers'] = csf_config.getint('csf_layers', 1)

        self.config_dict['atrophy'] = default_config.getboolean('atrophy', False)

        # Smoothing features
        self.config_dict['smooth'] = default_config.getboolean('smooth', False)

        if self.config_dict['smooth']:
            smooth_config = config['smooth']
            self.config_dict['iterations'] = smooth_config.getint('iterations')
            smooth_co_effs = smooth_config.get('co_effs')
            self.config_dict['co_effs'] = [float(x.strip()) for x in smooth_co_effs.split(",")]

            smooth_regions = smooth_config.get('smooth_regions')
            if smooth_regions == '':
                self.config_dict['smooth_regions'] = []
            else:
                self.config_dict['smooth_regions'] = [x.strip() for x in smooth_regions.split(",")]

            region_iterations = smooth_config.get('region_iterations')
            if region_iterations == '':
                self.config_dict['region_iterations'] = []
            else:
                self.config_dict['region_iterations'] = [int(x.strip()) for x in region_iterations.split(",")]

            region_co_effs_tmp = smooth_config.get('region_co_effs')
            if region_co_effs_tmp == '':
                self.config_dict['region_co_effs'] = []
            else:
                region_co_effs_tmp = [x.strip() for x in region_co_effs_tmp.replace("[", "").replace("]", "").split(",")]
                region_co_effs = []
                assert (len(region_co_effs_tmp) % 2 == 0)
                for r_count in range(0, len(region_co_effs_tmp), 2):
                    region_co_effs.append([float(region_co_effs_tmp[r_count]), float(region_co_effs_tmp[r_count+1])])
                self.config_dict['region_co_effs'] = region_co_effs

            if self.get('lesion') and 'lesion' not in ", ".join(self.get('smooth_regions')).lower():
                self.config_dict['smooth_regions'].append('lesion')
                self.config_dict['region_iterations'].append(4)
                self.config_dict['region_co_effs'].append([0.6, -0.4])

            self.config_dict['smooth_regions'].reverse()
            self.config_dict['region_iterations'].reverse()
            self.config_dict['region_co_effs'].reverse()

        # Boundary elements
        if 'boundary' in config:
            boundary_config = config['boundary']
            boundary_ele_numbers_tmp = boundary_config.get('boundary_element_numbers')
            if boundary_ele_numbers_tmp == '' or boundary_ele_numbers_tmp is None:
                boundary_element_numbers = []
            else:
                boundary_element_numbers = [int(x.strip()) for x in boundary_ele_numbers_tmp.split(",")]

            excluded_regions_tmp = boundary_config.get('excluded_regions')
            if excluded_regions_tmp == '':
                excluded_regions = []
            else:
                excluded_regions_list = []
                stack = []
                for s in excluded_regions_tmp:
                    s = s.strip()
                    if s == ']':
                        excluded_regions_list.append(stack)
                        stack = []
                    elif s != '[' and s != ',' and s != ' ':
                        stack.append(s)
                excluded_regions = excluded_regions_list

            boundary_tests_tmp = boundary_config.get('boundary_tests') # ['OpenBottomCSF', 'OnlyOnLabel-Ventricles', 'OnlyOnLabel-Lesion']
            if boundary_tests_tmp == '':
                boundary_tests = []
            else:
                boundary_tests = [x.strip() for x in boundary_tests_tmp.split(",")]

            if self.get('lesion') and 'lesion' not in ", ".join(boundary_tests).lower():
                boundary_num = 600
                while boundary_element_numbers.count(boundary_num):
                    boundary_num += 10
                boundary_element_numbers.append(boundary_num)
                excluded_regions.append([])
                boundary_tests.append('OnlyOnLabel-Lesion')

            boundary_element_numbers.reverse()
            self.config_dict['boundary_element_numbers'] = boundary_element_numbers
            excluded_regions.reverse()
            self.config_dict['excluded_regions'] = excluded_regions
            boundary_tests.reverse()
            self.config_dict['boundary_tests'] = boundary_tests

        # Refining features
        self.config_dict['refine'] = default_config.getboolean('refine')
        if self.get('refine'):
            refine_config = config['refine']
            point_refinement = refine_config.getboolean('point', False)
            self.config_dict['refine.point'] = point_refinement
            bounding_box_refinement = refine_config.getboolean('bounding_box', False)
            self.config_dict['refine.bounding_box'] = bounding_box_refinement
            element_refinement = refine_config.getboolean('elements', False)
            self.config_dict['refine.elements'] = element_refinement

            if point_refinement:
                centers_tmp = refine_config.get('centers')
                centers_tmp = centers_tmp.split(']')
                centers = []
                for c in centers_tmp:
                    c.replace("[","")
                    co_ords_string = c.split(",")
                    co_ords = [float(x.strip()) for x in co_ords_string]
                    centers.append(co_ords)
                radii_tmp = refine_config.get('radii')
                radii = [float(x.strip()) for x in radii_tmp.split(",")]
                self.config_dict['refine.centers'] = centers
                self.config_dict['refine.radii'] = radii

            if bounding_box_refinement:
                b_box_tmp = refine_config.get('centers')
                b_box_tmp = b_box_tmp.split(']')
                b_boxes = []
                for bb in b_box_tmp:
                    bb.replace("[", "")
                    co_ords_string = bb.split(",")
                    co_ords = [float(x.strip()) for x in co_ords_string]
                    b_boxes.append(co_ords)
                self.config_dict['refine.bounding_box'] = b_boxes

            if element_refinement:
                element_tmp = refine_config.get('element_numbers')
                elements = [int(x.strip()) for x in element_tmp.split(",")]
                self.config_dict['refine.element_numbers'] = elements


        # Material converter preference
        material_converter = default_config.get('converter_type')
        if material_converter == '1R':
            self.config_dict['converter_type'] = Heterogeneity.ONER
        elif material_converter == '2R':
            self.config_dict['converter_type'] = Heterogeneity.TWOR
        elif material_converter == '4R':
            self.config_dict['converter_type'] = Heterogeneity.FOURR
        else:
            self.converter_type = Heterogeneity.NINER

        # Material labels (PRESET)
        self.material_labels = Material_Label()
        self.material_labels.addLabelToMap('BrainStem', 16)
        self.material_labels.addLabelToMap('GreyMatter', [3, 42])  # Left, Right
        self.material_labels.addLabelToMap('WhiteMatter', [2, 41, 77])  # Left, Right, WM-hypointensities
        self.material_labels.addLabelToMap('CorpusCallosum', [251, 252, 253, 254,
                                                              255])  # CC_Posterior, CC_Mid_Posterior, CC_Central, CC_Mid_Anterior, CC_Anterior
        self.material_labels.addLabelToMap('BasalGanglia', [11, 50, 12, 51, 13, 52, 26, 58, 62,
                                                            30])  # Caudate(L&R), Putamen(L&R), Palladium(L&R), Accumbens Area(L&R), vessel(L&R)
        self.material_labels.addLabelToMap('Cerebellum', [7, 46, 8, 47])  # WM(L&R), GM(L&R)
        self.material_labels.addLabelToMap('Thalamus', [10, 49, 28, 60])  # Thalamus(L&R), Ventral DC(L&R)
        self.material_labels.addLabelToMap('Hippocampus', [17, 53])  # Left, Right
        self.material_labels.addLabelToMap('Amygdala', [18, 54])  # Left, Right
        self.material_labels.addLabelToMap('Ventricles',
                                           [4, 5, 43, 44, 14, 15])  # Lateral(L&R), 3rd, 4th, Inf-Lat-Vent(L&R)

        unused_labels = [0, 31, 63, 85]  # choroid-plexus (L&R), Optic-Chiasm
        if self.get('add_csf'):
            self.material_labels.addLabelToMap('CSF', [24])  # CSF
        else:
            unused_labels.append(24)

        if self.get('lesion'):
            self.material_labels.addLabelToMap('Lesion', [25, 57])  # Left, Right
        else:
            unused_labels += [25, 57]

        self.material_labels.addLabelToMap('Unused', unused_labels)

    def get(self, key):
        result = self.config_dict.get(key, None)
        if result is None:
            warnings.warn("There is no configuration setting with the name: {}".format(key))
        return result

    def set_materials_label(self, newLabels):
        self.material_labels = newLabels

    def add_data(self, importedData):
        """
        Adds external data to config file

        Parameters
        ----------
        importedData :
            Array of data to be used

        Raises
        -------
        Error raised if data is not a 3Dimensional array

        """
        self.config_dict['read_file'] = False
        self.config_dict['read_data'] = True
        assert len(array(importedData).shape) == 3
        self.data = importedData

    def open_config_file(self):
        """
        Opens and write data to config log file

        """
        self.f = open("/".join([self.get('file_out_path'), self.get('fileout')]) + ".txt", 'w')
        if self.config_dict.get('read_file', False):
            self.f.write("Input file: " + self.get('file_in_path') + self.get('file_in') + "\n")
        else:
            self.f.write("Data input from user")
        self.f.write("Write output to file: " + str(self.config_dict.get('write_to_file', False)) + "\n")
        if self.config_dict.get('write_to_file', False):
            self.f.write("Output written to: " + self.config_dict.get('file_out_path') + self.get('fileout') + "\n")
            self.f.write("Output file types: " + ", ".join(self.config_dict.get('fileout_types')) + "\n")

        self.f.write("Coarsen: " + str(self.config_dict.get('coarsen', False)) + "\n")
        self.f.write("Add CSF: " + str(self.get('add_csf')) + "\n")
        if self.config_dict.get('add_csf'):
            self.f.write("Layers of CSF: " + str(self.get('csf_layers')) + "\n")
        self.f.write("Smooth global mesh: " + str(self.get('smooth')) + "\n")
        if self.get('smooth'):
            self.f.write(
                "Iterations: " + str(self.get('iterations')) + ", co_effs: " + ", ".join([str(x) for x in self.get('co_effs')]) + "\n")
            for count, r in enumerate(self.get('smooth_regions')):
                self.f.write("Smooth region: " + r + "\n")
                self.f.write("Iterations: " + str(self.get('region_iterations')[count]) + ", co_effs: " + \
                             ", ".join([str(x) for x in self.get('region_co_effs')[count]]) + "\n")

    def write_to_config(self, name, value):
        """
        Write additional key-value pair data to config file.

        Parameters
        ----------
        name : string
            key string.
        value : string
            value string.

        """
        self.f.write("{}: {}\n".format(name, str(value)))

    def close_config_file(self):
        """
        Closes config log file

        """
        self.write_to_config("Regions included", "")
        for r in self.material_labels.get_homogenized_labels_map().keys():
            self.write_to_config("\t", r)
        self.f.write("COMPLETE\n")
        self.f.close()
