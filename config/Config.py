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
import re


class ConfigFile:
    """
    A class used to set up preferences for brain creation.
    """

    def __init__(self, configFilePath):
        self.MATERIAL_LABELS = None
        self.f = None
        self.data = []
        self.config_dict = {}

        config = configparser.ConfigParser()
        config.read(configFilePath)
        curr_config = config['DEFAULT']

        # File settings
        # Read Settings
        self.config_dict['read_file'] = curr_config.getboolean('read_file', False)
        self.config_dict['file_in_path'] = curr_config.get('file_in_path', '')
        if curr_config.get('file_in_path') == '':
            self.config_dict['file_in_path'] = os.getcwd().replace("\\", "/")
        self.config_dict['file_in'] = curr_config.get('file_in').replace("\\", "/")

        # Write settings
        self.config_dict['write_to_file'] = curr_config.getboolean('write_to_file', False)
        self.config_dict['file_out_path'] = curr_config.get('file_out_path', '').replace("\\", "/")
        if curr_config.get('file_out_path') == '':
            self.config_dict['file_out_path'] = self.config_dict.get('file_in_path')
            self.config_dict['file_out_path'] = "/".join([self.config_dict.get('file_in_path'), "Models"])
        if not os.path.exists(self.config_dict['file_out_path']):
            os.mkdir(self.get('file_out_path'))
        self.config_dict['fileout'] = curr_config.get('fileout', 'unnamed_test')

        types_string = curr_config.get('fileout_types')
        types = [x.strip() for x in types_string.split(",")]
        self.config_dict['fileout_types'] = types  # 'ucd' | 'vtk' | 'abaqus'

        # If corpus callosum is saved as an external file (must be stored as cc.mgz in mri folder)
        self.config_dict['external_cc'] = curr_config.getboolean('external_cc', False)

        # preprocessing options
        self.config_dict['basic_nocsf'] = curr_config.getboolean('basic_nocsf', False)
        self.config_dict['basic_partialcsf'] = curr_config.getboolean('basic_partialcsf', False)
        self.config_dict['basic_fullcsf'] = curr_config.getboolean('basic_fullcsf', False)
        self.config_dict['atrophy'] = curr_config.getboolean('atrophy', False)
        self.config_dict['lesion'] = curr_config.getboolean('lesion', False)

        if self.get('basic_nocsf'):
            curr_config = config['basic_nocsf']
        if self.get('basic_partialcsf'):
            curr_config = config['basic_partialcsf']
        if self.get('basic_fullcsf'):
            curr_config = config['basic_fullcsf']
        if self.get('atrophy'):
            curr_config = config['atrophy']
        if self.get('lesion'):
            curr_config = config['lesion']
            self.config_dict['lesion_layers'] = curr_config.getint('lesion_layers', 1)

        # CSF options
        self.config_dict['add_csf'] = curr_config.getboolean('add_csf', False)
        if self.get('add_csf'):
            self.config_dict['csf_type'] = curr_config.get('csf_type', 'full').lower()  # 'none' | 'full' | 'partial'
            self.config_dict['csf_layers'] = curr_config.getint('csf_layers', 1)

        self.MATERIAL_LABELS = Material_Label()
        if 'materials' in config:
            default_keys = list(config['DEFAULT'])
            materials_config = config['materials']
            for key in materials_config:
                if not default_keys.count(key):
                    label = key
                    values = [int(x.strip()) for x in materials_config[key].split(",")]
                    self.MATERIAL_LABELS.addLabelToMap(label, values)
        else:
            # Material labels (PRESET)
            self.MATERIAL_LABELS.addLabelToMap('brainStem', 16)
            self.MATERIAL_LABELS.addLabelToMap('greymatter', [3, 42])  # Left, Right
            self.MATERIAL_LABELS.addLabelToMap('whitematter', [2, 41, 77])  # Left, Right, WM-hypointensities
            self.MATERIAL_LABELS.addLabelToMap('corpuscallosum', [251, 252, 253, 254,
                                                             255])  # CC_Posterior, CC_Mid_Posterior, CC_Central, CC_Mid_Anterior, CC_Anterior
            self.MATERIAL_LABELS.addLabelToMap('basalganglia', [11, 50, 12, 51, 13, 52, 26, 58, 62,
                                                           30])  # Caudate(L&R), Putamen(L&R), Palladium(L&R), Accumbens Area(L&R), vessel(L&R)
            self.MATERIAL_LABELS.addLabelToMap('cerebellum', [7, 46, 8, 47])  # WM(L&R), GM(L&R)
            self.MATERIAL_LABELS.addLabelToMap('thalamus', [10, 49, 28, 60])  # Thalamus(L&R), Ventral DC(L&R)
            self.MATERIAL_LABELS.addLabelToMap('hippocampus', [17, 53])  # Left, Right
            self.MATERIAL_LABELS.addLabelToMap('amygdala', [18, 54])  # Left, Right
            self.MATERIAL_LABELS.addLabelToMap('ventricles',
                                          [4, 5, 43, 44, 14, 15])  # Lateral(L&R), 3rd, 4th, Inf-Lat-Vent(L&R)

        # unused_labels = [0, 31, 63, 85, 24]  # choroid-plexus (L&R), Optic-Chiasm
        # if self.get('add_csf'):
        #     self.MATERIAL_LABELS.addLabelToMap('csf', [24])  # CSF
        # else:
        #     unused_labels.append(24)

        if self.get('lesion'):
            self.MATERIAL_LABELS.addLabelToMap('lesion', [25, 57])  # Left, Right
            self.MATERIAL_LABELS.addLabelToMap('edemictissue', [29])  # Left, Right
        # else:
        #     unused_labels += [25, 57]

        # self.MATERIAL_LABELS.addLabelToMap('Unused', unused_labels)

        # Smoothing features
        self.config_dict['smooth'] = curr_config.getboolean('smooth', False)

        if self.config_dict['smooth']:
            self.config_dict['iterations'] = curr_config.getint('iterations')
            smooth_co_effs = curr_config.get('co_effs')
            self.config_dict['co_effs'] = [float(x.strip()) for x in smooth_co_effs.split(",")]

            smooth_regions = curr_config.get('smooth_regions')
            if smooth_regions == '':
                self.config_dict['smooth_regions'] = []
            else:
                self.config_dict['smooth_regions'] = [x.strip() for x in smooth_regions.split(",")]

            region_iterations = curr_config.get('region_iterations')
            if region_iterations == '':
                self.config_dict['region_iterations'] = []
            else:
                self.config_dict['region_iterations'] = [int(x.strip()) for x in region_iterations.split(",")]

            region_co_effs_tmp = curr_config.get('region_co_effs')
            if region_co_effs_tmp == '':
                self.config_dict['region_co_effs'] = []
            else:
                region_co_effs_tmp = [x.strip() for x in
                                      region_co_effs_tmp.replace("[", "").replace("]", "").split(",")]
                region_co_effs = []
                assert (len(region_co_effs_tmp) % 2 == 0)
                for r_count in range(0, len(region_co_effs_tmp), 2):
                    region_co_effs.append([float(region_co_effs_tmp[r_count]), float(region_co_effs_tmp[r_count + 1])])
                self.config_dict['region_co_effs'] = region_co_effs

            self.config_dict['smooth_regions'].reverse()
            self.config_dict['region_iterations'].reverse()
            self.config_dict['region_co_effs'].reverse()

        # Boundary elements
        boundary_ele_numbers_tmp = curr_config.get('boundary_element_numbers')
        if boundary_ele_numbers_tmp == '' or boundary_ele_numbers_tmp is None:
            boundary_element_numbers = []
        else:
            boundary_element_numbers = [int(x.strip()) for x in boundary_ele_numbers_tmp.split(",")]

        excluded_regions_tmp = curr_config.get('excluded_regions')
        excluded_regions_tmp = excluded_regions_tmp.strip().replace(" ", "")
        if excluded_regions_tmp == '':
            excluded_regions = []
        else:
            excluded_regions_list = []
            excluded_list = []
            inside_pa = ''
            bracketOpen = False
            for s in excluded_regions_tmp:
                if s == ']' and bracketOpen:
                    excluded_list.append(inside_pa)
                    inside_pa = ''
                    bracketOpen = False
                elif s == '[':
                    bracketOpen = True
                elif bracketOpen:
                    inside_pa += s
            # for entries in excluded_list:
            #     all_labels = self.MATERIAL_LABELS.get_homogenized_labels_map()
            #     if entries == '':
            #         all_labels.clear()
            #     else:
            #         for regions in entries.split(","):
            #             all_labels.pop(regions.strip().lower())
            #     excluded_regions_list.append(list(all_labels.values()))
            excluded_regions = excluded_list

        boundary_tests_tmp = curr_config.get(
            'boundary_tests')  # ['OpenBottomCSF', 'OnlyOnLabel-Ventricles', 'OnlyOnLabel-Lesion']
        boundary_tests = []
        if boundary_tests_tmp != '':
            for x in boundary_tests_tmp.split(","):
                x = x.strip()
                if x == '':
                    boundary_tests.append('none')
                else:
                    boundary_tests.append(x)

        assert (len(boundary_element_numbers) == len(excluded_regions) and
                len(excluded_regions) == len(boundary_tests)), \
            "Boundary mesh data incomplete. Please review config file"

        boundary_element_numbers.reverse()
        self.config_dict['boundary_element_numbers'] = boundary_element_numbers
        excluded_regions.reverse()
        self.config_dict['excluded_regions'] = excluded_regions
        boundary_tests.reverse()
        self.config_dict['boundary_tests'] = boundary_tests

        # Refining features
        self.config_dict['refine'] = curr_config.getboolean('refine')
        if self.get('refine'):
            point_refinement = curr_config.getboolean('point', False)
            self.config_dict['refine.point'] = point_refinement
            bounding_box_refinement = curr_config.getboolean('bounding_box', False)
            self.config_dict['refine.bounding_box'] = bounding_box_refinement
            element_refinement = curr_config.getboolean('elements', False)
            self.config_dict['refine.elements'] = element_refinement

            if point_refinement:
                centers_tmp = curr_config.get('centers')
                centers_tmp = centers_tmp.split(']')
                centers = []
                for c in centers_tmp:
                    c.replace("[", "")
                    co_ords_string = c.split(",")
                    co_ords = [float(x.strip()) for x in co_ords_string]
                    centers.append(co_ords)
                radii_tmp = curr_config.get('radii')
                radii = [float(x.strip()) for x in radii_tmp.split(",")]
                self.config_dict['refine.centers'] = centers
                self.config_dict['refine.radii'] = radii

            if bounding_box_refinement:
                b_box_tmp = curr_config.get('centers')
                b_box_tmp = b_box_tmp.split(']')
                b_boxes = []
                for bb in b_box_tmp:
                    bb.replace("[", "")
                    co_ords_string = bb.split(",")
                    co_ords = [float(x.strip()) for x in co_ords_string]
                    b_boxes.append(co_ords)
                self.config_dict['refine.bounding_box'] = b_boxes

            if element_refinement:
                element_tmp = curr_config.get('element_numbers')
                elements = [int(x.strip()) for x in element_tmp.split(",")]
                self.config_dict['refine.element_numbers'] = elements

        # Material converter preference
        material_converter = curr_config.get('converter_type')
        if material_converter == '1R':
            self.config_dict['converter_type'] = Heterogeneity.ONER
        elif material_converter == '2R':
            self.config_dict['converter_type'] = Heterogeneity.TWOR
        elif material_converter == '4R':
            self.config_dict['converter_type'] = Heterogeneity.FOURR
        else:
            self.converter_type = Heterogeneity.NINER

    def get(self, key):
        result = self.config_dict.get(key, None)
        if result is None:
            warnings.warn("There is no configuration setting with the name:" + key)
            return False
        return result

    def set_materials_label(self, newLabels):
        self.MATERIAL_LABELS = newLabels

    def get_material_value(self, name):
        name = name.lower()
        all_labels = self.MATERIAL_LABELS.get_homogenized_labels_map()
        return all_labels.get(name, -1000)


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
        # self.f = open("/".join([self.get('file_out_path'), self.get('fileout')]) + ".txt", 'w')
        self.f = open("/".join([self.get('file_out_path'), self.get('fileout')]) + ".txt", 'w')

    def write_preamble(self):
        if self.config_dict.get('read_file', False):
            self.f.write("Input file: " + self.get('file_in_path') + self.get('file_in') + "\n")
        else:
            self.f.write("Data input from user")
        self.f.write("Write output to file: " + str(self.config_dict.get('write_to_file', False)) + "\n")
        if self.config_dict.get('write_to_file', False):
            self.f.write("Output written to: " + self.config_dict.get('file_out_path') + self.get('fileout') + "\n")
            self.f.write("Output file types: " + ", ".join(self.config_dict.get('fileout_types')) + "\n")

        self.f.write("Coarsen: " + str(self.config_dict.get('coarsen', False)) + "\n")

        self.f.write("basic configuration settings without CSF: " + str(self.get('basic_nocsf')) + "\n")
        self.f.write("basic configuration settings with partial CSF: " + str(self.get('basic_partialcsf')) + "\n")
        self.f.write("basic configuration settings with full CSF: " + str(self.get('basic_fullcsf')) + "\n")
        self.f.write("atrophy configuration settings: " + str(self.get('atrophy')) + "\n")
        self.f.write("lesion configuration settings: " + str(self.get('lesion')) + "\n")
        if self.get('lesion'):
            self.f.write("layers of edemic tissue: " + str(self.get('lesion_layers')) + "\n")

        self.f.write("Add CSF: " + str(self.get('add_csf')) + "\n")
        if self.config_dict.get('add_csf'):
            self.f.write("Type of CSF added: " + str(self.get('csf_type')) + "\n")
            self.f.write("Layers of CSF: " + str(self.get('csf_layers')) + "\n")

        self.f.write("Smooth global mesh: " + str(self.get('smooth')) + "\n")
        if self.get('smooth'):
            self.f.write(
                "Iterations: " + str(self.get('iterations')) + ", co_effs: " + ", ".join(
                    [str(x) for x in self.get('co_effs')]) + "\n")
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
        for r, values in self.MATERIAL_LABELS.labelsMap.items():
            self.write_to_config("\t", r + ": " + ", ".join([str(x) for x in values]))
        self.f.write("COMPLETE\n")
        self.f.close()
