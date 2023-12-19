# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 13:55:28 2023

@author: grife
"""
import warnings

from config.Material_Label import Material_Label
from writers.HeterogeneityConverter import Heterogeneity
import configparser
import os
from definitions import MATERIALS_9R_PATH, MATERIALS_17R_PATH


class ConfigFile:
    """
    A class used to set up preferences for brain creation.
    """

    def __init__(self, file_in_path, file_in_name, file_out_path, file_out_name, configFilePath="./IOput/model_config.ini", model_type = ''):
        self.MATERIAL_LABELS = None
        self.f = None
        self.data = []
        self.config_dict = {}

        config = configparser.ConfigParser()
        config.read(configFilePath)
        curr_config = config['DEFAULT']

        # Material converter preference
        material_converter = curr_config.get('converter_type')
        if material_converter == '1R':
            self.converter_type = Heterogeneity.ONER
        elif material_converter == '2R':
            self.converter_type = Heterogeneity.TWOR
        elif material_converter == '4R':
            self.converter_type = Heterogeneity.FOURR
        elif material_converter == '19R':
            self.converter_type = Heterogeneity.NINETEENR
        else:
            self.converter_type = Heterogeneity.NINER

        self.MATERIAL_LABELS = Material_Label()
        try:
            if self.converter_type == Heterogeneity.NINETEENR:
                materials_config = configparser.ConfigParser()
                materials_config.read(MATERIALS_17R_PATH)
            else:
                materials_config = configparser.ConfigParser()
                materials_config.read(MATERIALS_9R_PATH)
            materials = materials_config['DEFAULT']
            for key in materials:
                label = key
                values = [int(x.strip()) for x in materials[key].split(",")]
                self.MATERIAL_LABELS.addLabelToMap(label, values)
        except KeyError:
            print("Default materials labels will be used")
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

        # File settings
        # Read Settings
        self.config_dict['read_file'] = True
        file_in_path = file_in_path.replace("\\", "/")
        while file_in_path[-1] == "/":
            file_in_path = file_in_path[:-1]
        self.config_dict['file_in_path'] = file_in_path
        self.config_dict['file_in'] = file_in_name

        # Write settings
        self.config_dict['write_to_file'] = True
        file_out_path = file_out_path.replace("\\", "/")
        while file_out_path[-1] == "/":
            file_out_path = file_out_path[:-1]

        files_to_create = []
        file_out_path_split = file_out_path.split("/")
        while not os.path.exists(file_out_path):
            files_to_create.append(file_out_path_split.pop())
            file_out_path = "/".join(file_out_path_split)

        files_to_create.reverse()
        for file in files_to_create:
            file_out_path = "/".join([file_out_path, file])
            os.mkdir(file_out_path)

        self.config_dict['file_out_path'] = file_out_path
        file_out_name = "unnamed_test" if file_out_name == '' else file_out_name
        self.config_dict['fileout'] = file_out_name

        types_string = curr_config.get('fileout_types')
        types = [x.strip() for x in types_string.split(",")]
        for t in types:
            assert ['ucd', 'vtk', 'abaqus'].count(t), "OUTPUT TYPE {} NOT SUPPORTED".format(t)
        self.config_dict['fileout_types'] = types  # 'ucd' | 'vtk' | 'abaqus'

        # If corpus callosum is saved as an external file (must be stored as cc.mgz in mri folder)
        self.config_dict['external_cc'] = curr_config.getboolean('external_cc', False)
        self.config_dict['segmented_brainstem'] = curr_config.getboolean('segmented_brainstem', False)

        # # preprocessing options
        self.config_dict['model_type'] = model_type
        if model_type in config:
            curr_config = config[model_type]
        if model_type == 'lesion':
            self.config_dict['lesion'] = True
            self.config_dict['lesion_layers'] = curr_config.getint('lesion_layers', 1)
            self.MATERIAL_LABELS.addLabelToMap('lesion', [25, 57])  # Left, Right
            self.MATERIAL_LABELS.addLabelToMap('edemictissue', [29])  # Left, Right

        # CSF options
        self.config_dict['add_csf'] = curr_config.getboolean('add_csf', False)
        if self.get('add_csf'):
            self.config_dict['csf_type'] = curr_config.get('csf_type', 'full').lower()  # 'none' | 'full' | 'partial'
            self.config_dict['csf_layers'] = curr_config.getint('csf_layers', 1)

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
        elif len(boundary_element_numbers) != 0:
            boundary_tests.append('none')

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

    def get(self, key):
        result = self.config_dict.get(key, None)
        if result is None:
            warnings.warn("There is no configuration setting with the name:" + key)
            return False
        return result

    def set(self, key, value):
        result = self.config_dict.get(key, None)
        if result is None:
            warnings.warn("There is no configuration setting with the name:" + key)
        self.config_dict.update({key:value})

    def set_materials_label(self, newLabels):
        self.MATERIAL_LABELS = newLabels

    def get_material_value(self, name):
        name = name.lower()
        all_labels = self.MATERIAL_LABELS.get_homogenized_labels_map()
        return all_labels.get(name, -1000)


    # def add_data(self, importedData):
    #     """
    #     Adds external data to config file
    #
    #     Parameters
    #     ----------
    #     importedData :
    #         Array of data to be used
    #
    #     Raises
    #     -------
    #     Error raised if data is not a 3Dimensional array
    #
    #     """
    #     self.config_dict['read_fil'] = False
    #     self.config_dict['read_dat'] = True
    #     assert len(array(importedData).shape) == 3
    #     self.data = importedData

    def open_config_file(self):
        """
        Opens and write data to config log file

        """
        self.f = open("/".join([self.get('file_out_path'), self.get('fileout')]) + ".txt", 'w')

    def write_preamble(self):
        self.f.write("Input file: " + self.get('file_in_path') + self.get('file_in') + "\n")
        self.f.write("Write output to file: " + str(self.config_dict.get('write_to_file', False)) + "\n")
        if self.config_dict.get('write_to_file', False):
            self.f.write("Output written to: " + self.config_dict.get('file_out_path') + self.get('fileout') + "\n")
            self.f.write("Output file types: " + ", ".join(self.config_dict.get('fileout_types')) + "\n")

        self.f.write("Coarsen: " + str(self.config_dict.get('coarsen', False)) + "\n")

        # self.f.write("basic configuration settings without CSF: " + str(self.get('basic_nocsf')) + "\n")
        self.f.write("model type configuration settings: " + str(self.get('model_type')) + "\n")
        if self.get('model_type') == 'lesion':
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
