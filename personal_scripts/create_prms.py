import csv
import writers.HeterogeneityConverter as hc
from common.brain_regions import BRAINREGIONNAMES, BRAINREGIONNUMBERS, NINEREGIONNUMBERS
from abc import ABC, abstractmethod
from readers.Reader import Reader
from definitions import ROOT_DIR


class MaterialDataBase:

    def __init__(self, type_name, number):
        self.type_name = type_name
        self.number = number
        self.parameter_values_map = {}

    def add_parameter(self, parameter_name, value):
        self.parameter_values_map[parameter_name] = value

    def material_data_string(self):
        start_string = "\tsubsection constitutive@[type={},instance=1,material_id={}]\n".format(self.type_name,
                                                                                                self.number)
        subsection_string = ""
        for p, v in self.parameter_values_map.items():
            # tmpv = v
            # tmpv *= 10
            # exp = 0
            # while (abs(int(tmpv)) < 1):
            #     exp += 1
            #     tmpv += 10
            subsection_string += "\t\t set {} = {}\n".format(p, v)

        end_string = "\tend\n"

        return start_string + subsection_string + end_string


class OneTermModifiedOgden(MaterialDataBase):

    def __init__(self, number):
        super().__init__("modified_one_term_ogden", number)
        self.add_parameter("empirical coefficient", -2)

    def material_data_string(self):
        assert self.parameter_values_map.get("alpha", False)
        assert self.parameter_values_map.get("poisson ratio", False)
        assert self.parameter_values_map.get("empirical coefficient", False)
        assert self.parameter_values_map.get("mu", False)
        return super().material_data_string()

    def set_alpha(self, value):
        self.add_parameter("alpha", value)

    def set_poisson(self, value):
        self.add_parameter("poisson ratio", value)

    def set_mu(self, value):
        self.add_parameter("mu", value)


class OneTermModifiedOgdenAtrophy(OneTermModifiedOgden):

    def __init__(self, number):
        MaterialDataBase.__init__(self, "modified_one_term_ogden_atrophy", number)
        self.add_parameter("empirical coefficient", -2)

    def material_data_string(self):
        assert self.parameter_values_map.get("alpha", False)
        assert self.parameter_values_map.get("poisson ratio", False)
        assert self.parameter_values_map.get("empirical coefficient", False)
        assert self.parameter_values_map.get("mu", False)
        assert self.parameter_values_map.get("normal atrophy rate", False)
        assert self.parameter_values_map.get("enhanced atrophy rate", False)
        return super().material_data_string()

    def set_normal_atrophy(self, value):
        self.add_parameter("normal atrophy rate", value)

    def set_enhanced_atrophy(self, value):
        self.add_parameter("enhanced atrophy rate", value)


class MaterialFactory:

    @staticmethod
    def get_material(material_type, number):
        if material_type == "one_term_ogden":
            return OneTermModifiedOgden(number)
        elif material_type == "one_term_ogden_atrophy":
            return OneTermModifiedOgdenAtrophy(number)
        else:
            raise NotImplementedError

    def material_types(self):
        return "one_term_ogden | one_term_ogden_atrophy"


class CreatePRM:
    TYPE = "tension_compression_testing_device"
    INPUT_DATA = "tension_1"
    all_regions_data_path = ROOT_DIR + "/personal_scripts/materials_data/region_params.csv"
    nine_regions_data_path = ROOT_DIR + "/personal_scripts/materials_data/gov_region_params.csv"

    def __init__(self, path_to_template):
        self.template = None
        self.mesh_material_volumes = None
        self.regions_dict = None
        self.complete_template = None
        self.materials = {}
        self.region_numbers = BRAINREGIONNUMBERS
        self.output_file = None

        self.open_template(path_to_template)

    def open_template(self, path_to_template):
        template_file = open(path_to_template, 'r')
        self.template = template_file.read()
        template_file.close()

    def open_nine_region_csv_data(self):
        self.regions_dict = {}
        f_nine_regions = open(self.nine_regions_data_path, 'r')

        # print("NINE REGIONS")
        nine_regions_reader = csv.DictReader(f_nine_regions, delimiter=";")
        for r in nine_regions_reader:
            inner_dict = {'preconditioned': {}, 'unconditioned': {}}
            region = r.pop('gov_region')
            if region == '':
                continue
            conditioning = r.pop('preconditioned')
            poissons = r.pop('nu')
            self.regions_dict.get(region, self.regions_dict.setdefault(region, inner_dict))
            self.regions_dict[region][conditioning][poissons] = r

        f_nine_regions.close()

    def open_nineteen_region_csv_data(self):
        self.regions_dict = {}
        f_all_regions = open(self.all_regions_data_path, 'r')

        # print("nineTEEN REGIONS")
        all_regions_reader = csv.DictReader(f_all_regions, delimiter=",")
        for r in all_regions_reader:
            inner_dict = {'preconditioned': {}, 'unconditioned': {}}
            region = r.pop('region')
            if region == '':
                continue
            conditioning = r.pop('preconditioned')
            poissons = r.pop('nu')
            self.regions_dict.get(region, self.regions_dict.setdefault(region, inner_dict))
            self.regions_dict[region][conditioning][poissons] = r

        f_all_regions.close()

    def calculate_mesh_volumes(self, mesh):
        self.mesh_material_volumes = {}
        total = 0
        for ele in mesh.elements.values():
            total += 1
            mat = ele.getMaterial()
            for m in mat:
                count = self.mesh_material_volumes.get(m, 0)
                count += 1
                self.mesh_material_volumes[m] = count

    def create_materials(self, mesh, selected_conditioning, seleceted_poissons, heterogeneity, constitutive_type="one_term_ogden"):
        print(heterogeneity.name)
        self.calculate_mesh_volumes(mesh)

        if not heterogeneity == hc.Heterogeneity.NINETEENR:
            self.open_nine_region_csv_data()
            self.region_numbers = NINEREGIONNUMBERS
        else:
            self.open_nineteen_region_csv_data()
            seleceted_poissons = seleceted_poissons.replace(",", ".")

        material_type = hc.MaterialsConverterFactory.get_converter(heterogeneity)
        material_values = [x.value for x in self.region_numbers]

        # Calculate volume ratios based on chosen heterogeneity
        materials_groups = {}
        total_volume_mat_groups = {}
        for i in material_values:
            if self.mesh_material_volumes.get(i):
                converted_type = material_type.converter[i]
                mats = materials_groups.get(converted_type, [])
                mats.append(i)
                materials_groups.update({converted_type: mats})

                volume = total_volume_mat_groups.get(converted_type, 0)
                volume += self.mesh_material_volumes[i]
                total_volume_mat_groups.update({converted_type: volume})

        # Calculate volume averaged values
        for mat, associated_mats in materials_groups.items():
            alpha = 0
            mu = 0
            total = total_volume_mat_groups[mat]

            for i in associated_mats:
                region_abbrev = self.region_numbers(i).name
                region_data = self.regions_dict[region_abbrev][selected_conditioning][seleceted_poissons]
                volume_average = self.mesh_material_volumes[i] / total

                curr_alpha = float(region_data.get('alpha_mean').replace(",", "."))
                alpha += curr_alpha * volume_average

                curr_mu = float(region_data.get('mu_mean').replace(",", "."))
                mu += curr_mu * volume_average

            # Create material objects for each region in volume averaged group
            for i in associated_mats:
                material = MaterialFactory.get_material(constitutive_type, i)
                material.set_alpha(alpha)
                material.set_mu("{:.6e}".format(mu / 1000000.))
                material.set_poisson(seleceted_poissons.replace(",", "."))
                self.materials[i] = material

    def write_materials(self):
        # Create material strings for prm files
        mat_string = ""
        for mat, material in self.materials.items():
            print(material.material_data_string())
            abbrev = self.region_numbers(mat).name
            fullname = BRAINREGIONNAMES[abbrev].value
            mat_string += "\t# {} ({}), {} elements\n".format(fullname.upper(), abbrev, self.mesh_material_volumes[mat])
            mat_string += material.material_data_string()
            mat_string += "\n"

        assert self.template.replace("%materials%", "") != self.template, ".inp file already fully filled in"
        self.template = self.template.replace("%materials%", mat_string)

    def complete_prm(self, configPath, simulationName, simulation_output_name):
        replacement_values = {}
        # Input filearameters
        input_data = self.INPUT_DATA
        type = self.TYPE
        replacement_values["%type%"] = type
        replacement_values["%testing_device%"] = input_data

        inp_file = simulationName + "_UCD.inp"
        replacement_values["%inp_file%"] = inp_file

        output_file = simulation_output_name
        replacement_values["%output%"] = output_file

        assert self.template.replace("%", "") != self.template, ".inp file already fully filled in"
        for key, value in replacement_values.items():
            self.template = self.template.replace(key, str(value))

    def write_prm(self, path_to_output):
        if path_to_output.split(".")[-1] != "prm":
            path_to_output += ".prm"
        self.output_file = open(path_to_output, 'w')
        self.output_file.write(self.template)

    def close_prm(self):
        self.output_file.close()


class CreateAtrophyPRM(CreatePRM):
    INPUT_DATA = "atrophy_1"
    TYPE = "atrophy"

    def create_materials(self, mesh, selected_conditioning, seleceted_poissons, heterogeneity,
                         constitutive_type="one_term_ogden"):
        super().create_materials(mesh, selected_conditioning, seleceted_poissons, heterogeneity,
                                 "one_term_ogden_atrophy")
        for mat in self.materials.values():
            mat_number = mat.number
            if [18, 26, 3, 7, 17, 10, 11, 12, 13].count(mat_number) or mat_number > 1000:
                mat.set_normal_atrophy(0.0015)
                mat.set_enhanced_atrophy(0.0035)
            elif [16, 251, 2].count(mat_number) or 173 <= mat_number <= 178:
                mat.set_normal_atrophy(0.001)
                mat.set_enhanced_atrophy(0.002)

    @staticmethod
    def __read_config_file__(path, simulationName):
        radius = 0;
        center = [0, 0, 0]
        with open("/".join([path, simulationName]) + ".txt") as configReader:
            line = configReader.readline()
            while line != '':
                if (line[:22] == "Concentration radius: "):
                    radius = int(float(line[22:].strip()))
                elif (line[:8] == "Center: "):
                    center = [int(float(x.strip())) for x in line[8:].split(", ")]
                line = configReader.readline()
        return [center, radius]

    def complete_prm(self, configPath, simulationName, simulation_output_name):
        super().complete_prm(configPath, simulationName, simulation_output_name)

        [center, radius] = CreateAtrophyPRM.__read_config_file__(configPath, simulationName)
        replacement_values = {"%center%": str(center).strip('[').strip(']'), "%radius%": radius}

        assert self.template.replace("%", "") != self.template, ".inp file already fully filled in"
        for key, value in replacement_values.items():
            self.template = self.template.replace(key, str(value))
        assert self.template.replace("%", "") == self.template, ".inp file not fully filled in"


class CreateTumorPRM(CreatePRM):

    TYPE = "tumor"
    INPUT_DATA = "growth_exp"


if __name__ == "__main__":
    material_type = hc.MaterialsConverterFactory.get_converter(hc.Heterogeneity.NINETEENR)
    material_values = [x.value for x in BRAINREGIONNUMBERS]

    filename = "csf_brain_centered"
    path = "../IOput/out/csf_brain_centered"

    reader = Reader('vtk')
    reader.openReader(filename + "_VTK", path)
    mesh = reader.getMesh()
    reader.closeReader()

    # filename = "21_region_brain_VTK"
    conditioning = 'preconditioned'
    poissons = '0,49'
    heterogeneity_model = hc.Heterogeneity.ONER
    atrophy_creator = CreateAtrophyPRM("./atrophy_template_folder/atrophy_template_V2.prm")
    atrophy_creator.create_materials(mesh, conditioning, poissons, heterogeneity_model)
    atrophy_creator.write_materials()
    atrophy_creator.complete_prm(path, filename, "{}_atrophy_{}R".format(filename, heterogeneity_model.value))
    output_prm = "/".join([path, "{}_atrophy_{}R".format(filename, heterogeneity_model.value)])
    atrophy_creator.write_prm(output_prm)
