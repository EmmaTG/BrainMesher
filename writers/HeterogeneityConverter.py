from enum import Enum

class Heterogeneity(Enum):
    ONER = 1
    TWOR = 2
    FOURR = 4
    NINER = 9
    NINETEENR = 19

class MaterialsConverter:

    def __init__(self):
        self.converter = {}
        self.materials = [2,3,4,7,10,11,16,17,18,24,25,29,251];

    def add_conversion(self, original_values, new_value):
        try:
            list(original_values)
        except TypeError:
            original_values = [original_values]
        for n in original_values:
            if not (self.converter.get(n), False):
                self.converter[n] = new_value
            else:
                print("This value already exists in converter")

    def convert_materials_labels(self, mesh):
        for e in mesh.elements.values():
            curr_mats = e.getMaterial()
            new_mats = []
            for m in curr_mats:
                if m > 1000:
                    new_mats.append(self.converter[3])
                else:
                    new_mats.append(self.converter[m])
            e.setMaterial(new_mats)


class NineRegionConverter(MaterialsConverter):

    def __init__(self):
        super().__init__()
        for n in self.materials:
            self.converter[n] = n

class NineteenRegionsConverter(MaterialsConverter):
    def __init__(self):
        super().__init__()
        for n in self.materials:
            self.converter[n] = n

    def convert_materials_labels(self, mesh):
        pass;

class HomogenousConverter(MaterialsConverter):

    def __init__(self):
        super().__init__()
        homogenous_materials = [2,3,7,10,11,16,17,18,251]
        for n in homogenous_materials:
            self.converter[n] = 2

        for n in self.materials:
            if not homogenous_materials.count(n):
                self.converter[n] = 24


class TwoRegionConverter(MaterialsConverter):

    def __init__(self):
        super().__init__()
        white_matter = [2,7,10,11,16,17,18,251]
        for n in self.materials:
            if n == 3:
                self.converter[n] = 3
            elif white_matter.count(n):
                self.converter[n] = 2
            else:
                self.converter[n] = 24


class FourRegionConverter(MaterialsConverter):

    def __init__(self):
        super().__init__()
        sub_cortical = [7,10,11,16,17,18]
        for n in self.materials:
            if n == 3:
                self.converter[n] = 3
            elif n == 2:
                self.converter[n] = 2
            elif n == 251:
                self.converter[n] = 251
            elif sub_cortical.count(n):
                self.converter[n] = 7
            else:
                self.converter[n] = 24


class MaterialsConverterFactory:

    @staticmethod
    def get_converter(converter_type):
        """Returns homogeneity converter

        Parameters
        ----------
        converter_type : Heterogeneity
            Enum of heterogeneity type. Options : ONER | TWOR | FOURR | NINER
        """
        if converter_type == Heterogeneity.ONER:
            return HomogenousConverter()
        elif converter_type == Heterogeneity.TWOR:
            return TwoRegionConverter()
        elif converter_type == Heterogeneity.FOURR:
            return FourRegionConverter()
        elif converter_type == Heterogeneity.NINETEENR:
            return NineteenRegionsConverter()
        else:
            return NineRegionConverter()








