class MaterialsConverter:

    def __init__(self):
        self.converter = {}
        self.materials = [2,3,4,7,10,11,16,17,18,24,25,29, 251];

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
                new_mats.append(self.converter[m])
            e.setMaterial(new_mats)


class NineRegionConverter(MaterialsConverter):

    def __init__(self):
        super().__init__()
        for n in self.materials:
            self.converter[n] = n


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
        sub_cortical = [7,10,11,16,17,18,251]
        for n in self.materials:
            if n == 3:
                self.converter[n] = 3
            elif n == 2:
                self.converter[n] = 2
            elif sub_cortical.count(n):
                self.converter[n] = 7
            else:
                self.converter[n] = 24


class MaterialsConverterFactory:

    @staticmethod
    def get_converter(converter_type):
        converter_type = converter_type.upper()
        if converter_type == "1R":
            return HomogenousConverter()
        elif converter_type == "2R":
            return TwoRegionConverter()
        elif converter_type == "4R":
            return FourRegionConverter()
        else:
            return NineRegionConverter()








