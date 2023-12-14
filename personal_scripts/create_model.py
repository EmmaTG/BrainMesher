import brain_creation
from config.Config import ConfigFile
import writers.HeterogeneityConverter as hc
from readers.Reader import Reader
from personal_scripts.create_prms import CreateAtrophyPRM, CreateTumorPRM

# Model type options: basic_fullcsf, basic_partilacsf, basic_nocsf, atrophy, lesion

# ####### ATROPHY CODE #######
# # Preferences are defined in ConfigFile
# path_in = "../IOput/in"
# file_name_in = "aseg.mgz"
# path_out = "../IOput/out/atrophy_example"
# file_name_out = "atrophy_example"
#
# config = ConfigFile(path_in, file_name_in, path_out, file_name_out,
#                     configFilePath="../IOput/model_config.ini", model_type='atrophy')
#
# config.set("smooth", False)
#
# mesh = brain_creation.run(config)
#
# conditioning = 'preconditioned'
# poissons = '0,49'
# for heterogeneity_model in [hc.Heterogeneity.ONER, hc.Heterogeneity.TWOR, hc.Heterogeneity.FOURR, hc.Heterogeneity.NINER]:
#     # heterogeneity_model = hc.Heterogeneity.ONER
#
#     atrophy_creator = CreateAtrophyPRM("./atrophy_template_folder/atrophy_template_V2.prm")
#     atrophy_creator.create_materials(mesh, conditioning, poissons, heterogeneity_model)
#     atrophy_creator.write_materials()
#     atrophy_creator.complete_prm(path_out, file_name_out, "{}_atrohpy_{}R".format(file_name_out, heterogeneity_model.value))
#     output_prm = "/".join([path_out, "{}_atrohpy_{}R".format(file_name_out, heterogeneity_model.value)])
#     atrophy_creator.write_and_close_prm(output_prm)
#
# brain_creation.write(mesh, path_out, file_name_out, ['vtk'])
#
#
# print("COMPLETE")
# print("Files written to",path_out)

####################################################################################################################
# ####### TUMOR LESION CODE #######
path_in = "../IOput/in"
file_name_in = "aseg.mgz"
path_out = "../IOput/out/tumor_example"
file_name_out = "tumor_example"

config = ConfigFile(path_in, file_name_in, path_out, file_name_out,
                    configFilePath="../IOput/model_config.ini", model_type='lesion')

config.set("smooth", False)

mesh = brain_creation.run(config)
brain_creation.write(mesh, path_out, file_name_out, ['vtk'])

conditioning = 'preconditioned'
poissons = '0,49'
for heterogeneity_model in [hc.Heterogeneity.ONER, hc.Heterogeneity.TWOR, hc.Heterogeneity.FOURR, hc.Heterogeneity.NINER]:
    # heterogeneity_model = hc.Heterogeneity.ONER

    tumor_creator = CreateTumorPRM("./tumor_template_folder/tumor_growth_template.prm")
    tumor_creator.create_materials(mesh, conditioning, poissons, heterogeneity_model)
    tumor_creator.write_materials()
    tumor_creator.complete_prm(path_out, file_name_out, "{}_{}R".format(file_name_out, heterogeneity_model.value))
    output_prm = "/".join([path_out, "{}_tumor_{}R".format(file_name_out, heterogeneity_model.value)])
    tumor_creator.write_and_close_prm(output_prm)

print("COMPLETE")
print("Files written to", path_out)
