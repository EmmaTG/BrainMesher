import brain_creation
from config.Config import ConfigFile
import writers.HeterogeneityConverter as hc
from readers.Reader import Reader
from personal_scripts.create_prms import CreateAtrophyPRM, CreateTumorPRM

# Model type options: basic_fullcsf, basic_partilacsf, basic_nocsf, atrophy, lesion

# ####################################################################################################################
# # ####### DEBUG CODE #######
# path_in = "../IOput/in"
# file_name_in = "aseg.mgz"
# path_out = "../IOput/out"
# file_name_out = "21_region_brain_VTK"
#
# reader = Reader('vtk')
# reader.openReader(file_name_out, path_out)
# mesh = reader.getMesh()
# reader.closeReader()
#
# conditioning = 'preconditioned'
# poissons = '0,49'
# for heterogeneity_model in [hc.Heterogeneity.NINER, hc.Heterogeneity.NINETEENR]:
#     # heterogeneity_model = hc.Heterogeneity.ONER
#
#     atrophy_creator = CreateAtrophyPRM("./atrophy_template_folder/atrophy_template_V2.prm")
#     atrophy_creator.create_materials(mesh, conditioning, poissons, heterogeneity_model)
#     atrophy_creator.write_materials()
#     output_prm = "/".join([path_out, "{}_atrohpy_{}R".format(file_name_out, heterogeneity_model.value)])
#     # atrophy_creator.write_prm(output_prm)
#     atrophy_creator.close_prm()
#
# # mesh.write(path_out, file_name_out, ['vtk'])
# print("Compete")
####################################################################################################################

# ####### ATROPHY CODE #######
# Preferences are defined in ConfigFile
path_to_oasis = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/OASIS/OASIS-1/oasis_cs_freesurfer_disc1.tar/oasis_cs_freesurfer_disc1/disc1"
path_in = "../IOput/in"
file_name_in = "aparc_DKTatlas+aseg.mgz"
path_to_out = "../IOput/out/atrophy_files"
# filenames = ["Silvia_brain"]
filenames = ["OAS1_0002_MR1"]
# filenames = ["OAS1_0002_MR1", "OAS1_0004_MR1", "OAS1_0005_MR1", "OAS1_0006_MR1", "OAS1_0007_MR1", "OAS1_0009_MR1"]
for name in filenames:
    path_in = "/".join([path_to_oasis, name, "mri"])
    path_out = "/".join([path_to_out, name])
    file_name_out = name

    config = ConfigFile(path_in, file_name_in, path_out, file_name_out,
                        configFilePath="../IOput/model_config.ini", model_type='atrophy')
    mesh = brain_creation.run(config)
    mesh.write(path_out, file_name_out, ['vtk','ucd'])

    # reader = Reader('vtk')
    # reader.openReader(file_name_out + "_VTK", path_out)
    # mesh = reader.getMesh()
    # reader.closeReader()

    conditioning = 'preconditioned'
    poissons = '0,49'
    for heterogeneity_model in [hc.Heterogeneity.ONER, hc.Heterogeneity.TWOR, hc.Heterogeneity.FOURR,
                                hc.Heterogeneity.NINER]:
    # for heterogeneity_model in [hc.Heterogeneity.NINETEENR]:

        atrophy_creator = CreateAtrophyPRM("./atrophy_template_folder/atrophy_template_V2.prm")
        atrophy_creator.create_materials(mesh, conditioning, poissons, heterogeneity_model)
        atrophy_creator.write_materials()
        atrophy_creator.complete_prm(path_out, file_name_out, "{}_atrohpy_{}R".format(file_name_out, heterogeneity_model.value))
        output_prm = "/".join([path_out, "{}_atrohpy_{}R".format(file_name_out, heterogeneity_model.value)])
        atrophy_creator.write_prm(output_prm)
        atrophy_creator.close_prm()

    print("COMPLETE")
    print("Files written to",path_out)

####################################################################################################################
# ####### TUMOR LESION CODE #######
# path_in = "../IOput/in"
# file_name_in = "aparc_DKTatlas+aseg.mgz"
# path_to_out = "../IOput/out/tumor_files"
# filenames = ["Silvia_brain"]
# # filenames = ["OAS1_0002_MR1"]
# # filenames = ["OAS1_0002_MR1", "OAS1_0004_MR1", "OAS1_0005_MR1", "OAS1_0006_MR1", "OAS1_0007_MR1", "OAS1_0009_MR1"]
# for name in filenames:
#     path_out = "/".join([path_to_out, name])
#     file_name_out = name
#
#     # config = ConfigFile(path_in, file_name_in, path_out, file_name_out,
#     #                     configFilePath="../IOput/model_config.ini", model_type='lesion')
#     # mesh = brain_creation.run(config)
#     # mesh.write(path_out, file_name_out, ['vtk'])
#
#     reader = Reader('vtk')
#     reader.openReader(file_name_out + "_VTK", path_out)
#     mesh = reader.getMesh()
#     reader.closeReader()
#
#
#     conditioning = 'preconditioned'
#     poissons = '0,49'
#     for heterogeneity_model in [hc.Heterogeneity.ONER, hc.Heterogeneity.TWOR, hc.Heterogeneity.FOURR,
#                                 hc.Heterogeneity.NINER]:
#         # heterogeneity_model = hc.Heterogeneity.ONER
#
#         tumor_creator = CreateTumorPRM("./tumor_template_folder/tumor_growth_template.prm")
#         tumor_creator.create_materials(mesh, conditioning, poissons, heterogeneity_model)
#         tumor_creator.write_materials()
#         tumor_creator.complete_prm(path_out, file_name_out, "{}_tumor_{}R".format(file_name_out, heterogeneity_model.value))
#         output_prm = "/".join([path_out, "{}_tumor_{}R".format(file_name_out, heterogeneity_model.value)])
#         tumor_creator.write_prm(output_prm)
#         tumor_creator.close_prm()
#
#     print("COMPLETE")
#     print("Files written to", path_out)
