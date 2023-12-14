from file_converters.Converter import Converter

# path = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/CoreFoamModels"
# path = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/PythonScripts/BrainMesher/IOput/out/"
path = "C:/Users/grife/OneDrive/Documents/PostDoc/Students/Yashasvi/meshes"

filename_in = "csf_brain_centered_VTK.vtk"
filename_out = "csf_brain_centered"

converter = Converter()
converter.convert_file(path, filename_in, "vtk", "ucd", filename_out)