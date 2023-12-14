from readers.Reader import Reader
from writers.Writers import Writer
import mesh.mesh_transformations as mt


def round_to_interval(value, interval):
    return interval*round(value/interval)

path = 'C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/CoreFoamModels/'
filenames = ["ellipse_in_cube_4mm_30mm_6mm_curve11"]
filenames = ["ellipse_curve11_4mm"]
for filename in filenames:
    # fileOut = filename
    fileOut = filename + "_fixed"

    reader = Reader('abaqus')
    reader.openReader(filename,path)
    mesh = reader.getMesh()

    # if "curve4" in filename or "curve10" in filename:
    mt.rotate_mesh(mesh.nodes,2,90)
    mt.rotate_mesh(mesh.nodes,0,90)
    # elif "curve5" in filename or "curve9" in filename or "curve6" in filename:
    #     mt.rotate_mesh(mesh.nodes,2,90)

    for num,node in mesh.nodes.items():
        coords = node.getCoords()
        if (coords[1] == -4.0) or (coords[1] == 4.0):
            coords[0] = round_to_interval(coords[0], 2 / 3.)
            coords[2] = round_to_interval(coords[2], 2 / 3.)

        if (coords[2] == -4.0) or (coords[2] == 4.0):
            coords[1] = round_to_interval(coords[1], 2 / 3.)
            coords[0] = round_to_interval(coords[0], 2 / 3.)

        if (coords[0] == -4.0):
            coords[2] = round_to_interval(coords[2], 2 / 3.)
            coords[1] = round_to_interval(coords[1], 2 / 3.)

    writer = Writer()
    writer.openWriter('vtk', fileOut, path)
    writer.writeMeshData(mesh)
    writer.closeWriter()

    writer.openWriter('abaqus', fileOut, path)
    writer.writeMeshData(mesh)
    writer.closeWriter()