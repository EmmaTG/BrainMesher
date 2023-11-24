
from readers.Reader import ABQReader, VTKReader
from writers.Writer import Writer


class Converter:

    @staticmethod
    def convert_file(path, file1, convert_from, convert_to, filename_out=''):
        split_file1 = file1.split(".")
        ext = ""
        if len(split_file1) < 2:
            if convert_from == 'abaqus':
                ext = "inp"
            elif convert_from == 'vtk':
                ext = "vtk"
        else:
            ext = "".join(split_file1[1].split())

        filename = "".join(split_file1[0].split())
        if filename_out == '':
            filename_out = filename

        if ext == "inp" and convert_from == "abaqus":
            reader = ABQReader()
        elif ext == "vtk" and convert_from == "vtk":
            reader = VTKReader()
        else:
            print("File type {} not support".format(convert_from))
            return

        reader.openReader(filename, path)
        mesh = reader.getMesh()

        writer = Writer()
        new_filename = "_".join([filename_out,"converted"])
        if convert_to == 'abaqus':
            writer.openWriter('abaqus', new_filename, path)
        elif convert_to == 'vtk':
            writer.openWriter('vtk', new_filename, path)
        elif convert_to == 'ucd':
            writer.openWriter('ucd', new_filename, path)

        writer.writeMeshData(mesh)
        writer.closeWriter()

