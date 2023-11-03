
from abc import ABC, abstractmethod

from mesh.Mesh import Mesh
from mesh.Node import Node
from mesh.Element import HexElement, QuadElement
from readers.vtkReader import readVtk
from vtk import vtkIdList

class IReader(ABC):
    """
    An abstract class used to define methods for mesh writers

    Abstract Methods
    -------
    openWriter(filename, path)
        open a file using the given path and filename
    saveAndClose()
        save and close the file
    writeNodes(renumber)
        write nodes and node data to file
    writeElements(renumber)
        write element and element data to file
    """

    @abstractmethod
    def openReader(self, filename, path):
        pass

    @abstractmethod
    def closeReader(self):
        pass

    @abstractmethod
    def readNodes(self):
        pass

    @abstractmethod
    def readElements(self):
        pass

    @abstractmethod
    def getMesh(self):
        pass


class BaseReader:
    """
    A base class used to define methods for mesh writers

    Methods
    -------
    openWriter(filename, path)
        open a file using the given path and filename
    saveAndClose()
        save and close the file
    initializeMesh(mesh)
        Initialize mesh object to be written
    mesh_statistics(renumber)
        print out data for the mesh written
    """

    def __init__(self, ext, tag):
        """
        Parameters
        ----------
        ext : str
            The extension of the file to be written
        tag : str
            The unique tag to identify the type of writer
        """
        self.__ext__ = ext
        self.__tag__ = tag
        self.__filename__ = ''
        self.mesh = Mesh()

    def openReader(self, filename, path):
        """Open file to be written to given path and filename.
        If the path does not end if separtor characters thes are added.
        If file name contains an extension, this is removed and replaced the
        tag of the writer '__tag__' and the correct extension '__ext__'

        Parameters
        ----------
        filename : string
            filename of writer
        path : string
            path to filename
        """
        if (path[-1] != "/"):
            path += "/"
        self.__path__ = path

        self.__filename__ = filename

        self.f = open(path + filename + "." + self.__ext__, 'r')

    def closeReader(self):
        """
        Saves and closed file.
        """
        self.f.close()
        print("Completed")


        return self.mesh


class ABQReader(BaseReader, IReader):
    """
    A class to write data in a format suitable for Abaqus.
    This class inherits methos from the Base writer class.
    This class implents the Iwriter interface

    Methods
    -------
    openWriter(filename, path)
        open a file using the given path and filename
    saveAndClose()
        save and close the file
    """

    def __init__(self):
        super().__init__("inp", "abq")
        self.node_lines = []
        self.element_lines = []
        self.elset_lines = {}

    def getMesh(self):
        self.readNodes()
        self.readElements()
        super().closeReader()
        return self.mesh

    def openReader(self, filename, path):
        """Open file to be written to given path and filename as defined
        by the super class BaseWriter. The first line of this
        file is also written

        Parameters
        ----------
        filename : string
            filename of writer
        path : string
            path to filename
        """
        super().openReader(filename, path)
        start_nodes = False
        start_elements = False
        start_elsets = False
        elset_name = None
        for line in self.f:
            line = "".join(line.split()).upper()
            if line != '':
                if "*NODE" in line:
                    start_nodes = True
                    start_elements = False
                    start_elsets = False

                elif "*ELEMENT" in line:
                    start_nodes = False
                    start_elements = True
                    start_elsets = False

                elif "*ELSET" in line:
                    start_nodes = False
                    start_elements = False
                    start_elsets = True
                    splitLine = line.split(",")
                    elset_name = splitLine[1].split("=")[1]
                    assert not self.elset_lines.get(elset_name,False), " elset {} already exists".format(elset_name)
                    self.elset_lines[elset_name] = []

                elif start_nodes:
                    self.node_lines.append(line)

                elif start_elements:
                    self.element_lines.append(line)

                elif start_elsets:
                    assert elset_name is not None
                    self.elset_lines[elset_name].append(line)

    def readNodes(self):
        node_map = {}
        for line in self.node_lines:
            split_line = line.split(",")
            node_number = int(split_line[0].strip())
            coords = [float(x) for x in split_line[1:]]
            new_node = Node(node_number, coords)
            node_map[node_number] = new_node
        self.mesh.nodes = node_map

    def readElements(self):
        element_map = {}
        boundary_elements = {}
        for line in self.element_lines:
            split_line = line.split(",")
            element_number = int(split_line[0].strip())
            ica = [int(x) for x in split_line[1:]]
            node_ica = []
            for n in ica:
                node_ica.append(self.mesh.nodes.get(n))
            if len(ica) == 8:
                new_element = HexElement(element_number, node_ica)
                element_map[element_number] = new_element
            elif len(ica) == 4:
                new_element = QuadElement(element_number, node_ica)
                boundary_elements[element_number] = new_element

        self.mesh.elements = element_map
        self.mesh.addBoundaryElements(boundary_elements)

        self.getMaterialsData()

    def getMaterialsData(self):
        for name, lines in self.elset_lines.items():
            mat = int(name.replace("E", ""))
            for line in lines:
                elements = [int(x.strip()) for x in line.split(",")]
                for e in elements:
                    element = self.mesh.elements.get(e, None)
                    if element is not None:
                        element.addMaterial(mat)

class VTKReader(BaseReader, IReader):

    def __init__(self):
        super().__init__("vtk", "vtk")
        self.grid = None

    def getMesh(self):
        self.readNodes()
        self.readElements()
        return self.mesh

    def openReader(self, filename, path):
        if (path[-1] != "/"):
            path += "/"
        self.grid = readVtk(path, filename)

    def readNodes(self):
        """
        Write nodes in the vtk format:
            coord1 coord2 ... coordn
        Parameters
        ----------
        renumber : boolean
            deteremines whetheer the node numbers shoudl be renumbered or not.
        """
        num_points = self.grid.GetNumberOfPoints()
        node_map = {}
        for x in range(num_points):
            point = self.grid.GetPoint(x)
            new_node = Node(x, point)
            node_map[x] = new_node
        self.mesh.nodes = node_map


    def readElements(self):
        element_map = {}
        boundary_elements = {}
        for k in range(self.grid.GetNumberOfCells()):
            cellIds = vtkIdList()
            self.grid.GetCellPoints(k, cellIds)
            node_ica = []
            for i in range(cellIds.GetNumberOfIds()):
                pointId = cellIds.GetId(i)
                node_ica.append(self.mesh.nodes[pointId])
            if cellIds.GetNumberOfIds() == 8:
                new_element = HexElement(k,node_ica)
                element_map[k] = new_element
            elif cellIds.GetNumberOfIds() == 4:
                new_element = QuadElement(k,node_ica)
                boundary_elements[k] = new_element
        self.mesh.elements = element_map
        self.mesh.addBoundaryElements(boundary_elements)

        self.readMaterialsData()

    def readMaterialsData(self):
        cell_data = self.grid.GetCellData()
        number_arrays =cell_data.GetNumberOfArrays()
        for i in range(number_arrays):
            arr = cell_data.GetArray(i)
            arr_name = arr.GetName()
            if "material" in arr_name:
                for e in range(arr.GetNumberOfValues()):
                    mat = arr.GetValue(e)
                    element = self.mesh.elements.get(e, None)
                    if element is not None:
                        element.addMaterial(mat)
            else:
                for e in range(arr.GetNumberOfValues()):
                    data = arr.GetValue(e)
                    element = self.mesh.elements.get(e, None)
                    if element is not None:
                        element.addProperty(arr_name, data)
                self.mesh.cellData.append(arr_name)
        self.readPointData()

    def readPointData(self):
        cell_data = self.grid.GetPointData()
        number_arrays = cell_data.GetNumberOfArrays()
        for i in range(number_arrays):
            arr = cell_data.GetArray(i)
            arr_name = arr.GetName()
            for n in range(arr.GetNumberOfValues()):
                data = arr.GetValue(n)
                node = self.mesh.nodes.get(n, None)
                if node is not None:
                    node.addData(arr_name, data)
            self.mesh.dataToWrite.append(arr_name)
