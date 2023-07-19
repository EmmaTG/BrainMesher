# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 08:18:00 2023

@author: grife
"""
from abc import ABC, abstractmethod, abstractproperty;
from Mesh import Mesh,Element
from Material_Label import Material_Label

class IWriter(ABC): #Product

    @abstractmethod
    def openWriter(self, filename, path):
        pass
    
    @abstractmethod
    def saveAndClose(self):
        pass
    
    @abstractmethod
    def writeNodes(self, renumber):
        pass
    
    @abstractmethod
    def writeElements(self, renumber):
        pass
    
class BaseWriter():
    
    def __init__(self,ext,tag):
        self.__ext__ = ext;
        self.__tag__ = tag;
        
    def openWriter(self, filename, path):
        self.__path__ = path;
        self.__filename__ = filename        
        if (path[-1] != "\\"):
            path += "\\"
        if len(filename.split('.'))>1:
            filename = filename.split('.')[0]
        filenameOUT = filename + "_" + self.__tag__.upper()
        self.f = open(path + filenameOUT + "." + self.__ext__, 'w')

    def saveAndClose(self):
        self.f.close()
        print("Completed")
        print("New {} file written to ".format(self.__tag__ .upper()) + self.__path__ + self.__filename__)   
    
    
    def initializeMesh(self, mesh):
        self._mesh = mesh;
        
    
    def mesh_statistics(self):
        print("MESH STATISTICS: ")
        print("\tNumber of nodes: " + str(len(self._mesh.nodes)))
        print("\tNumber of elements: " + str(len(self._mesh.elements)))
        if (len(self._mesh.boundaryElements)>0):  
            print("\tNumber of boundary surfaces: " + str(len(self._mesh.boundaryElements)))
    
class ABQWriter(BaseWriter,IWriter):
    
    def __init__(self):
        super().__init__("inp","abq")
        

    def openWriter(self, filename, path):
        super().openWriter(filename, path)
        from datetime import date
        firstLine = "# UCD SCRIPT\n" \
            + "# Inp file created using python script BrainHexMesh\n"\
            + "# Script developed by Emma Griffiths ca. 2022\n"\
            + "# UCD file created on " + date.today().isoformat() + "\n"
        self.f.write(firstLine)

    def saveAndClose(self):
        super().saveAndClose()        
    
    def writeNodes(self,  reNumber):
        nodeMap =  self._mesh.nodes;
        self.f.write("*NODE\n")
        self.oldNumToNewNum = {}
        print("Writing node data")
        for count,n in enumerate(nodeMap.keys()):
            if reNumber:
                nodeNum = count+1
            else:
                nodeNum = n
            self.oldNumToNewNum[n] = nodeNum
            self.f.write(str(nodeNum) + ",\t" + ",\t".join([str(round(i,6)) for i in nodeMap[n]])+"\n")
    
    def writeElements(self, reNumber):
        elementMap = self._mesh.elements;
        self.oldELemTonewELem = {}
        self.f.write("*ELEMENT, TYPE=C3D8, ELSET=ALL\n")
        print("Writing element data")
        for count,e in enumerate(elementMap.keys()):
            assert len(elementMap[e].ica) == 8
            if reNumber:
                elemNum = count+1
            else:
                elemNum = e
            self.oldELemTonewELem[e] = elemNum
            self.f.write(str(elemNum) + ",\t" + ", ".join([str(self.oldNumToNewNum[i]) for i in elementMap[e].ica])+"\n")           
            
        self.writeMaterialsData();
            
    def writeMaterialsData(self):
        elements = self._mesh.elements;
        materialToElements = {}
        for num,element in elements.items():
            materials = element.getMaterial()
            for material in materials:
                if materialToElements.get(material,False):
                    materialToElements[material].append(num)
                else:
                    materialToElements[material] = [num]
        for elsetName,elset_elements in materialToElements.items():
            print("Writing elset data: " + str(elsetName).upper())
            self.f.write("*ELSET, ELSET=" + str(elsetName).upper() + '\n')
            elset_elements.sort()
            for x in range(0,len(elset_elements),15):
                str_to_write = ", ".join([str(self.oldELemTonewELem[y]) for y in elset_elements[x:x+15]])
                self.f.write(str_to_write)
                self.f.write("\n")
        pass
    
class VTKWriter(BaseWriter,IWriter):
    
    def __init__(self):
        super().__init__("vtk","vtk")  

    def openWriter(self, filename, path):
        super().openWriter(filename, path)        
        from datetime import date
        firstLine = "# vtk DataFile Version 2.0\n" \
        + "VTK file created using python script BrainHexMesh script developed by Emma Griffiths"\
        + " ca. 2022 file created on " + date.today().isoformat() \
        + " to view: https://www.paraview.org/download/\n\n" \
        + "ASCII\n" \
        + "DATASET UNSTRUCTURED_GRID\n\n"
        self.f.write(firstLine)  

    def saveAndClose(self):
        super().saveAndClose()
    
    def writeNodes(self, renumber):
        nodeMap = self._mesh.nodes;
        numNodes = len(nodeMap);
        self.f.write(" ".join(["\nPOINTS",str(numNodes),"float"]) + "\n")
        nodeKeys = nodeMap.keys()
        nodeKeys = sorted(nodeKeys)
        count = 0
        self.node_num_map_old_to_new = {}
        for n in nodeKeys:
            nodeNum = count       
            self.node_num_map_old_to_new[n] = nodeNum
            self.f.write(" ".join([str(float(coord)) for coord in nodeMap[n]]) + "\n")
            count += 1
    
    def writeElements(self, renumber):
        elementMap = self._mesh.elements
        boundaryElementMap = self._mesh.boundaryElements;
        # Writing cell data
        numHexElements = len(elementMap)
        numQuadElements = len(boundaryElementMap)
        num_elements = numHexElements + numQuadElements
        self.f.write(" ".join(["\nCELLS",str(num_elements),
                          str(int(numHexElements*9) + int(numQuadElements*5))]) + "\n")
        # element_count = 0
        for element in elementMap.values():
            # element_count += 1
            # elements = list(elements[:4]) + list(elements[4:])        
            # elementNum = element_count
            renumber_ica = []    
            
            for ica_node in element.ica:
                renumber_ica.append(self.node_num_map_old_to_new[ica_node])
            
            self.f.write("8 " + " ".join([str(node) for node in renumber_ica]) + "\n")
        # element_count = 0
        for element in boundaryElementMap.values():
            # element_count += 1
            # elements = list(elements[:4]) + list(elements[4:])        
            # elementNum = element_count
            renumber_ica = []    
            
            for ica_node in element.ica:
                renumber_ica.append(self.node_num_map_old_to_new[ica_node])
            
            self.f.write("4 " + " ".join([str(node) for node in renumber_ica]) + "\n")
            
        #Writing cell types
        self.f.write("\nCELL_TYPES " + str(num_elements) + "\n")
        for cell in range(numHexElements):
            self.f.write("12\n")    
        for cell in range(numQuadElements):
            self.f.write("9\n")
        
        self.writeMaterialsData();
    
    def writeMaterialsData(self):
        elementMap = self._mesh.elements
        boundaryElementMap = self._mesh.boundaryElements;
        numHexElements = len(elementMap)
        numQuadElements = len(boundaryElementMap)
        num_elements = numHexElements + numQuadElements
        
        hexElemKeys = elementMap.keys()
        hexElemKeys = sorted(hexElemKeys)
        #Writing cell data
        self.f.write("\nCELL_DATA " + str(num_elements) + "\n")
        self.f.write("SCALARS material int 1\n")
        self.f.write("LOOKUP_TABLE default\n")
        for num,e in elementMap.items():
            material = e.getMaterial()[0]
            self.f.write(str(int(material)) + "\n")
        for num,e in boundaryElementMap.items():
            material = e.getMaterial()[0]
            self.f.write(str(int(material)) + "\n")            
            
        self.writeData(self._mesh.dataToWrite)
    
    def writeData(self, dataToWrite):
        from re import sub
        elementMap = self._mesh.elements
        #Writing point data
        for d in dataToWrite:
            self.f.write("\nPOINT_DATA " + str(len(elementMap)) + "\n")
            self.f.write("FIELD FieldData 1 \n")
            self.f.write("{} 3 {} float\n".format(d, len(elementMap)))
            for numE, e in elementMap.items():
                data = [0,0,0]
                if not (e.properties.get(d,False)):
                    data = e.properties.get(d)
                    line = sub("[\[\]\(\),]*", '', str(data))
                    self.f.write(line + "\n")

class UCDWriter(BaseWriter,IWriter):
    
    def __init__(self):
        super().__init__("inp","ucd") 

    def openWriter(self, filename, path):
        super().openWriter(filename, path)
        from datetime import date
        firstLine = "*HEADING \n" \
            + "Abaqus file created using python script BrainHexMesh\n"\
            + "Script developed by Emma Griffiths ca. 2022\n"\
            + "INP file created on " + date.today().isoformat() + "\n"\
            + "**\n** Model Definition\n**\n"
        self.f.write(firstLine)

    def saveAndClose(self):
        super().saveAndClose()
    
    def writeElements(self, renumber):
        pass
    
    def writeNodes(self, renumber):
        pass
    
class Writer():    
    
    def openWriter(self,filetype, filename,filePath):
        self.fileType = filetype
        if filetype == "abaqus":
            self.writer = ABQWriter();
        if filetype == "vtk":
            self.writer = VTKWriter();
        if filetype == "ucd":
            self.writer = UCDWriter();
        self.writer.openWriter(filename,filePath)
    
    def initializeWriter(self, mesh):
        if (hasattr(self, "writer")):
            self.writer.initializeMesh(mesh);
        else:
            raise Exception(("Writer has not been opened"));
    
    def writeMeshData(self):
        self.writer.writeNodes(True); 
        self.writer.writeElements(True); 
        self.writer.mesh_statistics()             
    
    def closeWriter(self):
        self.writer.saveAndClose();
        
        

# if __name__ == "__main__":
#     path = "C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\PythonScripts\\BrainMesher"
#     filename = "writer_test"
#     elements = {}
#     elements[1] = Element(1, [269992, 269991, 270077, 270078, 278506, 278505, 278591, 278592], mat=[1])
#     elements[2] = Element(2, [269993, 269992, 270078, 270079, 278507, 278506, 278592, 278593], mat=[1])
#     elements[3] = Element(3, [269994, 269993, 270079, 270080, 278508, 278507, 278593, 278594], mat=[1])
    
#     nodes = {}
#     nodes[269991] = [62.54884263644444, 140.53209928444446, 72.53209928444446]
#     nodes[269992] = [62.396075723333325, 140.40006344666668, 73.99759119666666]
#     nodes[269993] = [62.40417773333333, 140.40485708333335, 75.98512558333333]
#     nodes[269994] = [62.4101889, 140.4100377, 77.9990292]
#     nodes[270077] = [62.39607347333333, 141.99760019666667, 72.40006119666667]
#     nodes[270078] = [61.98393933666667, 141.97769256666666, 73.97768581666666]
#     nodes[270079] = [61.972787649999994, 141.98392084999998, 75.98560610000001]
#     nodes[270080] = [61.98429660000001, 141.9846935, 78.0003105]
#     nodes[278505] = [63.67354880266668, 139.99236300533332, 71.99235760533335]
#     nodes[278506] = [63.58423740400001, 139.59068357999996, 73.98168108]
#     nodes[278507] = [63.588033710000005, 139.57933936999996, 75.99799762]
#     nodes[278508] = [63.59028804, 139.58971149999996, 78.0008271]
#     nodes[278591] =     [63.584229304000004, 141.98168108, 71.59067008]
#     nodes[278592] = [64.0, 142.0, 74.0]
#     nodes[278593] = [64.0, 142.0, 76.0]
#     nodes[278594] = [64.0, 142.0, 78.0]
    
#     mesh = Mesh()
#     mesh.elements = elements
#     mesh.nodes = nodes
    
#     material_labels  = Material_Label()
#     material_labels.addLabelToMap('TestMaterial', 1)
    
#     writer = Writer()
#     writer.openWriter("vtk", filename, path)
#     writer.initializeWriter(mesh)
#     writer.writeMeshData()
#     writer.closeWriter();
    
    
    
    
    
    
    