# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 09:28:06 2023

@author: grife
"""

import vtk;
from os.path import exists
from Mesh import Node, Element, Mesh;

def readVtk(path,filename):
    fullFilename1 = path + filename + ".vtk"
    filenameData = {}
    if (exists(fullFilename1)):

        print("##Reading: " + filename)
        
        # Set up poly data reader for result set 1
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(fullFilename1)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        return reader.GetOutput()
        
    return None

inputPath = "C:\\Users\grife\OneDrive\Documents\PostDoc\BrainModels\Tumor_growth\\"
f = "OAS1_0004_MR1"
grid = readVtk(inputPath, f);
displacementArray = grid.GetPointData().GetArray("displacement")
materialsArray = grid.GetPointData().GetArray("material_ids")


## Get mesh data
mesh = Mesh();
disp_data = {} # Key = starting position, displacement
point_position_to_nodes = {} # Key= point position, value  = node value
cell_point_to_nodes = {} # Key= cell pointId, value  = node value
nodeMap = {}
elementsMap = {}
elementToMaterial = {}
for k in range(grid.GetNumberOfCells()):
    cellIds = vtk.vtkIdList() 
    grid.GetCellPoints(k, cellIds)
    ica = []
    mat = 0
    displacement_tot = [0,0,0];
    for i in range(8):
        pointId = cellIds.GetId(i)
        pointPosition = grid.GetPoint(pointId)        
        mat += materialsArray.GetValue(pointId)
        position_key = "-".join([str(round(x,6)) for x in pointPosition])
        nodeValue = point_position_to_nodes.get(position_key,-1)
        if nodeValue != -1:
            cell_point_to_nodes[pointId] = nodeValue
            node = nodeMap[int(nodeValue)]
            displacement = disp_data[nodeValue];
            assert [disp_data.get(nodeValue,-1) != -1]
        else:
            nodeValue = pointId
            cell_point_to_nodes[pointId] = nodeValue
            point_position_to_nodes[position_key] = pointId
            displacement = [round(y,6) for y in displacementArray.GetTuple(pointId)]
            disp_data[pointId] = displacement
            newPosition = [round(y,6) for y in pointPosition]
            for d in range(len(displacement)):
                newPosition[d] = pointPosition[d]+displacement[d]
            node = Node(int(nodeValue),newPosition)
            node.addData("displacement", displacement)
            nodeMap[int(nodeValue)] = node;
        for d in range(3):
            displacement_tot[d] += displacement[d];
        ica.append(node)
    element = Element(k,ica)
    element.setMaterial((mat/8))
    element.properties['displacement'] = displacement_tot
    elementsMap[k] = element
    
mesh.nodes = nodeMap
mesh.elements = elementsMap

my_vtk_dataset = vtk.vtkUnstructuredGrid();

## Create Points
points = vtk.vtkPoints()
num_Points = len(list(nodeMap.keys()))
for nId in range(num_Points):
    node = nodeMap[nId];
    points.InsertPoint(node.number, node.getCoords())
my_vtk_dataset.SetPoints(points)

## Connectivity
my_vtk_dataset.Allocate(len(elementsMap))
for element in elementsMap.values():
    point_ids = [n.number for n in element.ica]
    my_vtk_dataset.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, point_ids)

##  Create arrays
# Add displacement data 
array = vtk.vtkDoubleArray()
array.SetNumberOfComponents(3)
array.SetNumberOfTuples(len(nodeMap))
array.SetName('My Array')
for node in nodeMap.values():
    values = node.getData('displacement');
    array.SetTuple(node.number, values)
my_vtk_dataset.GetPointData().AddArray(array)

writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName("output.vtu")
writer.SetInputData(my_vtk_dataset)
writer.Write()