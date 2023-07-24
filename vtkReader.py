# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 10:40:35 2023

@author: grife
"""

def readVtk(path,filename):
    fullFilename1 = path + filename + ".vtk"
    filenameData = {}
    if (exists(fullFilename1)):

        print("##Reading: " + filename)
        
        # Set up poly data reader for result set 1
        reader = vtkUnstructuredGridReader()
        reader.SetFileName(fullFilename1)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        return reader.GetOutput()
        
    return None

from vtk import vtkIdList,vtkUnstructuredGridReader
from os.path import exists 
from Mesh import Node, Element, Mesh;
from Writer import Writer;

inputPath = "C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\Tumor_growth\\"
f = "output_test_slice_x"
dataMap = {}
grid = readVtk(inputPath, f);

mesh = Mesh();
mesh.dataToWrite = ["displacement"]
displacementArray = grid.GetPointData().GetArray("displacement")
materialsArray = grid.GetPointData().GetArray("material_ids")

disp_data = {} # Key = starting position, displacement
point_position_to_nodes = {} # Key= point position, value  = node value
cell_point_to_nodes = {} # Key= cell pointId, value  = node value
nodeMap = {}
elementsMap = {}
elementToMaterial = {}
for k in range(grid.GetNumberOfCells()):
    cellIds = vtkIdList() 
    grid.GetCellPoints(k, cellIds)
    ica = []
    mat = 0
    for i in range(8):
        pointId = cellIds.GetId(i)
        pointPosition = grid.GetPoint(pointId)        
        mat += materialsArray.GetValue(pointId)
        position_key = "-".join([str(x) for x in pointPosition])
        nodeValue = point_position_to_nodes.get(position_key,-1)
        if nodeValue != -1:
            cell_point_to_nodes[pointId] = nodeValue
            assert [disp_data.get(nodeValue,-1) != -1]
        else:
            nodeValue = pointId
            cell_point_to_nodes[pointId] = nodeValue
            point_position_to_nodes[position_key] = pointId
            displacement = displacementArray.GetTuple(pointId)
            disp_data[pointId] = displacement
            newPosition = [0,0,0]
            for d in range(len(displacement)):
                newPosition[d] = pointPosition[d]+displacement[d]
            node = Node(int(nodeValue),newPosition)
            node.addData("displacement", [round(y,6) for y in displacement])
            nodeMap[int(nodeValue)] = node;
        ica.append(nodeValue)
    element = Element(k,ica)
    element.setMaterial(int(mat/8))
    elementsMap[k] = element
    
mesh.nodes = nodeMap
mesh.elements = elementsMap

## Write deformed mesh to new vtk
path = inputPath
filename = "output_"+f
writer = Writer()
writer.openWriter("vtk", filename, path, mesh)
writer.writeMeshData()
writer.closeWriter()







# from vtk.util.numpy_support import vtk_to_numpy
# from vtkmodules.vtkCommonDataModel import (
#     vtkCellTypes,
#     vtkGenericCell
# )
# from ABQ_UCD_handling import writeVTK


# Read points.
# vtk_points = grid.GetPoints()
# xyz3d = vtk_to_numpy( vtk_points.GetData() )

# Read cells.
# cell_locations = grid.GetCellLocationsArray()
# cells = grid.GetCells().GetData()
# pointData = grid.GetPointData()
# numArraysInData = pointData.GetNumberOfArrays()
# for i in range(numArraysInData):
#     arr = grid.GetPointData().GetArray(i)
#     print(arr.GetName())
# von_mises = grid.GetPointData().GetArray(1)

# numPoints = grid.GetNumberOfPoints()
# avg_min_strain = np.zeros(numPoints)
        # if disp_data.get(pointId,False):
        #     print("a")
        # else:
        #     disp_data[pointId] = displacement.GetTuple(pointId)
        # avgValue += min_strain.GetValue(pointId)
        # print(grid.GetPoint(cellIds.GetId(i)))   
        # print(pointId)  
        # print("min strain: ",min_strain.GetValue(pointId))
        # print("von mises: ",von_mises.GetValue(pointId))
    # for i in range(8):
    #     pointId = cellIds.GetId(i)
        # avg_min_strain[pointId] = avgValue/8.;
        
# Create new unstrauctured grid
# output = vtk.vtkUnstructuredGridReader()
# output.SetFileName("C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\Data\\Results\\template.vtk")
# output.ReadAllVectorsOn()
# output.ReadAllScalarsOn()
# output.Update()
# output_data = grid
# output_pointData = output_data.GetPointData()

# newArray = vtk.vtkFloatArray();
# newArray.SetName("averaged_min_strain")
# newArray.SetNumberOfValues(numPoints)
# for x in range(numPoints):
#     newArray.SetValue(x,avg_min_strain[x])
# output_pointData.AddArray(newArray)

# writer = vtk.vtkUnstructuredGridWriter()
# writer.SetInputData(output_data)
# outputFullFilename = inputPath + "incision_full_phere_averaged.vtk"
# writer.SetFileName(outputFullFilename)
# writer.Update()
# writer.Write()

    # print(grid.GetPointData(cellIds.GetId(i)))
# for idx in range(0,numArraysInData):
#     arr_1 = pointData.GetArray(idx)
#     print(arr_1)

# it = grid.NewCellIterator()
# it.InitTraversal()
# count = 0
# while ((not it.IsDoneWithTraversal()) and (count <100)):
#     count += 1
#     cell = vtkGenericCell()
#     it.GetCell(cell)
#     id_list = vtk.vtkIdList() 
#     print(cell.GetPointIds())
    # for i in range(8):
        # print(cell.GetPoints())
    
#     cellName = vtkCellTypes.GetClassNameFromTypeId(cell.GetCellType())
#     # print(cellName, 'NumberOfPoints:', cell.GetNumberOfPoints(), 'CellDimension:', cell.GetCellDimension())
#     # legendValues.InsertNextValue(cellName)
#     it.GoToNextCell()