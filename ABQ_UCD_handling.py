# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 08:48:34 2022

@author: grife
"""

def printNodesABQ(nodes,f):
    """
    Prints nodes in style of Abaqus INP file

    Parameters
    ----------
    nodes : Array(Array)
        Array of arrays of node data (node number and coordinates)
    f : open file
        File to write data to.

    Returns
    -------
    None.

    """
    numberOfNodes = len(list(nodes))
    f.write("*Node\n")
    for i in range(numberOfNodes):
        n = nodes[i]
        f.write(str(int(n[0])) + ", " + str(round(n[1],6) )+ "," + str(round(n[2],6)) + "," + str(round(n[3],6))+"\n")
        
def printElementsABQ(elements,f):
    """
    Prints elements of type C3D8 in style of Abaqus INP file

    Parameters
    ----------
    elements : Array(Array)
        Array of arrays of element datat (element number and ica)
    f : open file
        File to write data to.

    Returns
    -------
    None.

    """
    numberOfELements = len(list(elements))
    f.write("**\n*Element, type=C3D8\n")
    for i in range(numberOfELements):
        elem = elements[i]
        string = ""
        for j in elem:
            string += str(j)
            if not j == elem[-1]:
                string += ","
        f.write(string+"\n")
    f.write("**\n")
    

def readABQ(path,filename):
    """
    Function to read in data from Abaqus .inp file and extract teh nodes, elements and elements sets

    Parameters
    ----------
    path : STRING
        full path to location of abaqus gfile
    filename : STRING
        Name of abaqus file w/out extension

    Returns
    -------
    nodeMap : Map(int,array)
        Map of node numbers (keys) to coordinates of nodes (values)
    elementMap : Map(int,array)
        Map of element numbers (keys) to ica (values)
    elementSetsMap : Map(string,array)
        Map of element set names (keys) to element numbers (values).

    """
    
    print("Reading in file " + filename)
    filename = remove_ext(filename)
    f = open(path + filename + ".inp", 'r')
    startedNodes = False
    startedElements = False
    startedElsets= False
    generate= False;
    nodeMap = {}
    elementMap = {}
    elset_name= ""
    elset_elements= []
    elementSetsMap = {}
    for line in f:
        if line.strip()[0:2].lower() == "**":
            pass
        elif line[0:6].lower() == "*elset":
            if line.split(",")[-1].strip().upper() == "GENERATE":
                generate = True
            if not elset_name == "":
                elementSetsMap[elset_name] = elset_elements
            startedElements = False
            startedElsets = True
            elset_name = line.split(",")[1].split("=")[1].strip().replace('"','')
            elset_elements = []
        elif line[0:8].lower() == "*element" and not startedElements:
            startedNodes = False  
            startedElements = True
        elif line[0:5].lower() == "*node" and not startedNodes:
            startedNodes = True
        elif line[0] == "*":
            startedNodes = False
            startedElements = False
            startedElsets = False
            pass
        elif startedElsets:
            if generate:
                [begin,end,interval] = [int(j.strip()) for j in line.strip().split(",")]
                for e in range(begin,end+1,interval):
                    elset_elements.append(int(e))                
                generate = False
            else:
               for e in [j for j in line.strip().split(",")]:
                if not e == '':
                    elset_elements.append(int(e))
        elif startedElements:
            e = [int(i.strip()) for i in line.split(",")]
            elementMap[int(e[0])] = e[1:]
        elif startedNodes:
            n = [i.strip() for i in line.split(",")]
            nodeMap[int(n[0])] = [round(float(k.strip()),6) for k in n[1:]]
    f.close()
    elementSetsMap[elset_name] = elset_elements
    return nodeMap,elementMap,elementSetsMap

def readINP(path,filename):
    from materialMappings import getMaterialMappings
    """
    Function to read in data from UCD file and extract the nodes, elements, elements sets and possible (2D) boundary elements

    Parameters
    ----------
    path : STRING
        full path to location of abaqus file
    filename : STRING
        Name of abaqus file w/out extension

    Returns
    -------
    nodeMap : Map(int,array)
        Map of node numbers (keys) to coordinates of nodes (values)
    elementMap : Map(int,array)
        Map of element numbers (keys) to ica (values)
    materialsMap : Map(string,array)
        Map of element set names (keys) to element numbers (values).
    boundaryElementMap : Map(int,array)
        Map of 2D element numbers (keys) to ica (values)

    """
    
    print("Reading in file " + filename)
    filename = remove_ext(filename)
    f = open(path + filename, 'r')
    [numNodes,numElems,*_] = [ int(x.strip()) for x in f.readline().split(" ")]
    nodeMap = {}
    for n in range(numNodes):
        splitLine =  f.readline().strip().split("  ")
        nodeNumber = int(splitLine[0])
        coords = [float(x.strip()) for x in splitLine[1].split(" ")]
        nodeMap[nodeNumber] = coords
    elementMap = {}
    boundaryElementMap = {}
    materialsMap = {}
    elementMaterialsMap = {}
    materialMappings = getMaterialMappings('number','9R')

    for line in f:
        splitLine = line.strip().split(" ")
        [elemNum,material,elemType,*ica] = list(filter(None, splitLine))
        elemNum = int(elemNum)
        ica = [int(x) for x in ica]
        if elemType == 'hex':
            elementMap[elemNum] = ica
            materialName = materialMappings[material]
            if materialsMap.__contains__(materialName):
                mat_elements = materialsMap[materialName]
            else:
                mat_elements = []
            mat_elements.append(elemNum)
            materialsMap[materialName] = mat_elements
        if elemType == 'quad':
            boundaryElementMap[elemNum] = ica        
            
    f.close()
    return nodeMap, elementMap, boundaryElementMap, materialsMap

def create_element_to_elset_map(elementSetsMap):
    """
    Creates a map of each element to teh elment set/s that it is apart of

    Parameters
    ----------
    elementSetsMap : Map(string,array)
        Map of element set names (keys) to element numbers (values).

    Returns
    -------
    elementToElsetMap : Map(int,string)
        Map of element numbers (keys) to element sets which they are a part of (values).

    """
    elementToElsetMap = {}
    for name,elements in elementSetsMap.items():
        for e in elements:
            if elementToElsetMap.__contains__(e):
                sets = elementToElsetMap[e]
            else:
                sets = []
            sets.append(name)
            elementToElsetMap[e] = sets
    return elementToElsetMap

def writeABQ(path, newfilename, nodeMap, elementMap, elsetsMap={}, nsetMap={}, reNumber = True):
    """
    Write mesh data in Abaqus file style

    Parameters
    ----------
    path : STRING
        full path to location of abaqus file
    newfilename : STRING
        Name of abaqus file w/out extension
    nodeMap : Map(int,array)
        Map of node numbers (keys) to coordinates of nodes (values)
    elementMap : Map(int,array)
        Map of element numbers (keys) to ica (values)
   elsetsMap : Map(string,array), optional
       Map of element set names (keys) to element numbers (values). The default is {}.
    nsetMap : Map(string,array), optional
        Map of node set names (keys) to node numbers (values). The default is {}.
    reNumber : Boolean, optional
        Option to renumber nodes and elements when writing. The default is True.

    Returns
    -------
    None.

    """
    newfilename = remove_ext(newfilename)   
        
    f = open(path + newfilename + ".inp", 'w')
    from datetime import date
    firstLine = "*HEADING \n" \
        + "Abaqus file created using python script BrainHexMesh\n"\
        + "Script developed by Emma Griffiths ca. 2022\n"\
        + "INP file created on " + date.today().isoformat() + "\n"\
        + "**\n** Model Definition\n**\n"
    f.write(firstLine)
    f.write("*NODE\n")
    oldNumToNewNum = {}
    print("Writing node data")
    for count,n in enumerate(nodeMap.keys()):
        if reNumber:
            nodeNum = count+1
        else:
            nodeNum = n
        oldNumToNewNum[n] = nodeNum
        f.write(str(nodeNum) + ",\t" + ",\t".join([str(round(i,6)) for i in nodeMap[n]])+"\n")
    f.write("*ELEMENT, TYPE=C3D8, ELSET=ALL\n")
    oldELemTonewELem = {}
    print("Writing element data")
    for count,e in enumerate(elementMap.keys()):
        if reNumber:
            elemNum = count+1
        else:
            elemNum = e
        oldELemTonewELem[e] = elemNum
        f.write(str(elemNum) + ",\t" + ", ".join([str(oldNumToNewNum[i]) for i in elementMap[e]])+"\n")
        
    
    for nset in nsetMap:
        print("Writing nset data: " + nset.upper())
        f.write("*NSET, NSET=" + nset.upper() + '\n')
        nset_elements = nsetMap[nset]
        for x in range(0,len(nset_elements),15):
            f.write(", ".join([str(oldNumToNewNum[y]) for y in nset_elements[x:x+15]]))
            f.write("\n")
            
    
    for elset in elsetsMap:
        print("Writing elset data: " + elset.upper())
        f.write("*ELSET, ELSET=" + elset.upper() + '\n')
        elset_elements = elsetsMap[elset]
        elset_elements.sort()
        for x in range(0,len(elset_elements),15):
            str_to_write = ", ".join([str(oldELemTonewELem[y]) for y in elset_elements[x:x+15]])
            # print(str_to_write)
            f.write(str_to_write)
            f.write("\n")
    f.close()
    print("Node and element data written to: " + path + newfilename)
    mesh_statistics(elementMap, nodeMap);
    
def mesh_statistics(elementMap, nodeMap, boundaryElements = {}):
    print("MESH STATISTICS: ")
    print("\tNumber of nodes: " + str(len(nodeMap)))
    print("\tNumber of elements: " + str(len(elementMap)))
    if (len(boundaryElements)>0):  
        print("\tNumber of boundary surfaces: " + str(len(boundaryElements)))

def remove_ext(filename):
    if len(filename.split('.'))>1:
        filename = filename.split('.')[0]
    return filename;


def writeElementUCD(f, e, material, elemType, element_ica):    
    f.write(str(e) + "\t" + str(material) + "\t " + elemType + "\t")
    f.write("\t".join([str(n) for n in element_ica])+ "\n")
    
def get_refined_mesh(base_model,length,y,z,local):
    import os.path
    refined_file = '_'.join([str(x) for x in [base_model,length,y,z,local]])+".inp"
    path = "C:\\Users\\grife\\OneDrive\\Documents\\PostDoc\\BrainModels\\MatLab\\INPFiles\\"
    if (os.path.isfile(path + refined_file)):
        print("Refined model found in cache")
        nodeMap, elementMap, elementSetsMap = readABQ(path, refined_file)
        return True,[nodeMap, elementMap, elementSetsMap]
    return False,[{},{},{}]
    
    
    

def writeUCD(path,filenameIN,nodeMap, elementMap, elementToElsetMap = {},
             elset_number_Mappings={},
             boundaryElementMap = {}, boundaryElementToElsetMap = {}, 
             boundary_elset_number_Mappings={'FIXED':'2', 'PRESCRIBED':'1'}, reNumber = True):
    """
    Function to write data to inp file for importing into deal.ii (and viewing on paraview)

    Parameters
    ----------
    path : STRING
        full path to location where file should be saved
    filenameIN : STRING
        Name of UCD file w/out extension
    nodeMap : Map(int,array)
        Map of node numbers (keys) to coordinates of nodes (values)
    elementMap : Map(int,array)
        Map of element numbers (keys) to ica (values)
    elementToElsetMap : Map(int,string), optional
        Map of element numbers (keys) to element sets which they are a part of (values).
        The default is an empty map
    elset_number_Mappings : Map(string,int), optional
        Map of elset names (keys) and the material model numbers to assign for the UCD file (values) 
        The default is is an empty map
    boundaryElementMap : Map(int,array), optional
        Map of 2D bounday element numbers (keys) to ica (values).
        The default is an empty map
    boundary_elset_number_Mappings : Map(string,int), optional
        Map of 2D elset names (keys) and the type of boundary to assign for the UCD file (values) 
        The default is {'FIXED':'2', 'PRESCRIBED':'1'}

    Returns
    -------
    None.

    """
    if (len(elset_number_Mappings) == 0):
        from materialMappings import getMaterialMappings
        elset_number_Mappings = getMaterialMappings("name", "9R")
    from datetime import date
    filenameIN = remove_ext(filenameIN)
    filenameOUT = filenameIN + "_UCD"
    f = open(path + filenameOUT + ".inp", 'w')    
    firstLine = "# UCD SCRIPT\n" \
        + "# Inp file created using python script BrainHexMesh\n"\
        + "# Script developed by Emma Griffiths ca. 2022\n"\
        + "# UCD file created on " + date.today().isoformat() + "\n"
    f.write(firstLine)
    numNodes = len(nodeMap)
    numElements = len(elementMap)+ len(boundaryElementToElsetMap)
    f.write("\t".join([str(numNodes),str(numElements),'0','0','0']) + "\n")
    nodeKeys = nodeMap.keys()
    nodeKeys = sorted(nodeKeys)
    count = 0
    node_num_map_old_to_new = {}
    for n in nodeKeys:
        count += 1
        if reNumber:
            nodeNum = count
        else:
            nodeNum = n        
        node_num_map_old_to_new[n] = nodeNum
        f.write(str(nodeNum) + "\t" + "\t".join([str(node) for node in nodeMap[n]]) + "\n")
        
    elemKeys = elementMap.keys()
    elemKeys = sorted(elemKeys)
    element_count = 0
    for e in elemKeys:        
        element_count += 1
        if not elementToElsetMap.__contains__(e):
            elsetnameList = []
        else:
            elsetnameList = elementToElsetMap[e]               
        if len(elsetnameList) == 0:
            elsetnameList = ['']
        if len(elsetnameList) > 1:
            # print("WARNING: element " + str(e) + "is part of two elsets.")
            nameFound = False;
            for name in elsetnameList:
                if (elset_number_Mappings.__contains__(name)):
                    elsetname = name.replace(" ", '')
                    # print("Element will be categorized at part of elset " + elsetname)
                    nameFound = True;
                    break
            if not nameFound:    
                elsetname = elsetnameList[0].replace(" ", '')
        elsetname = elsetnameList[0].replace(" ", '')
        if elset_number_Mappings.__contains__(elsetname):
            material = elset_number_Mappings[elsetname]
        else:
            material = 0;
        elements = elementMap[e]
        elements = list(elements[:4]) + list(elements[4:])        
        if reNumber:
            elementNum = element_count
        else:
            elementNum = e
        renumber_ica = []    
        
        for ica_node in elements:
            renumber_ica.append(node_num_map_old_to_new[ica_node])
            
        writeElementUCD(f, elementNum, material, 'hex', renumber_ica)
        
    print("Number of boundary Elements: " + str(len(boundaryElementMap)))
    for e in boundaryElementMap:
        if not boundaryElementToElsetMap.__contains__(e):
            elsetnameList = []
        else:
            elsetnameList = boundaryElementToElsetMap[e]               
        if len(elsetnameList) == 0:
            elsetnameList = ['']
        # if len(elsetnameList) > 1:
            # print("WARNING: element " + str(e) + "is part of two elsets.")
            # print("Element will be categorized at part of elset " + elsetnameList[0])
        elsetname = elsetnameList[0].replace(" ", '')
        if boundary_elset_number_Mappings.__contains__(elsetname):
            material = boundary_elset_number_Mappings[elsetname]
        else:
            material = 0;
        writeElementUCD(f, e, material, 'quad', boundaryElementMap[e])
        
    f.close()    
    print("Completed")
    print("New UCD file written to " + path + filenameOUT)
    mesh_statistics(elementMap, nodeMap, boundaryElementMap);


def writeVTK(path, filenameIN, nodeMap, elementMap, elementToMaterial = {},
             elset_number_Mappings={}):
    
    from datetime import date
    filenameIN = remove_ext(filenameIN)
    filenameOUT = filenameIN + "_VTK"
    f = open(path + filenameOUT + ".vtk", 'w')    
    firstLine = "# vtk DataFile Version 2.0\n" \
        + "VTK file created using python script BrainHexMesh script developed by Emma Griffiths"\
        + " ca. 2022 file created on " + date.today().isoformat() \
        + " to view: https://www.paraview.org/download/\n\n" \
        + "ASCII\n" \
        + "DATASET UNSTRUCTURED_GRID\n\n"
    f.write(firstLine)
    # Writing node/point data
    numNodes = len(nodeMap)
    f.write(" ".join(["\nPOINTS",str(numNodes),"float"]) + "\n")
    nodeKeys = nodeMap.keys()
    nodeKeys = sorted(nodeKeys)
    count = 0
    node_num_map_old_to_new = {}
    for n in nodeKeys:
        nodeNum = count       
        node_num_map_old_to_new[n] = nodeNum
        f.write(" ".join([str(float(coord)) for coord in nodeMap[n]]) + "\n")
        count += 1
    
    # Writing cell data
    numElements = len(elementMap)
    f.write(" ".join(["\nCELLS",str(numElements),str(int(numElements*9))]) + "\n")
    elemKeys = elementMap.keys()
    elemKeys = sorted(elemKeys)
    # element_count = 0
    for e in elemKeys:
        # element_count += 1
        elements = elementMap[e]
        # elements = list(elements[:4]) + list(elements[4:])        
        # elementNum = element_count
        renumber_ica = []    
        
        for ica_node in elements:
            renumber_ica.append(node_num_map_old_to_new[ica_node])
        
        f.write("8 " + " ".join([str(node) for node in renumber_ica]) + "\n")
        
    #Writing cell types
    f.write("\nCELL_TYPES " + str(numElements) + "\n")
    for cell in range(numElements):
        f.write("12\n")    
        
    #Writing cell data
    f.write("\nCELL_DATA " + str(numElements) + "\n")
    f.write("SCALARS material int 1\n")
    f.write("LOOKUP_TABLE default\n")
    for e in elemKeys:
        f.write(str(elementToMaterial[e]) + "\n")
        
        
    f.close()
    print("Completed")
    print("New VTK file written to " + path + filenameOUT)
    mesh_statistics(elementMap, nodeMap);


##########################################################################################














############### OLD CODE ######################

    # elif name =="REFINED_ELEMENTS":
    #     for e in elements:
    #         nodes = elementMap[e]
    #         ref_element_map[e] = nodes
    #         ref_element_list.append([e] + nodes)
    #         for n in nodes:
    #             if not ref_node_numbers.count(n):
    #                 ref_node_numbers.append(n)
    #                 ref_node_map[n] = nodeMap[n]
    #                 ref_node_list.append([n]+ nodeMap[n])
        
## Creating isolated refined area model in UCD   
# isolatedUCD = "refined_area_UCD.inp"   
# f = open(path + isolatedUCD, 'w')
# key_for_nodes = ref_node_map
# key_for_elems = ref_element_map
# materialMappings = {}
# firstLine = "# UCD SCRIPT\n" \
#     + "# Inp file created from Abaqus file (" + filenameABQ + ") using python script ABQtoUCDConversion.py\n"\
#     + "# Script developed by Emma Griffiths ca. 2002\n"\
#     + "# UCD file created on " + date.today().isoformat() + "\n"
# f.write(firstLine)
# numNodes = len(ref_node_numbers)
# numElements = len(ref_element_map)
# f.write("\t".join([str(numNodes),str(numElements),'0','0','0']) + "\n")
# nodeKeys = key_for_nodes.keys()
# nodeKeys = sorted(nodeKeys)
# for n in nodeKeys:
#     f.write(str(n) + "\t" + "\t".join([str(node) for node in key_for_nodes[n]]) + "\n")
    
# elemKeys = key_for_elems.keys()
# elemKeys = sorted(elemKeys)
# for e in elemKeys:
#     elsetname = "0"
#     assert len(elsetname) == 1
#     if materialMappings.__contains__(elsetname[0]):
#         material = materialMappings[elsetname[0]]
#     else:
#         material = 0;
#     f.write(str(e) + "\t" + str(material) + "\t hex\t")
#     elements = key_for_elems[e]
#     elements = elements[:4] + elements[4:]
#     f.write("\t".join([str(n) for n in elements])+ "\n")
    
# f.close()    
# print("Completed")
# print("New isolated UCD file written to " + path + isolatedUCD)  
      
 ## Creating isolated refined area model in ABQ
# isolatedABQ = "refined_area_ABQ.inp"
# f = open(path + isolatedABQ, 'w')
# firstline = "** ABAQUS file converted from UCD inp file\n"
# f.write(firstline)
# printNodesABQ(ref_node_list, f)
# printElementsABQ(ref_element_list, f)
# f.close()
# print("Completed")
# print("New isolated ABQ inp file written to " + path + isolatedABQ)

## Creating UCD for full model

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

