# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 11:10:39 2022

@author: grife

Creates the .prm file to be used in EFI1 simulations

Parameters
----------
fileOutput : string
    Name of output file 
sim_name : string
    Name of simulation for which file is to be created
    
Returns
-------
A .prm file
"""

import os

def readConfigFile(path,simulationName):
    radius = 0;
    center = [0,0,0]
    with open("/".join([path, simulationName]) + ".txt") as configReader:
        line = configReader.readline();
        while line != '':
            if (line[:22] == "Concentration radius: "):
                radius = int(float(line[22:].strip()))
            elif (line[:8] == "Center: "):
                center = [int(float(x.strip())) for x in line[8:].split(", ")]
            line = configReader.readline();
    return [center, radius]
        

# simulations = ["OAS1_0002_MR1","OAS1_0004_MR1","OAS1_0005_MR1","OAS1_0006_MR1","OAS1_0007_MR1","OAS1_0009_MR1","OAS1_0011_MR1"];
# for s in simulations:
def createPRM(configPath, pathout, simulationName, path_to_template, output_template_name=None):
    if output_template_name is None:
        tmp_name = path_to_template.split("/")[-1].split(".")
        full_name = path_to_template if len(tmp_name)==1 else "".join(tmp_name[:-1])
        output_template_name = full_name + "_Complete"
    [center,radius] = readConfigFile(configPath, simulationName)
    replacement_values = {}
    
    # Input filearameters
    testing_device = "atrophy_1"
    replacement_values["%testing_device%"] = testing_device
    
    replacement_values["%center%"] = str(center).strip('[').strip(']')
    replacement_values["%radius%"] = radius
    
    inp_file = simulationName + "_UCD.inp"
    replacement_values["%inp_file%"] = inp_file

    
    output_file = simulationName
    replacement_values["%output%"] = output_file

    filepathToTemplate = path_to_template

    template = open(filepathToTemplate, 'r')
    inp_tmp = template.read()
    template.close()    
                
    new_inp = inp_tmp
    assert new_inp.replace("%","") != new_inp, ".inp file already fully filled in"
    for key,value in replacement_values.items():
        new_inp = new_inp.replace(key, str(value))            
    assert new_inp.replace("%","") == new_inp, ".inp file not fully filled in"
    output_filename = pathout + "/" + output_template_name  + ".prm"
    outputFile = open(output_filename,'w')
    print(output_file + " WRITTEN")
    outputFile.write(new_inp)
    outputFile.close()

    return output_filename

if __name__ == "__main__":
    configPath = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/PythonScripts/BrainMesher/Atrophy/OAS1_0004_MR1"
    pathout = "C:/Users/grife/OneDrive/Documents/PostDoc/BrainModels/PythonScripts/BrainMesher/Atrophy/OAS1_0004_MR1"
    filename = "OAS1_0004_MR1"
    createPRM(configPath, pathout, filename)


                                                                                    
        
        

    

    


