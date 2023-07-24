# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 13:55:28 2023

@author: grife
"""

class ConfigFile():
    def __init__(self, inputPath, inputFileName):
        self.readFile = True
        self.fileInPath = inputPath
        self.fileIn = inputFileName
        self.readData = False
        self.data = []        
        self.fileoutPath = inputPath
        self.writeToFile = True
        self.fileout = inputFileName
        self.fileoutTypes = ['vtk'] # 'ucd' | 'vtk' | 'abaqus'
        self.Coarsen = True
        self.Add_CSF = True
        self.layers = 1
        self.Smooth = True
        self.iterations = 6
        self.coeffs = [0.6,-0.4]
        
        self.Smooth_regions = ['Lesion']
        self.region_iterations = [4]
        self.region_coeffs =[[0.6,-0.4]]
        
        self.material_labels  = Material_Label()
        self.material_labels.addLabelToMap('BrainStem', 16)
        self.material_labels.addLabelToMap('GreyMatter', [3,42]) # Left, Right
        self.material_labels.addLabelToMap('WhiteMatter' , [2,41,77]); # Left, Right, WM-hypointensities
        self.material_labels.addLabelToMap('Corpuscallosum' , [251,252,253,254,255]); # CC_Posterior, CC_Mid_Posterior, CC_Central, CC_Mid_Anterior, CC_Anterior
        self.material_labels.addLabelToMap('BasalGanglia' , [11,50,12,51,13,52,26,58,62,30]); # Caudate(L&R), Putamen(L&R), Palladium(L&R), Accumbens Area(L&R), vessel(L&R)
        self.material_labels.addLabelToMap('Cerebellum' , [7,46,8,47]); # WM(L&R), GM(L&R)
        self.material_labels.addLabelToMap('Thalamus' , [10,49,28,60]); # Thalamus(L&R), Ventral DC(L&R)
        self.material_labels.addLabelToMap('Hippocampus' , [17,53]); # Left, Right
        self.material_labels.addLabelToMap('Amygdala' , [18,54]); # Left, Right
        self.material_labels.addLabelToMap('Lesion' , [25,57]); # Left, Right
        self.material_labels.addLabelToMap('CSF' , [24]); # Left, Right, Lateral(L&R), 3rd, 4th ventricles
        self.material_labels.addLabelToMap('Unused' , [4,43,14,15,31,63,85]); # Lateral(L&R), 3rd, 4th NOTE: These labels will be removed to create empty spaces!!
        # OR
        # material_labels.addLabelToMap('Left-Lateral-Ventricle' , [4]);
        # material_labels.addLabelToMap('Right-Lateral-Ventricle' , [43]);
        # material_labels.addLabelToMap('Left-Inf-Lat-Vent' , [5]);
        # material_labels.addLabelToMap('Right-Inf-Lat-Vent ' , [44]);
        # material_labels.addLabelToMap('3rd-Ventricle' , [14]);
        # material_labels.addLabelToMap('4th-Ventricle' , [15]);
        #
        # material_labels.addLabelToMap('Left-choroid-plexus' , [31]);
        # material_labels.addLabelToMap('Right-choroid-plexus' , [63]);
        # material_labels.addLabelToMap('Optic-Chiasm' , [85]);
        
    def add_data(self, importedData):
        self.readFile = False
        self.readData = True
        self.data = importedData
        
    def openConfigFile(self):
        self.f = open(self.fileoutPath + self.fileout + ".txt", 'w')
        if self.readFile:
            self.f.write("Input file: " + self.fileInPath + self.fileIn + "\n")
        else:
            self.f.write("Data input fromn user")
        self.f.write("Write output to file: " + str(self.writeToFile) + "\n")
        if self.writeToFile:            
            self.f.write("Output written to: " + self.fileoutPath + self.fileout + "\n")
            self.f.write("Output file types: " + ", ".join(self.fileoutTypes) + "\n")
        
        self.f.write("Coarsen: " + str(self.Coarsen) + "\n")
        self.f.write("Add CSF: " + str(self.Add_CSF) + "\n")
        if (self.Add_CSF):
            self.f.write("Layers of CSF: " + str(self.layers) + "\n")
        self.f.write("Smooth global mesh: " + str(self.Smooth) + "\n")
        if (self.Smooth):            
            self.f.write("Iterations: " + str(self.iterations) + ", coeffs: " + ", ".join([str(x) for x in self.coeffs]) + "\n")
        for count,r in enumerate(self.Smooth_regions):
            self.f.write("Smooth region: " + r + "\n")
            self.f.write("Iterations: " + str(self.region_iterations[count]) + ", coeffs: " + \
                         ", ".join([str(x) for x in self.region_coeffs[count]]) + "\n")
        
    def writeToConfig(self, name, value):
        self.f.write("{}: {}\n".format(name, str(value)))
        
    def closeConfigFile(self):         
        self.f.write("COMPLETE\n")
        self.f.close();