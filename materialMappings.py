# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 10:38:33 2023

@author: grife
"""
def getMaterialMappings(start, model):
    if (start =="name" and model == "9R"):
        return {'AMYGDALA':'18', 'BASALGANGLIA':'11', 'BRAINSTEM':'16', 
                    'CEREBELLUM':'7', 'CORPUSCALLOSUM':'251', 'GREYMATTER':'3',
                    'HIPPOCAMPUS':'17', 'THALAMUS':'10', 'WHITEMATTER':'2', 'CSF':'24'
                    }
    if (start =="number" and model == "9R"):
        return {'18':'AMYGDALA', '11':'BASALGANGLIA', '16':'BRAINSTEM', 
                            '7':'CEREBELLUM', '251':'CORPUSCALLOSUM', '3':'GREYMATTER',
                            '17':'HIPPOCAMPUS', '10':'THALAMUS', '2':'WHITEMATTER',
                            '24': "CSF" }
    
    if (start =="name" and model == "4R"):
        return {'AMYGDALA':'4', 'BASALGANGLIA':'4', 'BRAINSTEM':'4', 
                    'CEREBELLUM':'4', 'CORPUSCALLOSUM':'251', 'GREYMATTER':'3',
                    'HIPPOCAMPUS':'4', 'THALAMUS':'4', 'WHITEMATTER':'2', 'CSF':'24'
                    }
    if (start =="number" and model == "4R"):
        return {'18':'INTERNAL_STRUCTURES', '11':'INTERNAL_STRUCTURES', '16':'INTERNAL_STRUCTURES', 
                            '7':'INTERNAL_STRUCTURES', '251':'CORPUSCALLOSUM', '3':'GREYMATTER',
                            '17':'INTERNAL_STRUCTURES', '10':'INTERNAL_STRUCTURES', '2':'WHITEMATTER',
                            '24': "CSF" }
    
    if (start =="name" and model == "2R"):
        return {'AMYGDALA':'2', 'BASALGANGLIA':'2', 'BRAINSTEM':'2', 
                    'CEREBELLUM':'2', 'CORPUSCALLOSUM':'2', 'GREYMATTER':'3',
                    'HIPPOCAMPUS':'2', 'THALAMUS':'2', 'WHITEMATTER':'2', 'CSF':'24'
                    }
    if (start =="number" and model == "2R"):
        return {'18':'WHITEMATTER', '11':'WHITEMATTER', '16':'WHITEMATTER', 
                            '7':'WHITEMATTER', '251':'WHITEMATTER', '3':'GREYMATTER',
                            '17':'WHITEMATTER', '10':'WHITEMATTER', '2':'WHITEMATTER',
                            '24': "CSF" }
    
    if (start =="name" and model == "1R"):
        return {'AMYGDALA':'3', 'BASALGANGLIA':'3', 'BRAINSTEM':'3', 
                    'CEREBELLUM':'3', 'CORPUSCALLOSUM':'3', 'GREYMATTER':'3',
                    'HIPPOCAMPUS':'3', 'THALAMUS':'3', 'WHITEMATTER':'3', 'CSF':'24'
                    }
    if (start =="number" and model == "1R"):
        return {'18':'GREYMATTER', '11':'GREYMATTER', '16':'GREYMATTER', 
                            '7':'GREYMATTER', '251':'GREYMATTER', '3':'GREYMATTER',
                            '17':'GREYMATTER', '10':'GREYMATTER', '2':'GREYMATTER',
                            '24': "CSF" }