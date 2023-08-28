# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 13:47:49 2023

@author: grife
"""


import nibabel as nb
import numpy as np
import VoxelDataUtils as bm 

def create_aseg_with_lesion(fileInPath, lesion_loc):
    
    filename = "\\mri\\aseg.mgz"
          
    # Step 1: Using freesurfer and 'recon-all' create mri outputs. Ensure aseg.mgz is created.
    t1_file = "\\".join([fileInPath,filename])
    t1 = nb.load(t1_file)
    # t1.orthoview()
    data = np.asarray(t1.dataobj)
    
    # Step 1: Using freesurfer and 'recon-all' create mri outputs. Ensure aseg.mgz is created.
    t1_file_cc = "\\".join([fileInPath,"mri\\cc.mgz"])
    t1_cc = nb.load(t1_file_cc)
    # t1.orthoview()
    cc_data = np.asarray(t1_cc.dataobj)
    
    cc_data = bm.create_binary_image(cc_data)
    bm.override_voxel_data(cc_data, data, 251)
    
    
    
    data_new = np.copy(data);
    for x in range(lesion_loc[0]-4,lesion_loc[0]+4):
        for y in range(lesion_loc[1]-4,lesion_loc[1]+4):
            for z in range(lesion_loc[2]-4,lesion_loc[2]+4):    
                    data_new[x,y,z] = 25;
                    
    # data_new = brainCreator.preprocess(data_new, lesion=True, edemicTissue=1, unusedLabel="Ventricles");
    # bm.clean_lesion(data_new,25)
    # bm.add_edemic_tissue(data_new, 1, 25, 29) 
    # Now we can save the changed data into a new NIfTI file
    new_img = nb.Nifti1Image(data_new, affine=t1.affine, header=t1.header)
    nb.save(new_img, fileInPath + "/tmp/lesion_data.mgz")


# fileInPath = "C:\\Users\grife\OneDrive\Documents\PostDoc\BrainModels\OASIS\OASIS-1" + \
#     "\oasis_cs_freesurfer_disc1.tar\oasis_cs_freesurfer_disc1\disc1\\OAS1_0004_MR1\\"  
# lesion_loc = [100, 118, 97];  
# create_aseg_with_lesion(fileInPath, lesion_loc)