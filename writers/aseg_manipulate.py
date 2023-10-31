# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 13:47:49 2023

@author: grife
"""

import os
import nibabel as nb
import numpy as np
from voxel_data import voxel_data_utils as bm


def add_cc_data(fileInPath,homogenized_data):
    # Step 1: Using freesurfer and 'recon-all' create mri outputs. Ensure aseg.mgz is created.
    t1_file_cc = "/".join([fileInPath, "mri/cc.mgz"])
    t1_cc = nb.load(t1_file_cc)
    # t1.orthoview()
    cc_data = np.asarray(t1_cc.dataobj)

    cc_data = bm.create_binary_image(cc_data)
    bm.override_voxel_data(cc_data, homogenized_data, 251)


def remove_ventricle_data(data_new):
    current_dimensions = data_new.shape
    for x in range(current_dimensions[0]):
        if np.sum(data_new[x, :, :]) > 0:
            for y in range(current_dimensions[1]):
                if np.sum(data_new[x, y, :]) > 0:
                    for z in range(current_dimensions[2]):
                        if data_new[x, y, z] == 4:
                            data_new[x, y, z] = 0


def create_aseg(fileInPath, filename, lesion_loc=[], add_CC=True, remove_ventricle=True):
          
    # Step 1: Using freesurfer and 'recon-all' create mri outputs. Ensure aseg.mgz is created.
    t1_file = "\\".join([fileInPath, filename])
    t1 = nb.load(t1_file)
    # t1.orthoview()
    data = np.asarray(t1.dataobj)
    
    if add_CC:
        add_cc_data(fileInPath,data)

    data_new = np.copy(data)
    if len(lesion_loc > 0):
        for x in range(lesion_loc[0]-4, lesion_loc[0]+4):
            for y in range(lesion_loc[1]-4, lesion_loc[1]+4):
                for z in range(lesion_loc[2]-4, lesion_loc[2]+4):
                    data_new[x, y, z] = 25

    if remove_ventricle:
        remove_ventricle_data(data_new)

    # Now we can save the changed data into a new NIfTI file
    new_img = nb.Nifti1Image(data_new, affine=t1.affine, header=t1.header)
    pathOut = fileInPath + "/tmp"
    file_in_split = filename.split(".")
    file_name_out = file_in_split[0] + "_new"
    if len(file_in_split) == 1:
        file_name_out += '.mgz'
    else:
        file_name_out += file_in_split[1]
    if not os.path.exists(pathOut):
        os.mkdir(pathOut)
    nb.save(new_img, fileInPath + "/tmp/" + file_name_out)
