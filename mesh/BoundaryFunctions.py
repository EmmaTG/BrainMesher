# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 08:51:41 2023

@author: grife
"""

import numpy as np
from abc import ABC, abstractmethod


class IBoundaryTest(ABC):
    """
    Interface for tests on CSF boundaries. If certain criteria have to be met for a boundary element to be added.
    """
    
    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'validElement') and 
                callable(subclass.validElement) or 
                NotImplemented)
    
    @abstractmethod
    def validElement(self, element_num):
        raise NotImplementedError 


class OnlyOnLabel(IBoundaryTest):
    """
    Boundary test to add elements that are only attached to CSF elements 
    """
    
    def __init__(self, mesh, label):
        self.mesh = mesh;
        self.label = label;
    
    def validElement(self, element_num):
        mat = self.mesh.elements[element_num].getMaterial();
        if mat.count(self.label):            
            return True;
        return False;

class BoundaryNormal(IBoundaryTest):

    def __init__(self, mesh, side, value):
        self.mesh = mesh
        self.side = side
        self.value = value
        e_centroids = np.zeros((max(mesh.elements.keys()) + 1, 4))
        count = 0
        for e_num, element in mesh.elements.items():
            element_centroid = element.calculate_element_centroid()
            e_centroids[e_num] = list(element_centroid) + [element.getMaterial()[0]]
            count += 1
        self.e_centroids = np.stack(e_centroids, axis=0)

    def validElement(self, element_num):
        [xc, yc, zc, m] = self.e_centroids[element_num]
        centroid_values = [xc, yc, zc]

            # elements_inline = elements_inline[elements_inline[:, 1].argsort()]
            # elements_inline = elements_inline[elements_inline[:, 2].argsort()]
            # elements_inline = elements_inline[np.where(elements_inline[:, 0] == xc)[0], :]
            # elements_inline = elements_inline[np.where(elements_inline[:, 2] == zc)[0], :]
            # elements_inline = elements_inline[elements_inline[:, 1].argsort()]
            # current_element_idx, = np.where(elements_inline[:, 1] >= yc)
            # if len(current_element_idx) != 0:
            #     current_element_idx = current_element_idx[0]
            #     elements_inline = elements_inline[:current_element_idx]
            # mats = list(elements_inline[:, 3])
            #
            # grey_matter_label = 3
            # white_matter_label = 2
            # csf_label = 24
            #
            # valid_element = (
            #             (mats.count(csf_label) + mats.count(grey_matter_label) + mats.count(white_matter_label)) == len(
            #         mats))
            # return valid_element
        return False
# class BrainStemBase (IBoundaryTest):
#     """
#     Boundary test to add elements that are only attached to CSF elements
#     """
#
#     def __init__(self, mesh, brain_stem_label):
#         self.mesh = mesh
#         self.label = brain_stem_label
#         self.longest_side = 0
#         self.limit
#
#     def ___lowest_value__(self):
#         bounds = [-10000000, -10000000, -10000000, 10000000, 10000000, 10000000]
#         for e in self.mesh.elements.values():
#             if e.getMaterial().count(self.label):
#                 centroid = e.calculate_element_centroid()
#                 for i in range(3):
#                     if centroid[i] > bounds[i]:
#                         bounds[i] = centroid[i]
#                 for i in range(2, 6):
#                     if centroid[i] < bounds[i]:
#                         bounds[i] = centroid[i]
#         difference = np.array(bounds[:3])- np.array(bounds[3:])
#         self.longest_side = np.where(difference == max(difference))[0][0]
#
#
#     def validElement(self, element_num):
#         centroid = volume_elem.calculate_element_centroid()
#         mat = volume_elem.getMaterial()
#         if mat.count(16) and centroid[0] < -57:
#             maxNum += 1
#             boundary_elements[maxNum] = QuadElement(maxNum, face_ica, mat=2)
#         e = self.mesh.elements[element_num]
#         mat = e.getMaterial()
#         centroid = e.calculate_element_centroid()
#         if mat.count(self.label) mat.count(16) and centroid[0] < -57::
#             return True;
#         return False;

class ExternalCSF(IBoundaryTest):
    """
    Boundary test to add elements that are attached to CSF elements and also are not directly below the subcortical structures
    """

    def __init__(self, mesh):
        e_centroids = np.zeros((max(mesh.elements.keys()) + 1, 4))
        count = 0
        for e_num, element in mesh.elements.items():
            element_centroid = element.calculate_element_centroid()
            e_centroids[e_num] = list(element_centroid) + [element.getMaterial()[0]]
            count += 1
        self.e_centroids = np.stack(e_centroids, axis=0)

    def validElement(self, element_num):
        [xc, yc, zc, m] = self.e_centroids[element_num]
        centroid_values = [xc, yc, zc]
        if m == 24:
            elements_inline = self.e_centroids
            for dim in range(3):
                elements_inline = elements_inline[np.where(elements_inline[:, dim] == centroid_values[dim])[0], :]
                elements_inline = elements_inline[elements_inline[:, dim].argsort()]
                current_element_idx, = np.where(elements_inline[:, dim] >= centroid_values[dim])
                if len(current_element_idx) != 0:
                    current_element_idx = current_element_idx[0]
                    elements_inline_A = elements_inline[:current_element_idx]
                    if len(elements_inline_A) == 0:
                        return True
                    elements_inline_B = elements_inline[current_element_idx+1:]
                    if len(elements_inline_B) == 0:
                        return True
            # elements_inline = elements_inline[elements_inline[:, 1].argsort()]
            # elements_inline = elements_inline[elements_inline[:, 2].argsort()]
            # elements_inline = elements_inline[np.where(elements_inline[:, 0] == xc)[0], :]
            # elements_inline = elements_inline[np.where(elements_inline[:, 2] == zc)[0], :]
            # elements_inline = elements_inline[elements_inline[:, 1].argsort()]
            # current_element_idx, = np.where(elements_inline[:, 1] >= yc)
            # if len(current_element_idx) != 0:
            #     current_element_idx = current_element_idx[0]
            #     elements_inline = elements_inline[:current_element_idx]
            # mats = list(elements_inline[:, 3])
            #
            # grey_matter_label = 3
            # white_matter_label = 2
            # csf_label = 24
            #
            # valid_element = (
            #             (mats.count(csf_label) + mats.count(grey_matter_label) + mats.count(white_matter_label)) == len(
            #         mats))
            # return valid_element
        return False

class OpenBottomCSF(IBoundaryTest):
    """
    Boundary test to add elements that are attached to CSF elements and also are not directly below the subcortical structures
    """
    
    def __init__(self, mesh):
        e_centroids = np.zeros((max(mesh.elements.keys())+1,4))
        count = 0
        for e_num,element in mesh.elements.items():
            element_centroid = element.calculate_element_centroid()
            e_centroids[e_num] = list(element_centroid) + [element.getMaterial()[0]]
            count += 1
        self.e_centroids = np.stack(e_centroids, axis = 0)
    
    def validElement(self, element_num):
        [xc,yc,zc,m] = self.e_centroids[element_num]
        if m == 24:            
            elements_inline = self.e_centroids
            elements_inline = elements_inline[np.where(elements_inline[:,0] == xc)[0],:]
            elements_inline = elements_inline[np.where(elements_inline[:,2] == zc)[0],:]
            elements_inline = elements_inline[elements_inline[:, 1].argsort()]
            current_element_idx, = np.where(elements_inline[:,1] >= yc)
            if len(current_element_idx) != 0:
                current_element_idx = current_element_idx[0]
                elements_inline = elements_inline[:current_element_idx]
            mats = list(elements_inline[:, 3])
            
            grey_matter_label = 3
            white_matter_label = 2
            csf_label = 24
            
            valid_element = ((mats.count(csf_label) + mats.count(grey_matter_label) + mats.count(white_matter_label)) == len(mats))
            return valid_element
        return False
