# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 15:02:33 2021

@author: Mingeon
"""


import glob, pylab, pandas as pd
import pydicom, numpy as np
from os import listdir
from os.path import isfile, join
import cv2 as cv

import matplotlib.pylab as plt
import os
import seaborn as sns


#import nibabel
from medpy.io import save, load
import dicom2nifti
#%%
path = 'C:/Users/User/Desktop/sample'
base_path = 'C:/Users/User/Desktop/sample/Phase'

if not os.path.isdir(base_path):
        os.makedirs(base_path)

#%%
#image example

images_path = path + '/p/'
images_list = [s for s in listdir(images_path) if isfile(join(images_path, s))]



dcm_brain = []
for i in range(len(images_list)):
    dcm_p = pydicom.dcmread(images_path + images_list[i], force = True)
    dcm_brain.append(dcm_p)

plt.imshow(dcm_brain[10].pixel_array, cmap = 'bone')

Acnum = dcm_brain[0].AcquisitionNumber - 1 #[0x0020, 0x0012] #다른 데이터 셋 통해서 확인 필요
slices_num = len(dcm_brain) / Acnum

i = 0
for i in range(int(Acnum)):
    savedir = os.path.join(base_path + '/dcm' + str(i))
    if not(os.path.exists(savedir)):
        os.mkdir(savedir)
    
    savedir_nii = os.path.join(base_path + '/nii' + str(i))
    if not(os.path.exists(savedir_nii)):
        os.mkdir(savedir_nii)
    
    for idx in range(i * int(slices_num), int(slices_num) + (i * int(slices_num))):
        dcm_brain[idx].save_as(savedir + '/' + str(i) + '_' + str(idx) + '.dcm')
        print(i, idx)
    
    dicom2nifti.convert_directory(savedir, savedir_nii, compression=False, reorient=True)





import SimpleITK as sitk
import dicom2nifti.settings as settings

def make_mips(image_path, output_dir):
    image = sitk.ReadImage(image_path)
    image_size = image.GetSize()

    basename = os.path.basename(image_path)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    for dim in range(3):
        projection = sitk.MaximumProjection(image, dim)

        if image_size[dim] % 2:  # odd number
            voxel = [0, 0, 0]
            voxel[dim] = (image_size[dim] - 1) / 2
            origin = image.TransformIndexToPhysicalPoint(voxel)
        else:  # even
            voxel1 = np.array([0, 0, 0], int)
            voxel2 = np.array([0, 0, 0], int)
            voxel1[dim] = image_size[dim] / 2 - 1
            voxel2[dim] = image_size[dim] / 2
            point1 = np.array(image.TransformIndexToPhysicalPoint(voxel1.tolist()))
            point2 = np.array(image.TransformIndexToPhysicalPoint(voxel2.tolist()))
            origin = np.mean(np.vstack((point1, point2)), 0)
        projection.SetOrigin(origin)
        projection.SetDirection(image.GetDirection())
        proj_basename = basename.replace('.nii.gz', '_mip_{}.nii.gz'.format(dim))
        sitk.WriteImage(projection, os.path.join(output_dir, proj_basename))
        



    nii_name = [s for s in listdir(savedir_nii) if isfile(join(savedir_nii, s))]
    make_mips(savedir_nii + '/' + nii_name[0], MIP_path + '/' + str(i))
    

imgg = []    
for i in range(10):
    imgg.append(dcm_brain[i].pixel_array)
    
imgg = np.array(imgg)    

def createMIP(np_img, slices_num = 15):
    ''' create the mip image from original image, slice_num is the number of 
    slices for maximum intensity projection'''
    img_shape = np_img.shape
    np_mip = np.zeros(img_shape)
    for i in range(img_shape[0]):
        start = max(0, i-slices_num)
        np_mip[i,:,:] = np.amax(np_img[start:i+1],0)
    return np_mip   

zzz = createMIP(imgg, 10)