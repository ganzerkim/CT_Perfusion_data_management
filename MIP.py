# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 15:02:33 2021

@author: github/fepegar
"""


import glob, pylab, pandas as pd
import pydicom, numpy as np
from os import listdir
from os.path import isfile, join
import cv2 as cv
import glob

import matplotlib.pylab as plt
import os

import SimpleITK as sitk

import dicom2nifti
import dicom2nifti.settings as settings
#%%

base_path = 'C:/Users/User/Desktop/sample/phase'

        
phase_s = input('start : ')
phase_e = input('end : ')



slice_path = 'C:/Users/User/Desktop/sample/slice_img'
MIP_path = 'C:/Users/User/Desktop/sample/MIP'

if not os.path.isdir(slice_path):
        os.makedirs(slice_path)
         
if not os.path.isdir(MIP_path):
        os.makedirs(MIP_path)



#%%        
for i in range(int(phase_s), int(phase_e) + 1):
    temp_path = base_path + '/dcm' + str(i)
    
   
    images_list = [s for s in listdir(temp_path) if isfile(join(temp_path, s))]
    
    temp = []
    for idx in range(len(images_list)):
        print(idx)
        db = pydicom.dcmread(temp_path + '/' + images_list[idx], force = True)
        if not os.path.isdir(slice_path + '/' + str(idx)):
            os.makedirs(slice_path + '/' + str(idx))
        db.save_as(slice_path + '/' + str(idx) + '/' + str(i) + '.dcm')
        
    #dicom2nifti.convert_directory(slice_path + '/' + str(idx), MIP_path + '/' + str(idx), reorient=True)  
        
        
#%%
for cnt in range(len(images_list)):
    slice_list = [s for s in listdir(slice_path + '/' + str(cnt)) if isfile(join(slice_path + '/' + str(cnt) , s))]
    dcm_tmp = []
    for cntt in range(len(slice_list)):
        ddd = pydicom.dcmread(slice_path + '/' + str(cnt) + '/' + slice_list[cntt], force = True)
        shape = np.shape(ddd.pixel_array)
        dcm_tmp.append(ddd.pixel_array)
    
    
    mip_final = dcm_tmp[0]    
    for iii in range(len(dcm_tmp)):
        
        mip_final = np.maximum(mip_final, dcm_tmp[iii])
        
        print(iii)
    img = sitk.GetImageFromArray(mip_final)
    sitk.WriteImage(img, MIP_path + '/' + str(cnt) + ".dcm")
        
    
    
#%%  
    
    
settings.disable_validate_orthogonal()
settings.enable_resampling()
settings.set_resample_spline_interpolation_order(1)
settings.set_resample_padding(-1000)
settings.disable_validate_slicecount()
        
def createMIP(np_img, slices_num = 15):
    ''' create the mip image from original image, slice_num is the number of 
    slices for maximum intensity projection'''
    img_shape = np_img.shape
    np_mip = np.zeros(img_shape)
    for i in range(img_shape[0]):
        start = max(0, i-slices_num)
        np_mip[i,:,:] = np.amax(np_img[start:i+1],0)
    return np_mip   
    
    
    
    
    
    
    
for ii in range(len(images_list)):
    slice_list = [s for s in listdir(slice_path + '/' + str(ii)) if isfile(join(slice_path + '/' + str(ii) , s))]
    mip_tmp = []
    dcm_tmp = []
    for iidx in range(len(slice_list)):
        ddd = pydicom.dcmread(slice_path + '/' + str(ii) + '/' + slice_list[iidx], force = True)
        dcm_tmp.append(ddd)
        mip_tmp.append(ddd.pixel_array)
        mip_data = np.array(mip_tmp)
    
    mip_cal = createMIP(mip_data, slices_num = int(len(slice_list))) #int(len(slice_list)
    
    result = np.uint16(mip_cal)
    
    
    
    img = sitk.GetImageFromArray(result[0, :, :])
    
    #hdr = sitk.ReadImage(slice_path + '/' + str(ii) + '/' + slice_list[0])
    #img.CopyInformation(hdr)  
    sitk.WriteImage(img, MIP_path + '/' + str(ii) + ".dcm")
                     
    #dcm_tmp[0].PixelData = img
    #dcm_tmp[0].save_as(MIP_path + '/' + str(ii) + ".dcm")
    print(ii)    


#%%

    
import SimpleITK as sitk
import dicom2nifti.settings as settings

def make_mips(image_path, output_dir):
    im = sitk.ReadImage(image_path)
    image = np.uint16(im)
    image_size = im.GetSize()

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




make_mips(slice_path + '/' + str(idx) + '/' +str(i) + '.dcm', MIP_path + '/' + str(idx))



image_path = slice_path + '/' + str(idx) + '/' +str(i) + '.dcm'





for ii in range(len(images_list)):
    slice_list = [s for s in listdir(slice_path + '/' + str(ii)) if isfile(join(slice_path + '/' + str(ii) , s))]
    mip_tmp = []
    for iidx in range(len(slice_list)):
        ddd = pydicom.dcmread(slice_path + '/' + str(ii) + '/' + slice_list[iidx], force = True)
        mip_tmp.append(ddd.pixel_array)
        mip_tmp = np.array(mip_tmp)
        result = createMIP(mip_tmp, int(len(slice_list)))
        
   
        


