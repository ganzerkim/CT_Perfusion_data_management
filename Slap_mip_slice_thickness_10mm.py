# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 16:13:25 2021

@author: mingeon
"""


import glob
import pydicom, numpy as np
from os import listdir
from os.path import isfile, join

import matplotlib.pylab as plt
import os

import dicom2nifti
import dicom2nifti.settings as settings

import SimpleITK as sitk

import natsort


def phase_split(folder_name):
    base_path = path + '/Phase/' + str(folder_name) 
    images_path = path + '/source_data/' + str(folder_name) + '/'

    if not os.path.isdir(base_path):
        os.makedirs(base_path)
        
    images_list = [s for s in listdir(images_path) if isfile(join(images_path, s))]
    
    dcm_brain = []
    for i in range(len(images_list)):
        dcm_p = pydicom.dcmread(images_path + images_list[i], force = True)
        dcm_brain.append(dcm_p)

    #plt.imshow(dcm_brain[10].pixel_array, cmap = 'bone')

    Acnum = dcm_brain[5].AcquisitionNumber - 1 #[0x0020, 0x0012] #다른 데이터 셋 통해서 확인 필요
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
        
    return Acnum


def slab_mip_10(folder_name, slab):

    base_path = path + '/phase/' + str(folder_name) + '/'

    MIP_path = path + '/MIP/' + str(folder_name)



    for i in range(Acnum):
        temp_path = base_path + '/dcm' + str(i)
        images_list = [s for s in natsort.natsorted(listdir(temp_path)) if isfile(join(temp_path, s))]

        temp = []
        temp_dcm = []
        for idx in range(len(images_list)):
            db = pydicom.dcmread(temp_path + '/' + images_list[idx], force = True)
            temp_dcm.append(db)
            temp.append(db.pixel_array)
        
        temp = np.array(temp)
    
        slab = slab
    
        mip_temp = []
        for lm in range(len(images_list)- slab + 1):
            print(lm)
            z_mip = []
        
            z_mip.append(temp[lm, :, :])
            z_mip.append(temp[lm + 1, :, :])
            z_mip.append(temp[lm + 2, :, :])
        
            zmax = np.max(z_mip, axis = 0)
        
            mip_temp.append(zmax)

        for ii in range(len(mip_temp)):
            img = sitk.GetImageFromArray(mip_temp[ii])
        
            if not os.path.isdir(MIP_path + '/temp/' + str(i)):
                os.makedirs(MIP_path + '/temp/' + str(i))
        
            sitk.WriteImage(img, MIP_path + '/temp/' + str(i) + '/' + str(ii) + ".dcm")


    mip_list = natsort.natsorted(listdir(MIP_path + '/temp'))

    hdr_tmp = []

    for mm in range(len(mip_list)):
        dcm_list = [s for s in natsort.natsorted(listdir(MIP_path + '/temp/' + mip_list[mm])) if isfile(join(MIP_path + '/temp/' + mip_list[mm], s))]
    
        hdr_tmp = []
        for mmm in range(len(dcm_list)):
            ddd = pydicom.dcmread(MIP_path + '/temp/' + mip_list[mm] + '/' + dcm_list[mmm], force = True)
            hdr_tmp.append(ddd)
        
        
            if not os.path.isdir(MIP_path + '/Results/' + str(mm)):
                os.makedirs(MIP_path + '/Results/' + str(mm))
            
            temp_dcm[mmm].PixelData = ddd.PixelData
            temp_dcm[mmm].save_as(MIP_path + '/Results/' + str(mm) + '/' + str(mmm) + '.dcm')
    
    print("=========   MIP 영상 계산완료! MIP\Results 폴더안을 확인하세요!  =============")  


#%%



bpath = input('Base root name ? : ')
path = 'C:/Users/User/Desktop/' + bpath
folder_name = input('Case folder name ? : ')
Acnum = phase_split(folder_name)
slab_mip_10(folder_name, 3)




