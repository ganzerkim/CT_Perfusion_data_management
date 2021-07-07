# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 12:26:51 2021

@author: Mingeon
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

#%%

def phase_split(folder_name):
    base_path = path + '/Phase/' + str(folder_name) 
    images_path = path + '/' + str(folder_name) + '/'

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



#%%

def MIP_maker(folder_name):
    base_path = path + '/phase/' + str(folder_name)

    phase_s = input('start timepoint: ')
    phase_e = input('end timepoint: ')

    slice_path = path + '/slice_img/' + str(folder_name)
    MIP_path = path + '/MIP/' + str(folder_name)

    if not os.path.isdir(slice_path):
        os.makedirs(slice_path)
         
    if not os.path.isdir(MIP_path):
        os.makedirs(MIP_path)

    print("==================== Slice 별 이미지 나누는 중입니다! 잠시만 기다려주세요! ==========================")
    for i in range(int(phase_s), int(phase_e) + 1):
        temp_path = base_path + '/dcm' + str(i)
    
   
        images_list = [s for s in natsort.natsorted(listdir(temp_path)) if isfile(join(temp_path, s))]
    
    
        temp = []
        for idx in range(len(images_list)):
            print(idx)
            db = pydicom.dcmread(temp_path + '/' + images_list[idx], force = True)
            temp.append(db)
            if not os.path.isdir(slice_path + '/' + str(idx)):
                os.makedirs(slice_path + '/' + str(idx))
            db.save_as(slice_path + '/' + str(idx) + '/' + str(i) + '.dcm')
    print("================    Slice 별 이미지 나누기 성공!  ===============================")       

    
    print("=================================   MIP 영상 계산 준비중입니다!   ==================================")
    for cnt in range(len(images_list)):
        slice_list = [s for s in natsort.natsorted(listdir(slice_path + '/' + str(cnt))) if isfile(join(slice_path + '/' + str(cnt) , s))]
  
        dcm_tmp = []
        for cntt in range(len(slice_list)):
            ddd = pydicom.dcmread(slice_path + '/' + str(cnt) + '/' + slice_list[cntt], force = True)
            #shape = np.shape(ddd.pixel_array)
            dcm_tmp.append(ddd.pixel_array)
    
        mipmax = np.array(dcm_tmp)
        mipmax_final = np.amax(mipmax, axis = 0)
        img = sitk.GetImageFromArray(mipmax_final)
        sitk.WriteImage(img, MIP_path + '/' + str(cnt) + ".dcm")

    if not os.path.isdir(MIP_path + '/MIP_Results'):
        os.makedirs(MIP_path + '/MIP_Results')

    mip_list = [s for s in natsort.natsorted(listdir(MIP_path)) if isfile(join(MIP_path, s))]
    for mm in range(len(images_list)):
 
        hdr_tmp = []
        for mmm in range(len(mip_list)):
            ddd = pydicom.dcmread(MIP_path + '/' + mip_list[mmm], force = True)
            hdr_tmp.append(ddd)
        hdr_p = hdr_tmp[mm]
    
        #temp[mm].pixel_array = hdr_p.pixel_array
        temp[mm].PixelData = hdr_p.PixelData
        temp[mm].save_as(MIP_path + '/MIP_Results/' + str(mm) + '+hdr.dcm')
  
    print("=========   MIP 영상 계산완료! MIP_Results 폴더안을 확인하세요!  =============")  

#%%
bpath = input('Base root name ? : ')
path = 'C:/Users/User/Desktop/' + bpath
folder_name = input('Case folder name ? : ')
phase_split(folder_name)
MIP_maker(folder_name)