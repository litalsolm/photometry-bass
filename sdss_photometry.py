# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 13:21:06 2019

@author: lital
"""

from astropy.io import fits
from photutils import CircularAperture
from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils import SkyCircularAperture
import numpy as np
from photutils import aperture_photometry
from astropy.wcs import WCS
import glob
import math
import pandas as pd
import csv
import os
from enum import Enum, auto

class Survey(Enum):
    sdss = auto()
    ps1 = auto()

error_agn = []

bands_ps1 = ['g','r','i','z','y']
bands_sdss = ['u','g','r','i','z']

# creates the data and header arrays for each star. num is the ID of the AGN that is near the star, 
# thus the name of the file with the fits
def arrays(num, survey):
    
    assert isinstance(survey, Enum)
    
    n=5
    data_arr = [0]*n
    hdr_arr = [0]*n
    wcs_arr = [0]*n
    path = ''
    for i in range(n):
        try: 
            if survey == Survey.sdss:
                bands = bands_sdss
                path='/home/litalsol/Documents/astro/fits/sdss/'+str(num).zfill(4)+'/BAT_ID_'+str(num).zfill(4)+'_sdss_'+bands[i]+'.fits'
            
            elif survey == Survey.ps1:
                bands = bands_ps1
                path = '/home/litalsol/Documents/astro/fits/ps1/'+str(num).zfill(4)+'/BAT_ID_'+str(num).zfill(4)+'_PS1_'+bands[i]+'_stack.fits'
            
            if os.path.exists(path):
                hdul = fits.open(path)
                hdr1=hdul[0].header
                data1=hdul[0].data           
                wcs_arr[i] = WCS(hdr1)
                data_arr[i]=data1
                hdr_arr[i]=hdr1
            
            else:
                wcs_arr[i] = np.array([-99])
                data_arr[i]= np.array([-99])
                hdr_arr[i]= np.array([-99])
        except:
            error_agn.append((num,bands[i]))
            print("wasn't able to read the fits")
    
    return (data_arr,hdr_arr,wcs_arr)
        
def find_rad(data_dict,j,survey): #psf aperture
    if survey == Survey.sdss:
        rad = 1.5
    if survey == Survey.ps1: 
        rad = (data_dict['psfMinorFWHM']+data_dict['psfMajorFWHM'])/2    #when using psf as radius
      #  rad = data_dict['ApRadius']
    return rad

def get_convertion(survey,hdr):
    convertion = 0
    if survey == Survey.sdss:
        convertion = np.sqrt((hdr[0]['CD1_2']*3600)**2+(hdr[0]['CD2_2']*3600)**2)
    if survey == Survey.ps1:
        convertion = hdr[0]['CDELT1']*3600
    
    
    return convertion

#find the med value of the background ans it's error in nanomaggy. error is calculated with std of the bg array.
def photBG(centerX,centerY,data,radii):
    bgColl=np.array([])
    for i in range(math.floor(centerX-radii[2]),math.ceil(centerX+radii[2])):
        for j in range(math.floor(centerY-radii[2]),math.ceil(centerY+radii[2])):
            distance=np.sqrt((i-centerX)**2+(j-centerY)**2)
            if distance<radii[2] and distance>radii[1]:
                bgColl=np.append(bgColl,data[j][i])
    medbg=np.median(bgColl)
    sd=np.std(bgColl)/np.sqrt(len(bgColl)) #inaccuracy
    #sd=np.sqrt(np.abs(medbg))*np.sign(medbg)
    return (medbg,sd)

#find the photometry using aperture_photometry, no bg substracted for one object in each band.
def Skyaperture_agn(coor,num,survey,data_dict): #data_dict = the dictionary {g:{}, i:{}...}
    data_arr, hdr_arr, wcs_arr = arrays(num,survey)  
    n = len(data_arr) # = 5 bands
    radius = 1.5
    arr=np.zeros(n)
    for j in range(n): 
        miss = data_arr[j] == np.array([-99])
        if miss.all():
            arr[j] = -99
        else:
            try:
                if survey == Survey.ps1:
                    radius = find_rad(data_dict[bands_ps1[j]],j, survey)
                position = SkyCoord(coor[0] , coor[1] , unit='deg', frame='icrs')
                aperture = SkyCircularAperture(position, radius*u.arcsec)
                phot_table = aperture_photometry(data_arr[j], aperture, wcs=wcs_arr[j])
                phot_table['aperture_sum'].info.format = '%.8g' 
                value = phot_table['aperture_sum'][0] #value in nanomaggy
                #value=(-2.5)*math.log((value*10**(-9)),10)
                if (survey == Survey.ps1):
                    #value = value/data_dict[bands_ps1[j]]['ApFillFac'] #when using const radius
                    value = value/0.9375 #when usinf psf as radius
                arr[j] = value
                #value = 3.631*(10**(-29))*value #value in erg/sec*cm^2*Hz
            except:
                print("wasn't able to perform sky apeture to AGN "+str(num))
                
    return arr


#subtructs the bg from the result of skyaperture_star
def phot_agn(coor,num,survey, data_dict):
    data_arr, hdr_arr, wcs_arr = arrays(num,survey)
    n = len(data_arr) #the number of bands = 5
    conv = get_convertion(survey, hdr_arr) #from arcsec to pixels in each survey
    radius = 1.5
    
    phot = Skyaperture_agn(coor,num,survey,data_dict)
    
    #getting a list of the coor in pixels
    pix_lst = [0]*n
    for i in range(n):
        miss = wcs_arr[i] == np.array([-99])
        if miss.all():
            pix_lst[i] = -99
        else:
            pix_lst[i] = wcs_arr[i].wcs_world2pix(coor[0],coor[1],0)   
    

    obj=np.zeros(n)
    bg_error=np.zeros(n)
    gain=np.ones(n) # we use it to calculate the error for sdss
    bg_sum=[]
    bg_e=[]
    
    try:
        for i in range(n): #iterating on the bands
            miss = data_arr[i] == np.array([-99])
            if not miss.all():
                if (survey == Survey.ps1):
                    radius = find_rad(data_dict[bands_ps1[i]],i, survey)
                radii = np.array([radius, radius+3, radius+4])/conv
                bg=photBG(pix_lst[i][0],pix_lst[i][1],data_arr[i],radii) #calculating bg
                obj[i] = bg[0]
                bg_error[i] = bg[1]
                if survey == Survey.sdss:
                    gain[i]=hdr_arr[i]['NMGY'] #gain = data units to photons
                if survey == Survey.ps1:
                    gain[i] = 1/hdr_arr[i]['HIERARCH CELL.GAIN']

            
        bg_sum=obj*np.pi*(radii[0]**2) #the bg within the inner aperture
        bg_e=bg_error*np.pi*(radii[0]**2) #the error of the bg
    except:
        print("wasn't able to perform photBG to AGN "+str(num))
    
    
    phot_s = np.zeros(n)
    phot_error = np.zeros(n)
    error = np.zeros(n)
    phot_s_mag = np.zeros(n)
    m_plus = np.zeros(n)
    m_minus = np.zeros(n)
    delta_minus = np.zeros(n)
    delta_plus = np.zeros(n)
    
    try:
    
        #this loop calculates the photometry-bg and all the errors. I used a loop because there is an exception when the fits file
        #of one of the bands is missing
        for i in range(n): 
            miss = data_arr[i] == np.array([-99])
            if not miss.all():
                phot_s[i]=phot[i]-bg_sum[i]
                phot_error[i] = gain[i] * np.sqrt(phot[i]/gain[i])
                error[i] = np.sqrt(phot_error[i]**2+bg_e[i]**2)
                
                if survey == Survey.sdss:
                    to_mag=lambda x:np.sign(x)*(-2.5)*np.log10(np.sign(x)*x*10**(-9))
                
                if survey == Survey.ps1:
                    exp_time = hdr_arr[i]['EXPTIME']
                    to_mag = lambda x : -2.5*np.log10(x)+25 + 2.5*np.log10(exp_time)
                
                phot_s_mag[i]=to_mag(phot_s[i])
                m_plus[i] = to_mag(phot_s[i]-error[i])
                m_minus[i] = to_mag(phot_s[i] + error[i])
                delta_plus[i] = np.absolute(m_plus[i]-phot_s_mag[i])
                delta_minus[i] = np.absolute(m_minus[i]-phot_s_mag[i])
            else: #if there was no info about this band
                phot_s_mag[i] = -99
                delta_minus[i] = -99
                delta_plus[i] = -99
    except:
        print("wasn't able to perform phot_agn to AGN "+str(num))
    
    return (phot_s_mag,delta_plus,delta_minus)


def photometry(coor_lst,img_lst,survey,data_array):
    
    if survey == Survey.sdss:
        bands = bands_sdss
        
    if survey == Survey.ps1:
        bands = bands_ps1
    
    n=len(coor_lst)
    arr=np.array([])
    arr_eplus=np.array([])
    arr_eminus=np.array([])
    for i in range(n): #iterating on the targets
        try:
            if (survey == Survey.ps1):    
                h = phot_agn(coor_lst[i],img_lst[i],survey, data_array[i])
            if (survey == Survey.sdss):
                h = phot_agn(coor_lst[i],img_lst[i],survey, {})
            temp1, temp2, temp3 = h
            arr = np.append(arr, temp1)
            arr_eplus = np.append(arr_eplus, temp2)
            arr_eminus = np.append(arr_eminus, temp3)
        except:
            error_agn.append(img_lst[i])
            print("wasn't able to perform photometry to AGN "+str(img_lst[i]))
    arr=np.reshape(arr,(-1,5))
    arr_eplus=np.reshape(arr_eplus,(-1,5))
    arr_eminus=np.reshape(arr_eminus,(-1,5))
    return (arr,arr_eplus,arr_eminus)


#combines the 3 tables to one dataframe
def create_table(data,plus,minus,survey):
    
    if survey == Survey.sdss:
        df = pd.DataFrame(data[:,[0]],columns=['u'])
        df = df.assign(u_minus=minus[:,[0]])
        df = df.assign(u_plus=plus[:,[0]])
        df=df.assign(g=data[:,[1]])
        df = df.assign(g_minus=minus[:,[1]])
        df = df.assign(g_plus=plus[:,[1]])
        df=df.assign(r=data[:,[2]])
        df = df.assign(r_minus=minus[:,[2]])
        df = df.assign(r_plus=plus[:,[2]])
        df=df.assign(i=data[:,[3]])
        df = df.assign(i_minus=minus[:,[3]])
        df = df.assign(i_plus=plus[:,[3]])
        df=df.assign(z=data[:,[4]])
        df = df.assign(z_minus=minus[:,[4]])
        df = df.assign(z_plus=plus[:,[4]])
    
    if survey == Survey.ps1:
        df = pd.DataFrame(data[:,[0]],columns=['g'])
        df = df.assign(g_minus=minus[:,[0]])
        df = df.assign(g_plus=plus[:,[0]])
        df=df.assign(r=data[:,[1]])
        df = df.assign(r_minus=minus[:,[1]])
        df = df.assign(r_plus=plus[:,[1]])
        df=df.assign(i=data[:,[2]])
        df = df.assign(i_minus=minus[:,[2]])
        df = df.assign(i_plus=plus[:,[2]])
        df=df.assign(z=data[:,[3]])
        df = df.assign(z_minus=minus[:,[3]])
        df = df.assign(z_plus=plus[:,[3]])
        df=df.assign(y=data[:,[4]])
        df = df.assign(y_minus=minus[:,[4]])
        df = df.assign(y_plus=plus[:,[4]])
    
    return df
    

def get_error_agn():
    return error_agn




#performs photometry on a list of coor using the func 'photometry', which calculates it manually
'''def flux_array (coor_lst,dir):
  
    data_arr, hdr_arr, wcs_arr = arrays(dir)
    n=len(coor_lst)
    pix_lst=[[wcs_arr[i].wcs_world2pix(p[0],p[1],0) for p in coor_lst] for i in range(5)]
    pix_lst2=[[pix_lst[i][j] for i in range(5)] for j in range(n)]
    
    #error sqrt flux div sqrt num of pixels
    radii=np.array([1.5,4.5,5.5])
    r_small=radii[0]/0.396
    r_medium=radii[1]/0.396
    r_large=radii[2]/0.369
    #0.396arcsec=1pix
    #pitagoras on  cd21*3600 and cd22*3600
    phot=np.ndarray(shape=(n,len(data_arr))) 
    for i in range(n):
        obj=np.array([])
        for j in range(5):
            obj = np.append(obj,photometry(pix_lst2[i][j][0],pix_lst2[i][j][1],data_arr[j],r_small,r_medium,r_large))
        phot[i] = obj
    return phot
    
#calculates manually photometry given point in pixels, subtructs the background
def photometry(centerX,centerY,data,small,medium,large):
    photColl=np.array([])
    for i in range(math.floor(centerX-small),math.ceil(centerX+small)):
        for j in range(math.floor(centerY-small),math.ceil(centerY+small)):
            distance=np.sqrt((i-centerX)**2+(j-centerY)**2)
            if distance<=small:
                photColl=np.append(photColl,data[j][i])          
            
    medbg=photBG(centerX,centerY,data,small,medium,large)*len(photColl)
    sum=float(np.sum(photColl)-medbg)
    #error=math.sqrt(float(np.sum(photColl))+((math.pi*(small**2))**2)*medbg)
    value=(-2.5)*math.log((sum*10**(-9)),10)
    #mag_error=(-2.5)*math.log((error*10**(-9)),10)
    #value = 3.631*(10**(-29))*value #value in erg/sec*cm^2*Hz
    return value


#returns the array of the bg. not being used in the code
def photBG_arr(centerX,centerY,data,small,medium,large):
    bgColl=np.array([])
    for i in range(math.floor(centerX-large),math.ceil(centerX+large)):
        for j in range(math.floor(centerY-large),math.ceil(centerY+large)):
            distance=np.sqrt((i-centerX)**2+(j-centerY)**2)
            if distance<large and distance>medium:
                bgColl=np.append(bgColl,data[j][i])
    return bgColl

#applies the whole picture as the bg. not being used in the code
def WholeBG(data):
    med=np.median(data)
    #med_error=np.std(data)/np.sqrt(len(data)*len(data[0]))
    med_error=np.std(data)/np.sqrt(1489*2048)
    return (med, med_error)'''
 


