# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 13:21:06 2019

@author: lital
"""

from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils import SkyCircularAperture
import numpy as np
from photutils import aperture_photometry
from astropy.wcs import WCS
import math
import pandas as pd
import os
from enum import Enum

#to make sure there no typos, the survey can only be given as an enum. in case of more surveys, add more enums and adjust each function.
class Survey(Enum):
    sdss = 1
    ps1 = 2

error_agn = []

bands_ps1 = ['g','r','i','z','y']
bands_sdss = ['u','g','r','i','z']

# creates the data, wcs and header arrays for each star. num is the ID of the AGN that is near the star, 
# thus the name of the file with the fits
# input: num is of type int, survey is of type enrm, as described above.
# output: (data_arr,hdr_arr,wcs_arr) - 3 arrays in size 5. data_arr = 5 sub-arrays containing the counts of the fits files in each band
# in the order of the bands. this array is the image, each sub-array is a matrix in the size of the image (x-pixels x y-pixels).
# hdr_arr = contains the header of each band.
# wcs_arr = contains wcs of each band.
# for either survey, fits files are extracted from well-orgenized files: fits/sdss/name-of-object/BAT_ID_name-of-object_sdss_band.fits 
# for sdss and a similar one for ps1 (described inside the function).
def arrays(num, survey, direc):
    
    assert isinstance(survey, Enum)
    
    n=5
    data_arr = [0]*n
    hdr_arr = [0]*n
    wcs_arr = [0]*n
    for i in range(n):
        try: 
            if survey == Survey.sdss:
                bands = bands_sdss
                path= direc + '/sdss/'+str(num).zfill(4)+'/BAT_ID_'+str(num).zfill(4)+'_sdss_'+bands[i]+'.fits'
            
            elif survey == Survey.ps1:
                bands = bands_ps1
                path = direc + '/ps1/'+str(num).zfill(4)+'/BAT_ID_'+str(num).zfill(4)+'_PS1_'+bands[i]+'_stack.fits'
            
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
        
# in case of a defined radius, such as fixed or psf, finds the radius from the data.
# input: data_dict = the dictionary containing the data obout the radius in the form: {g:{},i:{},...}. j = the int representing the 
# band (for sdss, j=0 represents g); survey = the enum representing the survey.
# output: double representing the radius for the specific object and band. 
def find_rad(data_dict,j,survey): #psf aperture
    if survey == Survey.sdss:
        rad = 1.5 #the psf of sdss is noe available so we use a fixed radius here
    if survey == Survey.ps1: 
        #rad = (data_dict['psfMinorFWHM']+data_dict['psfMajorFWHM'])/2    #when using psf as radius
        rad = data_dict['ApRadius']  # when using fixed radius given by the data
        if data_dict['psfMinorFWHM']==0 and data_dict['psfMajorFWHM']==0:
            rad = 1.5
    return rad

# find the convertion from sky units to pixels (1 sky-unit = 1/conv pixels)
# input: survey = enum representing the survey, hdr = the header file for the survey (the same for all bands).
# output: double representing the convertion.
def get_convertion(survey,hdr):
    convertion = 0
    if survey == Survey.sdss:
        convertion = np.sqrt((hdr[0]['CD1_2']*3600)**2+(hdr[0]['CD2_2']*3600)**2)
    if survey == Survey.ps1:
        convertion = hdr[0]['CDELT1']*3600
    
    return convertion


# find the med value of the background ans it's error in fits units. error is calculated with std of the bg array.
# input: centerX = the center pixel of the object in the image along the X coord; centerY = same in the Y coord; data = the data_arr
# of the image in a certain band (explained in arrays). radii = smal, medium and large radii around the object. the bg is calculated 
# as the median of all the pixel counts between the medium and the large radii. 
# output: medbg = the background to be subtracted from each pixel when calculating photometry. sd = the error in the background.
def photBG(centerX,centerY,data,radii):
    bgColl=np.array([])
    for i in range(math.floor(centerX-radii[2]),math.ceil(centerX+radii[2])):
        for j in range(math.floor(centerY-radii[2]),math.ceil(centerY+radii[2])):
            distance=np.sqrt(((i-centerX)**2)+((j-centerY)**2))
            if distance <= radii[2] and distance >= radii[1]:
                bgColl=np.append(bgColl,data[j][i])
    medbg = np.median(bgColl)
    sd = np.std(bgColl)/np.sqrt(len(bgColl)) #inaccuracy
    return (medbg,sd)

# finds the photometry using aperture_photometry, no bg substracted for one object in each band.
# input: coor = a tuple representing sky coordinates (ea, dec); num = in representing the object ID in bass; survey = as above; 
# data_dict = same as in find_rad.
# output: an array of size 5, each value representing the sum count of the object inside a radius (each survey has different count)
# in a band.
def Skyaperture_agn(coor,num,survey,data_dict, path): #data_dict = the dictionary {g:{}, i:{}...}
    data_arr, hdr_arr, wcs_arr = arrays(num,survey, path)  
    n = len(data_arr) # = 5 bands
    radius = 1.5 #for sdss we havn't determined a radius yet so we'll work with a fixed one
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
                if (survey == Survey.ps1):
                    value = value/data_dict[bands_ps1[j]]['ApFillFac'] #when using const radius
                    #value = value/0.9375 #when usinf psf as radius
                arr[j] = value
                #value = 3.631*(10**(-29))*value #value in erg/sec*cm^2*Hz
            except:
                print("wasn't able to perform sky apeture to AGN "+str(num))
                
    return arr


# calculates the photometry count of each star in all bands in Magnitudes after subtructing the background
# input: same as in Skyaperture_agn
# output: phot_s_mag = an array of size 5 containing the avaluation of the photometry on an object in Magnitudes in all bands; 
# delta_plus = an array representing the up error of the photometry in magnitudes for each band; down_error = same but with doen error.
def phot_agn(coor,num,survey, data_dict, path):
    data_arr, hdr_arr, wcs_arr = arrays(num,survey, path)
    n = len(data_arr) #the number of bands = 5
    conv = get_convertion(survey, hdr_arr) #from arcsec to pixels in each survey
    radius = 1.5
    
    phot = Skyaperture_agn(coor,num,survey,data_dict, path)
    
    #getting a list of the coor in pixels
    pix_lst = [0]*n
    for i in range(n):
        miss = wcs_arr[i] == np.array([-99]) #if a fits file of a band is missing
        if miss.all():
            pix_lst[i] = -99
        else:
            pix_lst[i] = wcs_arr[i].wcs_world2pix(coor[0],coor[1],0)    #the object's coordinates in pixels in the image
    

    med_bg = np.zeros(n) #the background median
    bg_error = np.zeros(n) #the error in the background median
    gain = np.ones(n) # we use it to calculate the error for sdss. gain = data units to photons
    bg_sum = [] # the bg within the inner aperture
    bg_e = [] #the error in the sum of the bg
    
    try:
        for i in range(n): #iterating on the bands
            miss = data_arr[i] == np.array([-99])
            if not miss.all(): # if the fits file exists
                if (survey == Survey.ps1):
                    radius = find_rad(data_dict[bands_ps1[i]],i, survey)/conv
                if survey == Survey.sdss:
                    gain[i]=1/hdr_arr[i]['NMGY'] #gain = data units to photons
                if survey == Survey.ps1:
                    gain[i] = hdr_arr[i]['HIERARCH CELL.GAIN']
                radii = np.array([radius, radius*2, radius*3])
                med_bg[i],bg_error[i] = photBG(pix_lst[i][0],pix_lst[i][1],data_arr[i],radii) #calculating bg and error in bg

            
        bg_sum = med_bg * np.pi*(radii[0]**2) # the bg within the inner aperture
        bg_e = bg_error*np.sqrt(np.pi*(radii[0]**2)) # the error of the bg
    except:
        print("wasn't able to perform photBG to AGN "+str(num))
    
    # we have the photometry in data units with bg, the bg and the error in the bg, now we need to find the photometry after subtructing
    # the bg in magnitudes and the doen and up error of each object photometry.
        
    phot_s = np.zeros(n) # the photometry in data units after subtructing the bg
    phot_error = np.zeros(n) # the error in the photometry before subtructing the bg, calculated as the error of a poisson distribution
    # - the sqrt of the number of photons (we use gain to calculate the photons)
    error = np.zeros(n) # the total error of the sky_aperture photometry and the bg error
    phot_s_mag = np.zeros(n) #the photometry in magnitudes
    m_plus = np.zeros(n) # the photometry in magnitudes after adding the error
    m_minus = np.zeros(n) # the photometry in magnitudes after subtructing the error
    delta_minus = np.zeros(n) # the up error 
    delta_plus = np.zeros(n) # the down error
    
    try:
    
        #this loop calculates the photometry-bg and all the errors. I used a loop because there is an exception when the fits file
        #of one of the bands is missing
        for i in range(n): 
            miss = data_arr[i] == np.array([-99])
            if not miss.all():
                phot_s[i] = phot[i]-bg_sum[i]
                phot_error[i] = 1/gain[i] * np.sqrt(phot[i]*gain[i])
                error[i] = np.sqrt(phot_error[i]**2+bg_e[i]**2)
                
                #the convertion from data units to magnitudes is different for each survey:
                if survey == Survey.sdss:
                    to_mag=lambda x:np.sign(x)*(-2.5)*np.log10(np.sign(x)*x*10**(-9)) 
                
                if survey == Survey.ps1:
                    exp_time = hdr_arr[i]['EXPTIME']
                    to_mag = lambda x : -2.5*np.log10(x)+25 + 2.5*np.log10(exp_time)  #25 is the zero point
                
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


# calculates photometry in all bands to a list of objects using phot_agn
# input: coor_lst = an array of tuples (ra, dec) with sky coordinates; img_lst = an array of objects IDs compatible with coor_lst; 
# survey = same as above; data_array = an array the size of the sample containingthe data_dict of every object.
# output: arr = a matrix containing the photometry of all the objects in all the bands [[object 1], [object 2]...]; arr_eplus = a matrix
# containing the photometry up error of all the objects in all the bands [[object 1], [object 2]...]; arr_eminus = same but with down
# eror.
def photometry(coor_lst,img_lst,survey,data_array, path):
        
    n=len(coor_lst)
    arr=np.array([])
    arr_eplus=np.array([])
    arr_eminus=np.array([])
    for i in range(n): # iterating on the targets
        try:
            if (survey == Survey.ps1):    
                phot_s_mag,delta_plus,delta_minus = phot_agn(coor_lst[i],img_lst[i],survey, data_array[i], path)
            if (survey == Survey.sdss):
                phot_s_mag,delta_plus,delta_minus = phot_agn(coor_lst[i],img_lst[i],survey, {}, path)
            arr = np.append(arr, phot_s_mag)
            arr_eplus = np.append(arr_eplus, delta_plus)
            arr_eminus = np.append(arr_eminus, delta_minus)
        except:
            error_agn.append(img_lst[i])
            print("wasn't able to perform photometry to AGN "+str(img_lst[i]))
    arr=np.reshape(arr,(-1,5))
    arr_eplus=np.reshape(arr_eplus,(-1,5)) # each sub array is in size 5
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




'''#performs photometry on a list of coor using the func 'photometry', which calculates it manually
def flux_array (coor_lst,dir):
  
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
 


