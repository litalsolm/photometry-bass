#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 16:00:37 2020

@author: litalsol
"""

import requests
import csv
import numpy as np
import os
import bz2
from bs4 import BeautifulSoup
import sys


# input: the output file from sdss cross-id and path to download the images to.
def download_sdss (file_dir, direc):
    with open(file_dir,'r') as csv_file: #opens the file that i get from cross-id
        csv_reader = csv.reader(csv_file)
        ra_lst=np.array([])
        dec_lst=np.array([])
        name_lst=np.array([])
        run_lst=np.array([])
        rerun_lst=np.array([])
        camcol_lst=np.array([])
        field_lst=np.array([])
        
        
        if next(csv_reader,None)[0] == '#Table1':
            next(csv_reader,None)
        for line in csv_reader:
            name_lst=np.append(name_lst,line[0])
            ra_lst=np.append(ra_lst,float(line[2])) 
            dec_lst=np.append(dec_lst,float(line[3])) 
            run_lst = np.append(run_lst,line[4])  #zfill(6) for the second part of the url
            rerun_lst = np.append(rerun_lst,line[5])
            camcol_lst = np.append(camcol_lst,line[6])
            field_lst = np.append(field_lst,line[7].zfill(4))
            
    n = len(name_lst)
    filter_lst = ['u','g','r','i','z']
    for i in range(n):
        path = direc + '/sdss/' +  name_lst[i].zfill(4)
        if not os.path.exists(path):
            os.mkdir(path)
        for fil in filter_lst:
            path2 = path +'/BAT_ID_'+name_lst[i].zfill(4)+'_sdss_'+fil+'.fits'
            if not os.path.exists(path2):
                try:
                    url = 'http://dr16.sdss.org/sas/dr15/eboss/photoObj/frames/'+ rerun_lst[i] +'/' + run_lst[i] + '/' + camcol_lst[i] +'/frame-'+ fil +'-'+ run_lst[i].zfill(6) +'-'+camcol_lst[i]+'-'+field_lst[i]+'.fits.bz2'
                    myfile = requests.get(url)
                    open(path +'/BAT_ID_'+name_lst[i].zfill(4)+'_sdss_'+fil+'.fits.bz2', 'wb').write(myfile.content)
                    with bz2.open(path +'/BAT_ID_'+ name_lst[i].zfill(4)+'_sdss_'+fil+'.fits.bz2', "rb") as f:
                        content = f.read()
                        open(path +'/BAT_ID_'+name_lst[i].zfill(4)+'_sdss_'+fil+'.fits', 'wb').write(content)
                    os.remove(path +'/BAT_ID_'+name_lst[i].zfill(4)+'_sdss_'+fil+'.fits.bz2')
                except:
                    print("did not download fits files of ANG "+name_lst[i]+" filter "+fil)
                
        
# input: the dir of a csv containing the columns: target, ra, dec, and a path to download thr files to. Downloads the fits files of the targets to
def download_ps1(file_dir, direc):
    url = "https://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=%s+%s&filter=color&filetypes=stack&auxiliary=data&size=480&output_size=0&verbose=0&autoscale=99.500000&catlist="
    ra_lst=np.array([])
    dec_lst=np.array([])
    name_lst=np.array([])
    with open(file_dir,'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader,None)
        for line in csv_reader:
            name_lst=np.append(name_lst,line[0])
            ra_lst=np.append(ra_lst,line[1]) 
            dec_lst=np.append(dec_lst,line[2])
        
    n = len(name_lst)
    for i in range(n):
        try:
            r = requests.get(url % (ra_lst[i], dec_lst[i]))
            soup = BeautifulSoup(r.text, 'html.parser')
            h = soup.find_all('h2')
            if len(h) == 1 or h[1].text.strip() != 'No PS1 3PI images were found at the search position':
                path = direc +'/ps1/' +  name_lst[i].zfill(4)
                filter_lst = find_indexes(soup)
                if not os.path.exists(path):
                    os.mkdir(path)
                for key in filter_lst:
                    if not os.path.exists(path+'/BAT_ID_'+name_lst[i].zfill(4)+'_PS1_'+filter_lst[key]+'_stack.fits'):
                        a = "https:" + soup.find_all('a')[key].get('href')
                        r2 = requests.get(a)
                        open(path+'/BAT_ID_'+name_lst[i].zfill(4)+'_PS1_'+filter_lst[key]+'_stack.fits' , 'wb').write(r2.content)
        except:
            print("could not download item: "+name_lst[i])


def find_indexes(soup): # finds the location of the cutout download fits in the html page
    a = soup.find_all('a')
    suffix = ['g.unconv.fits','r.unconv.fits','i.unconv.fits','z.unconv.fits','y.unconv.fits']
    n = len(a)
    indices = [0]*n
    for i in range(n):
        str = a[i].get('href')
        if str.endswith(suffix[0]):
            indices[0] = i
        if str.endswith(suffix[1]):
            indices[1] = i
        if str.endswith(suffix[2]):
            indices[2] = i
        if str.endswith(suffix[3]):
            indices[3] = i
        if str.endswith(suffix[4]):
            indices[4] = i
    filter_lst = {indices[0]:'g', indices[1]:'r',indices[2]:'i',indices[3]:'z',indices[4]:'y'}
    return filter_lst


def download_all(file_dir_ps1, file_dir_sdss, path):
    download_ps1(file_dir_ps1, path)
    download_sdss(file_dir_sdss, path)
    
#download_ps1('/home/litalsol/Documents/astro/tables/stars_coor.csv', '/home/litalsol/Documents/astro/fits')
download_sdss('/home/litalsol/Documents/astro/tables/Skyserver_SQL1_5_2021_7_27_01AM.csv', '/home/litalsol/Documents/astro/fits') 
#file_dir = '/home/litalsol/Documents/astro/tables/stars_coor.csv' #in the future when I want to download all the fits files, I will replace this file with 'BAT_catalog_for_cross_id.csv'
#file_dir2 = '/home/litalsol/Documents/astro/tables/Skyserver_SQL1_5_2021_7_27_01AM.csv'
#download_sdss(file_dir2)
#download_ps1(file_dir)
'''
if __name__ == "__main__":
 #   download_all(sys.argv[1], sys.argv[2], sys.argv[3])
    #download_ps1(sys.argv[1], sys.argv[3])
    download_sdss(sys.argv[2], sys.argv[3])

        '''