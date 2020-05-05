#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 16:00:37 2020

@author: litalsol
"""

# http://dr16.sdss.org/sas/dr15/eboss/photoObj/frames/301/3918/3/frame-u-003918-3-0213.fits.bz2
#http://dr12.sdss.org/sas/dr12/boss/photoObj/frames/RERUN/RUN/CAMCOL/frame-FILTER-RUN6-CAMCOL-FIELD.fits.bz2

import requests
import csv
import numpy as np
import os
import bz2


'''  very important - if the fits file already exists, there's no need to download it again. I need to add a condition for that,
cause right now the code just downloads everything'''

def download (file_dir):
    with open(file_dir,'r') as csv_file: #opens the file that i get from cross-id
        csv_reader = csv.reader(csv_file)
        ra_lst=np.array([])
        dec_lst=np.array([])
        name_lst=np.array([])
        run_lst=np.array([])
        rerun_lst=np.array([])
        camcol_lst=np.array([])
        field_lst=np.array([])
        
        next(csv_reader,None)
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
        path = '/home/litalsol/Documents/astro/fits/sdss/'+name_lst[i].zfill(4)
        if not os.path.exists(path):
            os.mkdir(path)
        for fil in filter_lst:
            path2 = path + '/'+name_lst[i].zfill(4)+'_'+fil+'.fits'
            if not os.path.exists(path2):
                try:
                    url = 'http://dr16.sdss.org/sas/dr15/eboss/photoObj/frames/'+ rerun_lst[i] +'/' + run_lst[i] + '/' + camcol_lst[i] +'/frame-'+ fil +'-'+ run_lst[i].zfill(6) +'-'+camcol_lst[i]+'-'+field_lst[i]+'.fits.bz2'
                    myfile = requests.get(url)
                    open(path +'/'+name_lst[i].zfill(4)+'_'+fil+'.fits.bz2', 'wb').write(myfile.content)
                    with bz2.open(path +'/'+ name_lst[i].zfill(4)+'_'+fil+'.fits.bz2', "rb") as f:
                        content = f.read()
                        open(path +'/'+name_lst[i].zfill(4)+'_'+fil+'.fits', 'wb').write(content)
                    os.remove(path +'/'+ name_lst[i].zfill(4)+'_'+fil+'.fits.bz2')
                except:
                    print("did not download fits files of ANG "+name_lst[i]+" filter "+fil)
                
        
file_dir = 'Skyserver_SQL4_21_2020 8_37_29 AM.csv'
download(file_dir)


        