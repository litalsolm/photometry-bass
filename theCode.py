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
import photometry
import pandas as pd
import csv

# [g,i,r,u,z]

with open('Skyserver_SQL4_21_2020 8_37_29 AM.csv','r') as csv_file: #opens the file that i get from cross-id
    csv_reader = csv.reader(csv_file)
    coor_lst=np.array([])
    img_lst=np.array([])
    next(csv_reader,None)
    next(csv_reader,None)
    for line in csv_reader:
        coor_lst=np.append(coor_lst,(float(line[2]),float(line[3])))
        img_lst=np.append(img_lst,line[0])

coor_lst = np.reshape(coor_lst,(-1,2))
#phot=photometry.flux_array(coor_lst,dir)
#check=photometry.flux_check(coor_lst,dir)
#c=phot-check
h=photometry.photometry(coor_lst,img_lst)

'''sdss=np.array([[21.78,19.2,20.54,23.29,18.49],[18.04,17.29,17.48,19.48,17.2],
      [20.45,19.36,19.66,22.11,19.29],[18.55,17.82,18.01,19.9,17.77],
      [19.27,17.12,17.85,23.81,16.71],[21.12,18.52,19.71,23.12,17.86],
      [20.42,19.94,18.98,22.26,17.38],[20.4,17.28,18.85,23.96,16.43],[19.03,17.77,18.25,20.07,17.36]])'''

    
'''diff=np.abs(h[0])-sdss
dif = pd.DataFrame(diff,columns=['g','i','r','u','z'])
dif = dif.assign(ra=coor_lst[:,[0]],dec=coor_lst[:,[1]],ID=img_lst)'''

df = pd.DataFrame(h[0],columns=['u','g','r','i','z'])
df = df.assign(ra=coor_lst[:,[0]],dec=coor_lst[:,[1]],ID=img_lst)

df_eplus = pd.DataFrame(h[1],columns=['u','g','r','i','z'])
df_eplus = df_eplus.assign(ra=coor_lst[:,[0]],dec=coor_lst[:,[1]],ID=img_lst)

df_eminus = pd.DataFrame(h[2],columns=['u','g','r','i','z'])
df_eminus = df_eminus.assign(ra=coor_lst[:,[0]],dec=coor_lst[:,[1]],ID=img_lst)

ph=photometry.create_table(h[0],h[1],h[2])
ph = ph.assign(ra=coor_lst[:,[0]],dec=coor_lst[:,[1]],ID=img_lst)
print(ph)

print(df)
print(df_eplus)
print(df_eminus)
#print(dif)


#df_list=[df,df_eplus,df_eminus]
#con=pd.concat(df_list,axis=0)
df.to_csv('stars.csv')
#dif.to_csv('differences.csv')
df_eplus.to_csv('up_error2.csv')
df_eminus.to_csv('down_error2.csv')
ph.to_csv('photometry.csv')







