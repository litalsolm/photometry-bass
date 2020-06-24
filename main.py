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
import sdss_photometry
import pandas as pd
import csv
from sdss_photometry import Survey
from collections import defaultdict

# [g,i,r,u,z]

sdss_file = '/home/litalsol/Documents/astro/Skyserver_SQL5_18_2020 9_05_37 AM.csv'
ps1_file = '/home/litalsol/Documents/astro/ps1_crossid.csv'

columns_sdss = defaultdict(list)
with open(sdss_file,'r') as csv_file: #opens the file that i get from cross-id
    csv_reader = csv.DictReader(csv_file)
    for row in csv_reader:
        for (k,v) in row.items():
            columns_sdss[k].append(v)
        
        
columns_ps1 = defaultdict(list)
with open(ps1_file,'r') as csv_file: 
    csv_reader = csv.DictReader(csv_file)
    for row in csv_reader:
        for (k,v) in row.items():
            columns_ps1[k].append(v)

ra = [float(i) for i in columns_sdss['ra']]
dec = [float(i) for i in columns_sdss['dec']]
zipped = zip(ra,dec)
coor_sdss = list(zipped)
h_sdss = sdss_photometry.photometry(coor_sdss,columns_sdss['col0'],Survey.sdss)

ra = [float(i) for i in columns_ps1['_ra_']]
dec = [float(i) for i in columns_ps1['_dec_']]
zipped = zip(ra,dec)
coor_ps1 = list(zipped)
h_ps1 = sdss_photometry.photometry(coor_ps1,columns_ps1['target'],Survey.ps1)

'''sdss=np.array([[21.78,19.2,20.54,23.29,18.49],[18.04,17.29,17.48,19.48,17.2],
      [20.45,19.36,19.66,22.11,19.29],[18.55,17.82,18.01,19.9,17.77],
      [19.27,17.12,17.85,23.81,16.71],[21.12,18.52,19.71,23.12,17.86],
      [20.42,19.94,18.98,22.26,17.38],[20.4,17.28,18.85,23.96,16.43],[19.03,17.77,18.25,20.07,17.36]])'''

    
'''diff=np.abs(h[0])-sdss
dif = pd.DataFrame(diff,columns=['g','i','r','u','z'])
dif = dif.assign(ra=coor_lst[:,[0]],dec=coor_lst[:,[1]],ID=img_lst)'''

df = pd.DataFrame(h_sdss[0],columns=['u','g','r','i','z'])
df = df.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['col0'])

df_eplus = pd.DataFrame(h_sdss[1],columns=['u','g','r','i','z'])
df_eplus = df_eplus.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['col0'])

df_eminus = pd.DataFrame(h_sdss[2],columns=['u','g','r','i','z'])
df_eminus = df_eminus.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['col0'])

ph=sdss_photometry.create_table(h_sdss[0],h_sdss[1],h_sdss[2],Survey.sdss)
ph = ph.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['col0'])
print(ph)

print(df)
print(df_eplus)
print(df_eminus)
#print(dif)
ph.to_csv('/home/litalsol/Documents/astro/photometry_sdss.csv')

df = pd.DataFrame(h_ps1[0],columns=['g','r','i','z','y'])
df = df.assign(ra=columns_ps1['_ra_'],dec=columns_ps1['_dec_'],ID=columns_ps1['target'])

df_eplus = pd.DataFrame(h_ps1[1],columns=['g','r','i','z','y'])
df_eplus = df_eplus.assign(ra=columns_ps1['_ra_'],dec=columns_ps1['_dec_'],ID=columns_ps1['target'])

df_eminus = pd.DataFrame(h_ps1[2],columns=['g','r','i','z','y'])
df_eminus = df_eminus.assign(ra=columns_ps1['_ra_'],dec=columns_ps1['_dec_'],ID=columns_ps1['target'])

ph=sdss_photometry.create_table(h_ps1[0],h_ps1[1],h_ps1[2],Survey.ps1)
ph = ph.assign(ra=columns_ps1['_ra_'],dec=columns_ps1['_dec_'],ID=columns_ps1['target'])
print(ph)

print(df)
print(df_eplus)
print(df_eminus)


#df_list=[df,df_eplus,df_eminus]
#con=pd.concat(df_list,axis=0)
#df.to_csv('stars.csv')
#dif.to_csv('differences.csv')
#df_eplus.to_csv('up_error2.csv')
#df_eminus.to_csv('down_error2.csv')
ph.to_csv('/home/litalsol/Documents/astro/photometry_ps1.csv')







