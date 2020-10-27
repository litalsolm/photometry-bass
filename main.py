from astropy.io import fits
from photutils import CircularAperture
from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils import SkyCircularAperture
import numpy as np
from photutils import aperture_photometry
from astropy.wcs import WCS
import bass_photometry
import pandas as pd
import csv
from bass_photometry import Survey
from collections import defaultdict
from astropy.coordinates import SkyCoord


# [g,i,r,u,z]

sdss_file = '/home/litalsol/Documents/astro/tables/Skyserver_SQL5_18_2020 9_05_37 AM.csv' # the sdss cross-id output 
ps1_file = '/home/litalsol/Documents/astro/tables/stars_coor_csv_25_10_2020.csv' # the ps1 catalog search output file
ps1_targets_file = "/home/litalsol/Documents/astro/tables/stars_coor.csv" # the input file for ps1 catalog search

#data array = [{g:{psfMajorFWHM:..., psfMinorFWHM:...},i:{},...},{}] -> to get data of the first target data_array[0]
#coor_ps1 = list of coor tuples [(ra,dec),(ra,dec)...]
#columns_ps1 = the rest of the columns from the file. columns_ps1[_ra_] = all the ra's.
def extract_data_from_ps1(ps1_file,ps1_targets_file):  
    
    lines = list()

    with open(ps1_file, 'r') as csv_file:
        reader = csv.reader(csv_file)
        curr = next(reader)
        lines.append(curr)
        for row in reader:
            lines.append(row)
            if row[2] == curr[2]:
                if row[6] != curr[6]:
                    if float(row[6]) < float(curr[6]):
                        lines.remove(curr)
                    else:
                        lines.remove(row)
                        
                else:
                    if float(row[162]) < float(curr[162]):
                        lines.remove(curr)
                    else:
                        lines.remove(row)
            curr = row
            
    with open(ps1_file, 'w') as writeFile:
        writer = csv.writer(writeFile)
        writer.writerows(lines)
    
    columns_ps1 = defaultdict(list)
    with open(ps1_file,'r') as csv_file: 
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            for (k,v) in row.items():
                columns_ps1[k].append(v)
    
    columns_ps1_targets = defaultdict(list)
    with open(ps1_targets_file,'r') as csv_file: 
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            for (k,v) in row.items():
                columns_ps1_targets[k].append(v)
                
    ra = [float(i) for i in columns_ps1['_ra_']]
    dec = [float(i) for i in columns_ps1['_dec_']] 
    targets = [int(i) for i in columns_ps1['_searchID_']]
    bass_ids = [int(i) for i in columns_ps1_targets['target']]
    
    for i in range(len(targets)):
        targets[i] = str(bass_ids[targets[i]])
    
    zipped = zip(ra,dec)
    coor_ps1 = list(zipped)

    bands = ['g','r','i','z','y']
    columns = ['psfMajorFWHM','psfMinorFWHM','ApFillFac','ApRadius','PSFMag','ApMag'] 
    data_array = [0]*len(ra)
    for i in range(len(ra)):
        data_dict = {}
        for band in bands:
            band_dict={}
            for column in columns:
                band_dict[column] = float(columns_ps1[band+column][i])
            data_dict[band] = band_dict
        data_array[i] = data_dict
        
    
    return (coor_ps1, columns_ps1, data_array, targets)

def extract_data_from_sdss(sdss_file):
    columns_sdss = defaultdict(list)
    with open(sdss_file,'r') as csv_file: #opens the file that i get from cross-id
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            for (k,v) in row.items():
                columns_sdss[k].append(v)
    
    ra = [float(i) for i in columns_sdss['ra']]
    dec = [float(i) for i in columns_sdss['dec']]
    zipped = zip(ra,dec)
    coor_sdss = list(zipped)
    return(coor_sdss, columns_sdss)
    

#coor_sdss, columns_sdss = extract_data_from_sdss(sdss_file)
#data_array_sdss = []
#h_sdss = bass_photometry.photometry(coor_sdss,columns_sdss['col0'],Survey.sdss, data_array_sdss)

'''df = pd.DataFrame(h_sdss[0],columns=['u','g','r','i','z'])
df = df.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['col0'])

df_eplus = pd.DataFrame(h_sdss[1],columns=['u','g','r','i','z'])
df_eplus = df_eplus.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['col0'])

df_eminus = pd.DataFrame(h_sdss[2],columns=['u','g','r','i','z'])
df_eminus = df_eminus.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['col0'])

ph=bass_photometry.create_table(h_sdss[0],h_sdss[1],h_sdss[2],Survey.sdss)
ph = ph.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['col0'])

print(ph)
ph.to_csv('/home/litalsol/Documents/astro/photometry_sdss.csv')'''

def calc_phot_ps1():
    coor_ps1, columns_ps1,data_array_ps1, targets = extract_data_from_ps1(ps1_file,ps1_targets_file)
    h_ps1 = bass_photometry.photometry(coor_ps1,targets,Survey.ps1,data_array_ps1)

    df = pd.DataFrame(h_ps1[0],columns=['g','r','i','z','y'])
    df = df.assign(ra=columns_ps1['_ra_'],dec=columns_ps1['_dec_'],ID=targets)
    
    df_eplus = pd.DataFrame(h_ps1[1],columns=['g','r','i','z','y'])
    df_eplus = df_eplus.assign(ra=columns_ps1['_ra_'],dec=columns_ps1['_dec_'],ID=targets)
    
    df_eminus = pd.DataFrame(h_ps1[2],columns=['g','r','i','z','y'])
    df_eminus = df_eminus.assign(ra=columns_ps1['_ra_'],dec=columns_ps1['_dec_'],ID=targets)
    
    ph=bass_photometry.create_table(h_ps1[0],h_ps1[1],h_ps1[2],Survey.ps1)
    ph = ph.assign(ra=columns_ps1['_ra_'],dec=columns_ps1['_dec_'],ID=targets)
    
    print(ph)
    ph.to_csv('/home/litalsol/Documents/astro/photometry_ps1.csv')
    return (h_ps1, data_array_ps1)


def create_distribution(h,data_array):
    bands = ['g', 'r', 'i', 'z', 'y']
    sample = np.array([])
    for i in range(len(h[0])):
        obj = np.zeros(5)
        for j in range (5):
            obj[j] = data_array[i][bands[j]]['ApMag']
        sample = np.append(sample, obj)
    sample = np.reshape(sample,(-1,5))
        
    diff =  sample - h[0]
    
    obj_diff_avg = np.zeros(len(h[0])) # the avarage of the differnce per object
    band_diff_avg = np.zeros(5) # the avg of the differnce per band
    obj_diff_var = np.zeros(len(h[0]))
    band_diff_var = np.zeros(5)
    
    for i in range(len(h[0])):
        obj_diff_avg[i] = np.average(diff[i])
        obj_diff_var[i] = np.var(diff[i])
    
    for j in range(5):
        band_diff_avg[j] = np.average(diff[:, j])
        band_diff_var[j] = np.var(diff[:,j])
    
    d = {"avarage":band_diff_avg, "variance":band_diff_var}
    df = pd.DataFrame(data=d)
    df.to_csv('/home/litalsol/Documents/astro/diff_const_rad_band.csv')
    
    return (obj_diff_avg,obj_diff_var,band_diff_avg,band_diff_var)



    





