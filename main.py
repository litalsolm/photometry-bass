
import numpy as np
import bass_photometry
import pandas as pd
import csv
from bass_photometry import Survey
from collections import defaultdict
import matplotlib.pyplot as plt 


# [g,i,r,u,z]

sdss_file = '/home/litalsol/Documents/astro/tables/Skyserver_SQL5_18_2020 9_05_37 AM.csv' # the sdss cross-id output 
ps1_file = '/home/litalsol/Documents/astro/tables/stars_coor_csv_25_10_2020.csv' # the ps1 catalog search output file
ps1_targets_ps1_file = "/home/litalsol/Documents/astro/tables/stars_coor.csv" # the input file for ps1 catalog search

#data array = [{g:{psfMajorFWHM:..., psfMinorFWHM:...},i:{},...},{}] -> to get data of the first target data_array[0]
#coor_ps1 = list of coor tuples [(ra,dec),(ra,dec)...]
#columns_ps1 = the rest of the columns from the file. columns_ps1[_ra_] = all the ra's.
def extract_data_from_ps1(ps1_file,ps1_targets_ps1_file):  
    
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
    
    columns_ps1_targets_ps1 = defaultdict(list)
    with open(ps1_targets_ps1_file,'r') as csv_file: 
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            for (k,v) in row.items():
                columns_ps1_targets_ps1[k].append(v)
                
    ra = [float(i) for i in columns_ps1['_ra_']]
    dec = [float(i) for i in columns_ps1['_dec_']] 
    targets_ps1 = [int(i) for i in columns_ps1['_searchID_']]
    bass_ids = [int(i) for i in columns_ps1_targets_ps1['target']]
    
    for i in range(len(targets_ps1)):
        targets_ps1[i] = str(bass_ids[targets_ps1[i]])
    
    zipped = zip(ra,dec)
    coor_ps1 = list(zipped)

    bands = ['g','r','i','z','y']
    columns = ['psfMajorFWHM','psfMinorFWHM','ApFillFac','ApRadius','PSFMag','ApMag','ApMagErr'] 
    data_array = [0]*len(ra)
    for i in range(len(ra)):
        data_dict = {}
        for band in bands:
            band_dict={}
            for column in columns:
                band_dict[column] = float(columns_ps1[band+column][i])
            data_dict[band] = band_dict
        data_array[i] = data_dict
        
    
    return (coor_ps1, columns_ps1, data_array, targets_ps1)

def extract_data_from_sdss(sdss_file):
    columns_sdss = defaultdict(list)
    
    lines = list()
    st = "#Table1"    
    with open(sdss_file, 'r') as readFile:
        reader = csv.reader(readFile)
        for row in reader:
            lines.append(row)
            if row[0] == st:
                lines.remove(row)
    
    with open(sdss_file, 'w') as writeFile:
        writer = csv.writer(writeFile)
        writer.writerows(lines)
            
    
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
    
def calc_phot_sdss():
    coor_sdss, columns_sdss = extract_data_from_sdss(sdss_file)
    data_array_sdss = []
    h_sdss = bass_photometry.photometry(coor_sdss,columns_sdss['col0'],Survey.sdss, data_array_sdss)

    df = pd.DataFrame(h_sdss[0],columns=['u','g','r','i','z'])
    df = df.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['col0'])
    
    df_eplus = pd.DataFrame(h_sdss[1],columns=['u','g','r','i','z'])
    df_eplus = df_eplus.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['col0'])
    
    df_eminus = pd.DataFrame(h_sdss[2],columns=['u','g','r','i','z'])
    df_eminus = df_eminus.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['col0'])
    
    ph=bass_photometry.create_table(h_sdss[0],h_sdss[1],h_sdss[2],Survey.sdss)
    ph = ph.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['col0'])
    
    print(ph)
    ph.to_csv('/home/litalsol/Documents/astro/photometry_sdss.csv')
    
    return h_sdss

def calc_phot_ps1():
    coor_ps1, columns_ps1,data_array_ps1, targets_ps1 = extract_data_from_ps1(ps1_file,ps1_targets_ps1_file)
    h_ps1 = bass_photometry.photometry(coor_ps1,targets_ps1,Survey.ps1,data_array_ps1)

    df = pd.DataFrame(h_ps1[0],columns=['g','r','i','z','y'])
    df = df.assign(ra=columns_ps1['_ra_'],dec=columns_ps1['_dec_'],ID=targets_ps1)
    
    df_eplus = pd.DataFrame(h_ps1[1],columns=['g','r','i','z','y'])
    df_eplus = df_eplus.assign(ra=columns_ps1['_ra_'],dec=columns_ps1['_dec_'],ID=targets_ps1)
    
    df_eminus = pd.DataFrame(h_ps1[2],columns=['g','r','i','z','y'])
    df_eminus = df_eminus.assign(ra=columns_ps1['_ra_'],dec=columns_ps1['_dec_'],ID=targets_ps1)
    
    ph=bass_photometry.create_table(h_ps1[0],h_ps1[1],h_ps1[2],Survey.ps1)
    ph = ph.assign(ra=columns_ps1['_ra_'],dec=columns_ps1['_dec_'],ID=targets_ps1)
    
    print(ph)
    ph.to_csv('/home/litalsol/Documents/astro/photometry_ps1.csv')
    return (h_ps1, data_array_ps1)

def create_common_table():
    coor_sdss, columns_sdss = extract_data_from_sdss(sdss_file)
    data_array_sdss = []
    h_sdss = bass_photometry.photometry(coor_sdss,columns_sdss['col0'],Survey.sdss, data_array_sdss)
    targets_sdss = columns_sdss['col0']
    
    #calculating ps1 photometry and sorting the data
    coor_ps1, columns_ps1,data_array_ps1, targets_ps1 = extract_data_from_ps1(ps1_file,ps1_targets_ps1_file)
    int_targets_ps1 = [int(targets_ps1[i]) for i in range(len(targets_ps1))]
    for i in range(len(data_array_ps1)):
        data_array_ps1[i]['target'] = int_targets_ps1[i]
    data_array_ps1.sort(key = lambda x:x['target'])
    coor_ps1 = [x for _,x in sorted(zip(int_targets_ps1,coor_ps1))]
    columns_ps1 = [x for _,x in sorted(zip(int_targets_ps1,columns_ps1))]
    targets_ps1 = [x for _,x in sorted(zip(int_targets_ps1,targets_ps1))]
    h_ps1 = bass_photometry.photometry(coor_ps1,targets_ps1,Survey.ps1,data_array_ps1)
    
    #combining the 2 target lists
    set_sdss = set(targets_sdss)
    set_ps1 = set(targets_ps1)
    sub_list = list(set_ps1-set_sdss)
    all_targets = targets_sdss+sub_list
    all_targets.sort(key= lambda x:int(x))
    
    bands = ['g','r','i','z']
    
    n = len(all_targets)
    
    #u is only in sdss
    u_sdss = np.zeros(n)
    u_plus_sdss = np.zeros(n)
    u_minus_sdss = np.zeros(n)
    for i in range(len(targets_sdss)):
        index = all_targets.index(targets_sdss[i])
        u_sdss[index] = h_sdss[0][i][0]
        u_plus_sdss[index] = h_sdss[1][i][0]
        u_minus_sdss[index] = h_sdss[2][i][0]
        
    #y is only in ps1
    y_ps1 = np.zeros(n)
    y_plus_ps1 = np.zeros(n)
    y_minus_ps1 = np.zeros(n)
    for i in range(len(targets_ps1)):
        index = all_targets.index(targets_ps1[i])
        y_ps1[index] = h_ps1[0][i][4]
        y_plus_ps1[index] = h_ps1[1][i][4]
        y_minus_ps1[index] = h_ps1[2][i][4]
    
    band_dict = {}
    for i in bands:
        m = 0
        band_ps1 = np.zeros(n)
        band_plus_ps1 = np.zeros(n)
        band_minus_ps1 = np.zeros(n)
        band_sdss = np.zeros(n)
        band_plus_sdss = np.zeros(n)
        band_minus_sdss = np.zeros(n)
        
        for j in range(len(targets_sdss)):
            index = all_targets.index(targets_sdss[j])
            band_sdss[index] = h_sdss[0][j][m+1]
            band_plus_sdss[index] = h_sdss[1][j][m+1]
            band_minus_sdss[index] = h_sdss[2][j][m+1]
        for j in range(len(targets_ps1)):
            index = all_targets.index(targets_ps1[j])
            band_ps1[index] = h_ps1[0][j][m]
            band_plus_ps1[index] = h_ps1[1][j][m]
            band_minus_ps1[index] = h_ps1[2][j][m]
        band_dict[i+'_sdss'] = band_sdss
        band_dict[i+'_plus_sdss'] = band_sdss
        band_dict[i+'_minus_sdss'] = band_sdss
        band_dict[i+'_ps1'] = band_sdss
        band_dict[i+'_plus_ps1'] = band_sdss
        band_dict[i+'_minus_ps1'] = band_sdss
        m+=1
        
    df = pd.DataFrame(np.array([all_targets,u_sdss,u_plus_sdss,u_minus_sdss,band_dict['g_sdss'],band_dict['g_plus_sdss']]))
    
    return(df)

    
    
    


def create_distribution(h,data_array):
    bands = ['g', 'r', 'i', 'z', 'y']
    sample = np.array([])
    error = np.array([])
    for i in range(len(h[0])):
        val_obj = np.zeros(5)
        err_obj = np.zeros(5)
        for j in range (5):
            val_obj[j] = data_array[i][bands[j]]['ApMag']
            err_obj[j] = data_array[i][bands[j]]['ApMagErr']
        sample = np.append(sample, val_obj)
        error = np.append(error, err_obj)
    sample = np.reshape(sample,(-1,5))
    error = np.reshape(error, (-1,5))
    
    result = h[0]
    upper_errors_result = h[1]
    lower_errors_result = h[2]
    diff =  sample - result
    
    obj_diff_med = np.zeros(len(result)) # the avarage of the differnce per object
    band_diff_med = np.zeros(5) # the avg of the differnce per band
    obj_diff_var = np.zeros(len(result))
    band_diff_var = np.zeros(5)
    
    for i in range(len(result)):
        obj_diff_med[i] = np.median(diff[i])
        obj_diff_var[i] = np.std(diff[i])
    
    for j in range(5):
        band_diff_med[j] = np.median(diff[:, j])
        band_diff_var[j] = np.std(diff[:,j])
    
    d = {"median":band_diff_med, "variance":band_diff_var}
    df = pd.DataFrame(data=d)
    df.to_csv('/home/litalsol/Documents/astro/diff_const_rad_band.csv')
       
    for i in range(5):
        y = diff[:,i]
        x = result[:,i]
        upper_error = upper_errors_result[:,i]
        lower_error = lower_errors_result[:,i]
        avg_error = (upper_error + lower_error)/2
        
        '''
        fig, axs = plt.subplots(2)
        fig.suptitle('the difference as a func of the magnitude in band %s' % (bands[i])) 
        axs[0].set(xlabel='ABMag', ylabel='difference')
        axs[0].scatter(x, y, label= "stars", color = "blue")
        
        n = len(x)
        j = 0
        
        #this part is responsible for slicing
        while j<n:
            if x[j] >= 20.5:
                x = np.delete(x,j,0) 
                y = np.delete(y,j,0)
                upper_error = np.delete(upper_error,j,0)
                lower_error = np.delete(lower_error,j,0)
                n -= 1
                j -= 1
            j+=1;
        
        axs[1].scatter(x, y, label= "stars", color = "red")
        axs[1].set(xlabel='ABMag', ylabel='difference')
        #end of part
        
        plt.show() 

        
        y = upper_error
        plt.figure()
        plt.xlabel('ABMag') 
        plt.ylabel('up_error') 
        plt.title('the up errors as a func of the magnitude in band %s' % (bands[i])) 
        plt.scatter(x, y, label= "stars", color= "green",  
            marker= "*", s=30)
        plt.show()'''

        x = avg_error
        y = error[:,i]
        plt.figure()
        plt.xlabel('my calculated error [Mag]') 
        plt.ylabel('PS1 calculated error [Mag]') 
        plt.title('PS1 error vs my error in the band %s' % (bands[i])) 
        plt.scatter(x, y, label= "stars", color= "green",  
            marker= "*", s=30)
        plt.show()
   
    return (obj_diff_med,obj_diff_var,band_diff_med,band_diff_var)


df = create_common_table()
#stars that in sdss only dont appear in the final table
#the colums are rows? and they have no names



    





