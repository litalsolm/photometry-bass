
import numpy as np
import bass_photometry
import pandas as pd
import csv
from bass_photometry import Survey, Aperture
from collections import defaultdict
#import matplotlib.pyplot as plt 
import sys

#sdss_file = '/home/litalsol/Documents/astro/tables/Skyserver_SQL1_5_2021_7_27_01AM.csv' # the sdss cross-id output 
#ps1_file = '/home/litalsol/Documents/astro/tables/stars_coor_csv_16_11_2020.csv' # the ps1 catalog search output file
#ps1_targets_file = "/home/litalsol/Documents/astro/tables/stars_coor.csv" # the input file for ps1 catalog search

# to activate from terminal: >> photometry-bass/photometry-bass >> python main.py /home/litalsol/Documents/astro/tables//home/litalsol/Documents/astro/tables/Skyserver_SQL1_5_2021_7_27_01AM.csv /home/litalsol/Documents/astro/tables/stars_coor_csv_16_11_2020.csv /home/litalsol/Documents/astro/tables/stars_coor.csv /home/litalsol/Documents/astro/fits


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
    
    ps1_targets_columns = defaultdict(list)
    with open(ps1_targets_file,'r') as csv_file: 
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            for (k,v) in row.items():
                ps1_targets_columns[k].append(v)
                
    ra = [float(i) for i in columns_ps1['_ra_']]
    dec = [float(i) for i in columns_ps1['_dec_']] 
    targets_ps1 = [int(i) for i in columns_ps1['_searchID_']]
    bass_ids = [int(i) for i in ps1_targets_columns['target']]
    coor_file = [(ps1_targets_columns['ra'][i], ps1_targets_columns['dec'][i]) for i in range(len(ps1_targets_columns['ra']))]
    
    for i in range(len(targets_ps1)):
        targets_ps1[i] = str(bass_ids[targets_ps1[i]])
    
    zipped = zip(ra,dec)
    coor_ps1 = list(zipped)

    bands = ['g','r','i','z','y']
    columns = ['psfMajorFWHM','psfMinorFWHM','ApFillFac','ApRadius','PSFMag','PSFMagErr','ApMag','ApMagErr'] 
    data_array = [0]*len(ra)
    for i in range(len(ra)):
        data_dict = {}
        for band in bands:
            band_dict={}
            for column in columns:
                band_dict[column] = float(columns_ps1[band+column][i])
            data_dict[band] = band_dict
        data_array[i] = data_dict
        
    #data_array is not in the order of bass id
    return (coor_ps1, columns_ps1, data_array, targets_ps1,coor_file,bass_ids)

def extract_data_from_sdss(sdss_file):
    columns_sdss = defaultdict(list)
    
    # copying the contant of the cross-id file to delete the first row if needed
    lines = list()
    st = "#Table1"    
    with open(sdss_file, 'r') as readFile:
        reader = csv.reader(readFile)
        for row in reader:
            lines.append(row)
            if row[0] == st:
                lines.remove(row)
    
    #writing it back without the #table1 line
    with open(sdss_file, 'w') as writeFile:
        writer = csv.writer(writeFile)
        writer.writerows(lines)
            
    
    # copys the columns of the cross-id file
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
    
# for finding photometry only for sdss
def calc_phot_sdss(sdss_file, path):
    coor_sdss, columns_sdss = extract_data_from_sdss(sdss_file)
    data_array_sdss = []
    h_sdss = bass_photometry.photometry(coor_sdss,columns_sdss['target'],Survey.sdss, data_array_sdss, path)

    df = pd.DataFrame(h_sdss[0],columns=['u','g','r','i','z'])
    df = df.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['target'])
    
    df_eplus = pd.DataFrame(h_sdss[1],columns=['u','g','r','i','z'])
    df_eplus = df_eplus.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['target'])
    
    df_eminus = pd.DataFrame(h_sdss[2],columns=['u','g','r','i','z'])
    df_eminus = df_eminus.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['target'])
    
    ph=bass_photometry.create_table(h_sdss[0],h_sdss[1],h_sdss[2],Survey.sdss, path)
    ph = ph.assign(ra=columns_sdss['ra'],dec=columns_sdss['dec'],ID=columns_sdss['target'])
    
    print(ph)
    ph.to_csv('/home/litalsol/Documents/astro/photometry_sdss.csv')
    
    return h_sdss

# for finding photometry only for ps1
def calc_phot_ps1(ps1_file, ps1_targets_file, path):
    coor_ps1, columns_ps1,data_array_ps1, targets_ps1, coor_file, bass_ids = extract_data_from_ps1(ps1_file,ps1_targets_file)
    h_ps1 = bass_photometry.photometry(coor_ps1,targets_ps1,Survey.ps1,data_array_ps1, path)

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

def create_common_table(sdss_file, ps1_file,ps1_targets_file, path, aperture):
    #calculating sdss photometry
    coor_sdss, columns_sdss = extract_data_from_sdss(sdss_file)
    data_array_sdss = []
    h_sdss = bass_photometry.photometry(coor_sdss,columns_sdss['target'],Survey.sdss, data_array_sdss, path)
    targets_sdss = columns_sdss['target']
    
    #calculating ps1 photometry
    coor_ps1, columns_ps1,data_array_ps1, targets_ps1,coor_file,bass_ids = extract_data_from_ps1(ps1_file,ps1_targets_file)
    h_ps1 = bass_photometry.photometry(coor_ps1,targets_ps1,Survey.ps1,data_array_ps1, path, aperture)
    
    #combining the 2 target lists
    set_sdss = set(targets_sdss)
    set_ps1 = set(targets_ps1)
    sub_list = list(set_ps1-set_sdss)
    all_targets = targets_sdss+sub_list
    all_targets.sort(key= lambda x:int(x))
    
    n = len(all_targets)
    
    #combining the 2 coordibates lists
    all_ra = np.array([])
    all_dec = np.array([])
    for i in range(len(bass_ids)):
        if str(bass_ids[i]) in all_targets:
            all_ra = np.append(all_ra, coor_file[i][0])
            all_dec = np.append(all_dec, coor_file[i][1])
    
    

    info_dict = {}
    
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


    for j in ['g','r','i','z','y']:
        psfMajorFWHM = np.zeros(n)
        psfMinorFWHM = np.zeros(n)
        ApFillFac = np.zeros(n)
        ApRadius = np.zeros(n)
        for i in range(len(targets_ps1)):
            index = all_targets.index(targets_ps1[i])
            psfMajorFWHM[index] = data_array_ps1[i][j]['psfMajorFWHM']
            psfMinorFWHM[index] = data_array_ps1[i][j]['psfMinorFWHM']
            ApFillFac[index] = data_array_ps1[i][j]['ApFillFac']
            ApRadius[index] = data_array_ps1[i][j]['ApRadius']
        info_dict[j+'_psfMajorFWHM'] = psfMajorFWHM
        info_dict[j+'_psfMinorFWHM'] = psfMinorFWHM
        info_dict[j+'_ApFillFac'] = ApFillFac
        info_dict[j+'_ApRadius'] = ApRadius
    
    #for the rest of the bands
    bands = ['g','r','i','z']
    band_dict = {}
    m = 0
    for i in bands:
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
        band_dict[i+'_plus_sdss'] = band_plus_sdss
        band_dict[i+'_minus_sdss'] = band_minus_sdss
        band_dict[i+'_ps1'] = band_ps1
        band_dict[i+'_plus_ps1'] = band_plus_ps1
        band_dict[i+'_minus_ps1'] = band_minus_ps1
        m+=1
        
    d = {'targets':all_targets,'ra':all_ra, 'dec':all_dec,
         'g_ps1':band_dict['g_ps1'],'g_plus_ps1':band_dict['g_plus_ps1'],'g_minus_ps1':band_dict['g_minus_ps1'], 
         'g_psfMajorFWHM':info_dict['g_psfMajorFWHM'],'g_psfMinorFWHM':info_dict['g_psfMinorFWHM'], 'g_ApFillFac':info_dict['g_ApFillFac'], 'g_ApRadius':info_dict['g_ApRadius'],
         'r_ps1':band_dict['r_ps1'],'r_plus_ps1':band_dict['r_plus_ps1'],'r_minus_ps1':band_dict['r_minus_ps1'],
         'r_psfMajorFWHM':info_dict['r_psfMajorFWHM'],'r_psfMinorFWHM':info_dict['r_psfMinorFWHM'], 'r_ApFillFac':info_dict['r_ApFillFac'], 'r_ApRadius':info_dict['r_ApRadius'],
         'i_ps1':band_dict['i_ps1'],'i_plus_ps1':band_dict['i_plus_ps1'],'i_minus_ps1':band_dict['i_minus_ps1'],
         'i_psfMajorFWHM':info_dict['i_psfMajorFWHM'],'i_psfMinorFWHM':info_dict['i_psfMinorFWHM'], 'i_ApFillFac':info_dict['i_ApFillFac'], 'i_ApRadius':info_dict['i_ApRadius'],
         'z_ps1':band_dict['z_ps1'],'z_plus_ps1':band_dict['z_plus_ps1'],'z_minus_ps1':band_dict['z_minus_ps1'],
         'z_psfMajorFWHM':info_dict['z_psfMajorFWHM'],'z_psfMinorFWHM':info_dict['z_psfMinorFWHM'], 'z_ApFillFac':info_dict['z_ApFillFac'], 'z_ApRadius':info_dict['z_ApRadius'],
         'y_ps1':y_ps1, 'y_plus_ps1':y_plus_ps1,'y_minus_ps1':y_minus_ps1,
         'y_psfMajorFWHM':info_dict['y_psfMajorFWHM'],'y_psfMinorFWHM':info_dict['y_psfMinorFWHM'], 'y_ApFillFac':info_dict['y_ApFillFac'], 'y_ApRadius':info_dict['y_ApRadius'],
         'u_sdss':u_sdss, 'u_plus_sdss':u_plus_sdss,'u_minus_sdss':u_minus_sdss,
         'g_sdss':band_dict['g_sdss'],'g_plus_sdss':band_dict['g_plus_sdss'],'g_minus_sdss':band_dict['g_minus_sdss'],
         'r_sdss':band_dict['r_sdss'],'r_plus_sdss':band_dict['r_plus_sdss'],'r_minus_sdss':band_dict['r_minus_sdss'],
         'i_sdss':band_dict['i_sdss'],'i_plus_sdss':band_dict['i_plus_sdss'],'i_minus_sdss':band_dict['i_minus_sdss'],
         'z_sdss':band_dict['z_sdss'],'z_plus_sdss':band_dict['z_plus_sdss'],'z_minus_sdss':band_dict['z_minus_sdss']}

    df = pd.DataFrame(data=d)
    df.to_csv('~/photometry.csv')
    
    return(df, h_ps1, data_array_ps1, h_sdss)

'''

def create_distribution(h, data_array, aperture):
    bands = ['g', 'r', 'i', 'z', 'y']
    sample = np.array([])
    error = np.array([])
    for i in range(len(h[0])):
        val_obj = np.zeros(5)
        err_obj = np.zeros(5)
        if aperture == Aperture.fillfac:
            for j in range (5):
                val_obj[j] = data_array[i][bands[j]]['ApMag']
                err_obj[j] = data_array[i][bands[j]]['ApMagErr']
        else:
            for j in range (5):
                val_obj[j] = data_array[i][bands[j]]['PSFMag']
                err_obj[j] = data_array[i][bands[j]]['PSFMagErr']
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
    
    d = {"median":band_diff_med, "std":band_diff_var}
    df = pd.DataFrame(data=d)
    #df.to_csv('/home/litalsol/Documents/astro/diff_psf_band.csv')
       
    for i in range(5):
        y = diff[:,i]
        x = result[:,i]
        upper_error = upper_errors_result[:,i]
        lower_error = lower_errors_result[:,i]
        avg_error = (upper_error + lower_error)/2
        
        
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
        
        band_diff_med[i] = np.median(y)
        band_diff_var[i] = np.std(y)
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
        plt.show()

        fig, ax = plt.subplots()
        x = avg_error
        y = error[:,i]

        fig.suptitle('PS1 error vs my error in the band %s' % (bands[i])) 
        ax.set(xlabel = 'my calculated error [Mag]', ylabel = 'PS1 calculated error [Mag]')
        ax.scatter(x, y, label= "stars", color= "green",  
            marker= "*", s=30)
        lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
        ]
    
    
        # now plot both limits against eachother
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        plt.show() 
   
    d = {"median":band_diff_med, "std":band_diff_var}
    df = pd.DataFrame(data=d)
    df.to_csv('~/diff_const_rad_band_after_slicing.csv')
    return (obj_diff_med,obj_diff_var,band_diff_med,band_diff_var)

'''

if __name__ == "__main__":
    joint_table, h_ps1, data_array_ps1, h_sdss = create_common_table(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

'''

joint_table, h_ps1, data_array_ps1, h_sdss = create_common_table('/home/litalsol/Documents/astro/tables/Skyserver_SQL1_5_2021_7_27_01AM.csv', 
                    '/home/litalsol/Documents/astro/tables/stars_coor_csv_16_11_2020.csv', 
                    '/home/litalsol/Documents/astro/tables/stars_coor.csv', '/home/litalsol/Documents/astro/fits',
                    Aperture.psf)

a,b,c,d = create_distribution(h_ps1, data_array_ps1, Aperture.psf)


'''