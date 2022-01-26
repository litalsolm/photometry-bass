# photometry-bass

Step 1: Create csv files that contain information about the objects you wish to study.
ps1: The csv file must contain 3 columns: target (the bass id), ra and dec (example: stars_coor.csv). Let's reffer to this file as "targets file".
sdss: You can insert the targets file to the sdss cross-id tool (http://skyserver.sdss.org/dr16/en/tools/crossid/crossid.aspx). The output cross-id file will help us in the next step. Let's reffer to this file as "cross-id file" (example: Skyserver_SQL11_10_2020 12_19_47_PM.csv).

step 2: Download fits.
Run with: python m arg1 arg2 arg3
Example: python m /home/litalsol/Documents/astro/tables/stars_coor.csv /home/litalsol/Documents/astro/tables/Skyserver_SQL1_5_2021_7_27_01AM.csv /home/litalsol/Documents/astro/fits
Input:
arg1 = the tergets file
arg2 = the cross-id file
arg3 = the path to download the images to

step 3: Extracting data about the ps1 targets.
We insert the targets file to the ps1 catalog search (https://catalogs.mast.stsci.edu/panstarrs/). We can leave the default values checked, they have all the info we need. We hit "search catalog" and download the results as a csv file (example:stars_coor_csv_16_11_2020.csv).

step 4: Running the main script.
Run with: python main.py arg1 arg2 arg3 arg4 arg5
Example: python main.py /home/litalsol/Documents/astro/tables/Skyserver_SQL1_5_2021_7_27_01AM.csv /home/litalsol/Documents/astro/tables/stars_coor_csv_16_11_2020.csv /home/litalsol/Documents/astro/tables/stars_coor.csv /home/litalsol/Documents/astro/fits Aperture.fillfac
Input:
arg1 = the cross-id file
arg2 = the dir to the csv containing the ps1 data on the targets (the one we extracted in step 3)
arg3 = the targets file
arg4 = the path to the dir containing the fits files
arg5 = the type of aperture we want to use (there are two options here: Aperture.fillfac and Aperture.psf)

step 5: view the results.
There is an output csv file called photometry.csv containing the joined data for both surveys.
If you want information on one survey alone, you can always run one of both functions from main: calc_phot_ps1(ps1_file, ps1_targets_file, path) or calc_phot_sdss(sdss_file, path), which provide the following output files: photometry_ps1.csv and photometry_sdss.csv. The path for the output files is hardcoded and it is '~/photometry.csv'.

Notes:
1. There is a table called BAT_catalog_for_crossid.csv that is the targets file for all the targets in bass, you can use in throughtout steps 1-5
2. Some of the sites may produce an output file that contains spaces in the name. For further use, it is recommended to replace the spaces with under scores.
3. The git directory contains 3 python modules: download_fits.py (responsible for downloading the fits files), main.py (contains the script that calculates the photometry) and bass_photometry.py (the back end of the main script).



