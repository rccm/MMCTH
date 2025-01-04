'''
   Group the coindident files between MISR and MODIS for each MODIS 5-minute granule
   
'''
import re
from datetime import datetime,timedelta
import os,sys
import numpy as np
sys.path.append(os.path.abspath('../'))
from data_readers.mod_read import MODL1Granule
import glob
import fnmatch
import logging

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    handlers=[logging.FileHandler("file_processing.log"),
                              logging.StreamHandler()])

mod06_basedir = '/data/gdi/d/MOD06/'
mod03_basedir = '/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/'
tccloud_basedir = '/data/gdi/e/MISR/TCCloud/'
agp_basedir = '/data/gdi/e/MISR/AGP/'
era_basedir = '/data/gdi/e/ERA5/'

def find_nite_files(mod_l1b_fullpath):
    mod_1b_file = MODL1Granule(mod_l1b_fullpath)
    dayflag = mod_1b_file.get_daynight()
    if dayflag == 'Night':
        logging.warning('Night time granule')
        with open('night_files.txt', 'a') as f:
            f.write(f"{mod_l1b_fullpath}\n")
        return mod_1b_file
    return None

def find_file(basedir, subdir, file_id, prelx):
    pattern = f'{prelx}*{file_id}*'
    matches = []
    for root, _, filenames in os.walk(os.path.join(basedir, subdir)):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches

def find_group_files(mod_l1b_fullpath):
    '''' Find coincident MISR and EAR5 granule for a given MODIS L1B granule'''
    mod_1b_file = MODL1Granule(mod_l1b_fullpath)
    dayflag = mod_1b_file.get_daynight()
    if dayflag == 'Night':
        logging.warning('Night time granule')
        # with open('night_files.txt', 'a') as f:
        #     f.write(f"{mod_l1b_fullpath}\n")
        return None
    
    mod_id = mod_1b_file.get_id()
    orbit = mod_1b_file.get_orbit()
    doy, dd, tt = mod_1b_file.get_datetime()
    mod_subdir = os.path.join(doy[0:4], doy[4:])
#   find mod 03 file
    mod_03 = find_file(mod03_basedir, mod_subdir, mod_id, 'MOD')
    if not mod_03:
        logging.error(f"No matching MOD03 files found for ID {mod_id} in directory {os.path.join(mod03_basedir, mod_subdir)}")
        with open('missing_mod03.txt', 'a') as f:
            f.write(f"{mod_l1b_fullpath}\n")
        return None
    else:
        mod_03 = mod_03[0]
#   find mod 06 file
    mod_06 = find_file(mod06_basedir, mod_subdir, mod_id, 'MOD')
    if not mod_06:
        logging.error(f"No matching MOD06 files found for ID {mod_id} in directory {os.path.join(mod06_basedir, mod_subdir)}")
        # with open('missing_mod06.txt', 'a') as f:
        #     f.write(f"{mod_l1b_fullpath}\n")
        mod_06 = 'X'
    else:
        mod_06 = mod_06[0]

#   find MISR TC_Cloud  file
    misr_subdir = mod_1b_file.get_misrtime()
    tc_cloud = find_file(tccloud_basedir, misr_subdir, str(orbit), 'MISR_AM1_TC_CLOUD_')
    if not tc_cloud:
        date_obj = datetime.strptime(misr_subdir, '%Y.%m.%d')
        # Calculate the day before
        day_before = date_obj - timedelta(days=1)
        # Convert back to the desired format
        day_before_str = day_before.strftime('%Y.%m.%d')
        # Assign to a variable
        previous_misr_subdir = day_before_str
        tc_cloud_second = find_file(tccloud_basedir, previous_misr_subdir, str(orbit), 'MISR_AM1_TC_CLOUD_')
        if not tc_cloud_second:
            logging.error(f"No {misr_subdir} matching TC_Cloud files found for orbit {orbit} and for {mod_l1b_fullpath}")
            # with open('missing_tccloud.txt', 'a') as f:
            #     f.write(f"{orbit}\n")
            if mod_06 == 'X':
                return None    
            tc_cloud = 'X'
            path_num = (176 + int(orbit) * 16) % 233 + 1
            agp_file = find_file(agp_basedir, '', f'P{str(path_num).zfill(3)}', 'MISR_AM1_A')
            if not agp_file:
                agp_file ='X'
            else:
                agp_file = agp_file[0]
        else:
            tc_cloud = tc_cloud_second[0]
            match = re.search(r'P(\d{3})', tc_cloud)
            if match:
                path_num = match.group(1)
                agp_file = find_file(agp_basedir, '', f'P{path_num}', 'MISR_AM1_A')
                if not agp_file:
                    agp_file ='X'
                else:
                    agp_file = agp_file[0]
            else:
                logging.error(f"No path number found in TC_Cloud file name for orbit {orbit}")
    else:
        tc_cloud = tc_cloud[0]
#   find AGP file 
        match = re.search(r'P(\d{3})', tc_cloud)
        if match:
            path_num = match.group(1)
            agp_file = find_file(agp_basedir, '', f'P{path_num}', 'MISR_AM1_A')
            if not agp_file:
                agp_file ='X'
            else:
                agp_file = agp_file[0]
        else:
            logging.error(f"No path number found in TC_Cloud file name for orbit {orbit}")
    
#   find ERA5 single layer file
    
    datetime_obj = datetime.strptime(str(dd), "%Y-%m-%d %H:%M:%S")
    
    era_date_str = datetime_obj.strftime("%Y_%m_%d")
 
    era5_singlefile =   f'{era_basedir}/single/era5_single_{str(era_date_str)}.nc'
  
#   find ERA5 Multi layer file
       
    era5_multifile =   f'{era_basedir}/multi/era5_profile_{str(era_date_str)}.nc'

    datetime_obj_extra = datetime_obj  
    datetime_obj_extra  += timedelta(days=1)
    era_date_str_extra = datetime_obj_extra.strftime("%Y_%m_%d")
    era5_singlefile_extra = f'{era_basedir}single/era5_single_{str(era_date_str_extra)}.nc'
    era5_multifile_extra = f'{era_basedir}multi/era5_profile_{str(era_date_str_extra)}.nc'
  
    return {
        'mod_06': mod_06,
        'mod_03': mod_03,
        'tc_cloud': tc_cloud,
        'agp': agp_file,
        'era5_single': era5_singlefile,
        'era5_multi': era5_multifile,
        'mod_id': mod_id,
        'orbit': orbit,
        'date': dd,
        'path_number': path_num,
        'era5_single_extra': era5_singlefile_extra,
        'era5_multi_extra': era5_multifile_extra,

    }
if __name__ == "__main__":
    # pass
    mod_l1b_fullpath = '/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2002/337/MOD021KM.A2002337.0215.061.2017185135616.hdf'
    a = find_group_files(mod_l1b_fullpath)
    print('.........\n')
    print( mod_l1b_fullpath)
    print(a['mod_06'])
    print(a['mod_03'])
    print(a['tc_cloud'])
    print(a['era5_multi'])
    print(a['agp'])
    print(a['era5_single'])
    print(a['era5_single_extra'])
    print(a['era5_multi_extra'])

    
