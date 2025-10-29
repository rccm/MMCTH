import os
import glob
import time
import modis
import scipy
import datetime
import numpy as np
import pandas as pd
import xarray as xr
from PIL import Image
import MisrToolkit as Mtk
from scipy import spatial
from pyhdf.SD import SD,SDC
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#------------------------------------------MODIS Parameters---------------------------------------
#NOTE TO SELF: MODIFY YEAR

P_levels = modis.transmission.pstd
#Parameters to be changed ---->
year = 2016;

#-------------------------------------------------------------------------------------------------

#-------------------------Local Directories and Reanalysis Files----------------------------------

main_directory   = '/data/keeling/a/arkam2/c/Reanalysis/era5_files/'+str(year)+'/'
main_directory2  = '/data/keeling/a/arkam2/c/Reanalysis/standard/'
geode            = np.loadtxt('/data/gdi/c/arkam2/Reanalysis/standard/OUTPUT.DAT')
land_mask        = np.loadtxt('/data/keeling/a/arkam2/LandOceanMask/ANCILLARY_1200/data/' + \
                                'land_ocean_masks_xdeg/land_ocean_mask2_hd.asc',skiprows=6)
A                = np.loadtxt('/data/keeling/a/arkam2/b/For_Fusion_'+str(year)+'.txt')

#these are modified versions with the first commented line removed so that read_fwf works okay.
#there is still a dead first column that is removed when gas_mixing_ratios are turned into a xarray dataset.
afglms = pd.read_fwf(os.path.join(main_directory2,'afglms.dat'))
afglmw = pd.read_fwf(os.path.join(main_directory2,'afglmw.dat'))
afglsw = pd.read_fwf(os.path.join(main_directory2,'afglsw.dat'))
afglss = pd.read_fwf(os.path.join(main_directory2,'afglss.dat'))
afglt  = pd.read_fwf(os.path.join(main_directory2,'afglt.dat'))

#-----------------------------Contains Functions--------------------------------------------

#round to nearest multiple of 5 - for CTP solutions
def myround(x, base=5):
    return base * np.round(x/base)

#find indices for the MODIS geolocations:
def find_coord_index(lat,lon):
    nlat = int((90-lat)/0.5)
    nlon = int((lon)/0.5)
    return nlat,nlon

#find indices for the ERA-5 geolocations:
def find_coord_index_reanalysis(d,lat,lon):
    nlat = np.where(np.abs(d.latitude.values-lat)<0.25)[0][0]
    nlon = np.where(np.abs(d.longitude.values-np.abs(lon+180.))<0.25)[0][0]
    return nlat,nlon

#Nearest-neighbour search routine
def do_kdtree(A, pt):
    mytree = scipy.spatial.cKDTree(A,balanced_tree=False,compact_nodes=False)
    dist, indexes = mytree.query(pt,distance_upper_bound=np.inf)
    return dist, indexes

def get_timelimits(h):
    if    h>=0 and h<6:
        h1=0; h2=6
    elif  h>=6 and h<12:
        h1=6; h2=12
    elif  h>=12 and h<18:
        h1=12; h2=18
    elif h>=18 and h<24:
        h1=18; h2=24
    return h1, h2

#-----------------------------Interpolation Routines----------------------------------

def interp3d(data, lat, lon):
    dset = xr.Dataset(
        data_vars={
            'data': (['level','latitude', 'longitude'], data)
        },
        coords={
            'level'    : data.level.data,
            'latitude' : data.latitude.data,
            'longitude': data.longitude.data
        }
    )

    coordinates_to_interpolate_to = xr.Dataset(
        data_vars={
            'level': (['npressures'], P_levels[np.where(P_levels>=1.)]),
            'latitude' : (['point'],      [lat]),
            'longitude': (['point'],      [lon])
        }
    )

    int_data = dset.data.interp(level    =coordinates_to_interpolate_to.level,
                                latitude =coordinates_to_interpolate_to.latitude,
                                longitude=coordinates_to_interpolate_to.longitude
                                )
    
    return int_data

def interp2d(data, lat, lon):
    dset = xr.Dataset(
        data_vars={
            'data': (['latitude', 'longitude'], data)
        },
        coords={
            'latitude' : data.latitude.data,
            'longitude': data.longitude.data
        }
    )

    coordinates_to_interpolate_to = xr.Dataset(
        data_vars={
            'latitude': (['npoints'] , [lat]),
            'longitude': (['npoints'], [lon])
        }
    )

    int_data = dset.data.interp(latitude=coordinates_to_interpolate_to.latitude,
                                longitude=coordinates_to_interpolate_to.longitude
                                )
    
    return int_data

def interp_time(d1, d2, t1, t2, t):
    d = d1 + (d2 - d1)/(t2 - t1) * (t - t1) 
    return d

#-------------------------------------------------------------------------------------------------
#--------------------------Routines to read 6h ERA-5 Reanalysis-------------------------------------
def read_reanalysis_files0(doy):
    date      = xr.open_dataset(os.path.join(main_directory,'era5_0UTC_'+str(year)+'_sst.nc')).time[doy].data
    sst0      = xr.open_dataset(os.path.join(main_directory,'era5_0UTC_'+str(year)+'_sst.nc')).sst[doy,:,:]
    skin0     = xr.open_dataset(os.path.join(main_directory,'era5_0UTC_'+str(year)+'_skin.nc')).skt[doy,:,:]
    sfp0      = xr.open_dataset(os.path.join(main_directory,'era5_0UTC_'+str(year)+'_surfP.nc')).sp[doy,:,:]
    temp0     = xr.open_dataset(os.path.join(main_directory,'era5_0UTC_'+str(year)+'_temp.nc')).t[doy,:,:,:]
    shum0     = xr.open_dataset(os.path.join(main_directory,'era5_0UTC_'+str(year)+'_shum.nc')).q[doy,:,:,:]
    geo0      = xr.open_dataset(os.path.join(main_directory,'era5_0UTC_'+str(year)+'_geo.nc')).z[doy,:,:,:]
    data_2mT0 = xr.open_dataset(os.path.join(main_directory,'era5_0UTC_'+str(year)+'_2mT.nc')).t2m[doy,:,:]
    data_2mD0 = xr.open_dataset(os.path.join(main_directory,'era5_0UTC_'+str(year)+'_2mD.nc')).d2m[doy,:,:]
#     emiss     = Dataset(os.path.join(main_directory2,'surface_emissivity_Fu-Liou_0.5x0.5_53deg.nc'),'r')
    month     = int(str(date).split('-')[1]); month = month - 1
    
    return month, sst0, skin0, sfp0, temp0, shum0, geo0, data_2mT0, data_2mD0

def read_reanalysis_files6(doy):
    date      = xr.open_dataset(os.path.join(main_directory,'era5_6UTC_'+str(year)+'_sst.nc')).time[doy].data
    sst6      = xr.open_dataset(os.path.join(main_directory,'era5_6UTC_'+str(year)+'_sst.nc')).sst[doy,:,:]
    skin6     = xr.open_dataset(os.path.join(main_directory,'era5_6UTC_'+str(year)+'_skin.nc')).skt[doy,:,:]
    sfp6      = xr.open_dataset(os.path.join(main_directory,'era5_6UTC_'+str(year)+'_surfP.nc')).sp[doy,:,:]
    temp6     = xr.open_dataset(os.path.join(main_directory,'era5_6UTC_'+str(year)+'_temp.nc')).t[doy,:,:,:]
    shum6     = xr.open_dataset(os.path.join(main_directory,'era5_6UTC_'+str(year)+'_shum.nc')).q[doy,:,:,:]
    geo6      = xr.open_dataset(os.path.join(main_directory,'era5_6UTC_'+str(year)+'_geo.nc')).z[doy,:,:,:]
    data_2mT6 = xr.open_dataset(os.path.join(main_directory,'era5_6UTC_'+str(year)+'_2mT.nc')).t2m[doy,:,:]
    data_2mD6 = xr.open_dataset(os.path.join(main_directory,'era5_6UTC_'+str(year)+'_2mD.nc')).d2m[doy,:,:]
#     emiss     = Dataset(os.path.join(main_directory2,'surface_emissivity_Fu-Liou_0.5x0.5_53deg.nc'),'r')
    month     = int(str(date).split('-')[1]); month = month - 1
    
    return month, sst6, skin6, sfp6, temp6, shum6, geo6, data_2mT6, data_2mD6

def read_reanalysis_files12(doy):
    date       = xr.open_dataset(os.path.join(main_directory,'era5_12UTC_'+str(year)+'_sst.nc')).time[doy].data
    sst12      = xr.open_dataset(os.path.join(main_directory,'era5_12UTC_'+str(year)+'_sst.nc')).sst[doy,:,:]
    skin12     = xr.open_dataset(os.path.join(main_directory,'era5_12UTC_'+str(year)+'_skin.nc')).skt[doy,:,:]
    sfp12      = xr.open_dataset(os.path.join(main_directory,'era5_12UTC_'+str(year)+'_surfP.nc')).sp[doy,:,:]
    temp12     = xr.open_dataset(os.path.join(main_directory,'era5_12UTC_'+str(year)+'_temp.nc')).t[doy,:,:,:]
    shum12     = xr.open_dataset(os.path.join(main_directory,'era5_12UTC_'+str(year)+'_shum.nc')).q[doy,:,:,:]
    geo12      = xr.open_dataset(os.path.join(main_directory,'era5_12UTC_'+str(year)+'_geo.nc')).z[doy,:,:,:]
    data_2mT12 = xr.open_dataset(os.path.join(main_directory,'era5_18UTC_'+str(year)+'_2mT.nc')).t2m[doy,:,:]
    data_2mD12 = xr.open_dataset(os.path.join(main_directory,'era5_12UTC_'+str(year)+'_2mD.nc')).d2m[doy,:,:]
#     emiss      = Dataset(os.path.join(main_directory2,'surface_emissivity_Fu-Liou_0.5x0.5_53deg.nc'),'r')
    month      = int(str(date).split('-')[1]); month = month - 1
    
    return month, sst12, skin12, sfp12, temp12, shum12, geo12, data_2mT12, data_2mD12

def read_reanalysis_files18(doy):
    date       = xr.open_dataset(os.path.join(main_directory,'era5_18UTC_'+str(year)+'_sst.nc')).time[doy].data
    sst18      = xr.open_dataset(os.path.join(main_directory,'era5_18UTC_'+str(year)+'_sst.nc')).sst[doy,:,:]
    skin18     = xr.open_dataset(os.path.join(main_directory,'era5_18UTC_'+str(year)+'_skin.nc')).skt[doy,:,:]
    sfp18      = xr.open_dataset(os.path.join(main_directory,'era5_18UTC_'+str(year)+'_surfP.nc')).sp[doy,:,:]
    temp18     = xr.open_dataset(os.path.join(main_directory,'era5_18UTC_'+str(year)+'_temp.nc')).t[doy,:,:,:]
    shum18     = xr.open_dataset(os.path.join(main_directory,'era5_18UTC_'+str(year)+'_shum.nc')).q[doy,:,:,:]
    geo18      = xr.open_dataset(os.path.join(main_directory,'era5_18UTC_'+str(year)+'_geo.nc')).z[doy,:,:,:]
    data_2mT18 = xr.open_dataset(os.path.join(main_directory,'era5_18UTC_'+str(year)+'_2mT.nc')).t2m[doy,:,:]
    data_2mD18 = xr.open_dataset(os.path.join(main_directory,'era5_18UTC_'+str(year)+'_2mD.nc')).d2m[doy,:,:]
#     emiss      = Dataset(os.path.join(main_directory2,'surface_emissivity_Fu-Liou_0.5x0.5_53deg.nc'),'r')
    month      = int(str(date).split('-')[1]); month = month - 1
    
    return month, sst18, skin18, sfp18, temp18, shum18, geo18, data_2mT18, data_2mD18

#-------------------------------------------------------------------------------------------------

#----------------------PICK CORRECT 6h ERA-5 FILES FOR INTERPOLATION------------------------------

def find_reanalysis(hr,doy):
    if    hr>=0 and hr<6:
        month1,sst1,skin1,sfp1,temp1,shum1,geo1,data_2mT1,data_2mD1 = read_reanalysis_files0(doy)
        month2,sst2,skin2,sfp2,temp2,shum2,geo2,data_2mT2,data_2mD2 = read_reanalysis_files6(doy)
    elif  hr>=6 and hr<12:
        month1,sst1,skin1,sfp1,temp1,shum1,geo1,data_2mT1,data_2mD1 = read_reanalysis_files6(doy)
        month2,sst2,skin2,sfp2,temp2,shum2,geo2,data_2mT2,data_2mD2 = read_reanalysis_files12(doy)
    elif  hr>=12 and hr<18:
        month1,sst1,skin1,sfp1,temp1,shum1,geo1,data_2mT1,data_2mD1 = read_reanalysis_files12(doy)
        month2,sst2,skin2,sfp2,temp2,shum2,geo2,data_2mT2,data_2mD2 = read_reanalysis_files18(doy)
    elif hr>=18 and hr<24:
        month1,sst1,skin1,sfp1,temp1,shum1,geo1,data_2mT1,data_2mD1 = read_reanalysis_files18(doy)
        month2,sst2,skin2,sfp2,temp2,shum2,geo2,data_2mT2,data_2mD2 = read_reanalysis_files0(doy+1)
    
    return month1,sst1,skin1,sfp1,temp1,shum1,geo1,data_2mT1,data_2mD1,sst2, \
            skin2,sfp2,temp2,shum2,geo2,data_2mT2,data_2mD2

#------------------------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!! ----------- find 
#------------------------------------CREATE ATMOSPHERE AT 101 MOD06 LEVELS-----------------------------------------
def find_reanalysis_vars(day, month, sstemp, skint, sfpres, temp, shum, geo, data_2mT, data_2mD, zouts = 'all',
                        afglms=afglms, afglmw=afglmw, afglsw=afglsw, afglss=afglss, afglt=afglt):
    
    column_lat = sstemp.latitude.data[0]
    column_lon = sstemp.longitude.data[0]

    #random logic to choose the best reference atmosphere for trace gases. . .
    if np.abs(column_lat) < 30.0:
        profile = afglt
        profile_flag = 4
    elif ((column_lat >= 30.0) & (month in (4,5,6,7,8,9))) | ((column_lat <= -30.0) & (month in (1,2,3,10,11,12))):
        if np.abs(column_lat) >= 60.0:
            profile = afglss
        else:
            profile = afglms
            if column_lat < 0:
                profile_flag = 3
            else:
                profile_flag = 1
    else:
        if np.abs(column_lat) >= 60.0:
            profile = afglsw
        else:
            profile = afglmw
            if column_lat < 0:
                profile_flag = 2
            else:
                profile_flag = 0

    T = temp.data[:]
    skt = skint.data
    q = shum.data[:]
    pressures = temp.level.data.astype(np.float32)
    sfp = sfpres.data/100.0
    z_levels = geo.data[:]/9.8e3

    t2 = data_2mT.data
    d2 = data_2mD.data
    spec_hum_2m = 622.0 *6.113e2 * np.exp(5423.0 * (d2 - 273.15)/(d2 * 273.15))/(1000.0 * sfp)

    #remove all levels that are below surface pressure.
    bad_levels = np.where(pressures > sfp)[0]
    if len(bad_levels) > 0:
        min_index = bad_levels.min()
        z_levels = z_levels[:min_index]
        pressures = pressures[:min_index]
        T = T[:min_index]
        q = q[:min_index]
        
    #hack to avoid duplicated z_levels in atmosphere file due to rounding.
    T1 = T[-1]
    if np.round(z_levels[-1],3) == np.array([0.0]):
        z_levels = z_levels[:-1]
        pressures = pressures[:-1]
        T = T[:-1]
        q = q[:-1]
        
    #add surface level data.
    T = T[np.where(z_levels>0.)]
    q = q[np.where(z_levels>0.)]
    pressures = pressures[np.where(z_levels>0.)[0]]
    z_levels = z_levels[np.where(z_levels>0.)[0]]
    z_levels = np.append(z_levels, np.array([0.0]))
    pressures = np.append(pressures, sfp)

    #### what is the difference between sst and 2mt
    if not np.isnan(sstemp.data):
        T = np.append(T, sstemp.data)
    else:
        T = np.append(T, t2)
    if np.isnan(T[len(T)-1]) or np.abs(T[len(T)-1]-T1)>5.:
        T[len(T)-1] = T1+(z_levels[len(T)-2]*6.5)
    q = np.append(q, spec_hum_2m)



    #calculate density
    h2o = q*1e3

    high_col    = profile.iloc[np.logical_and(profile['p(mb)'].values<1., profile['p(mb)'].values>0.0005)]
    low_col     = profile.iloc[profile['p(mb)'].values>=sfp]
    pressures   = np.append(high_col['p(mb)'], pressures) ; z_levels = np.append(high_col['z(km)'], z_levels)
    h2o         = np.append(high_col['h2o(cm-3)'].values/high_col['air(cm-3)'].values*18.016/28.966*1e3, h2o)
    pressures   = np.append(pressures, low_col['p(mb)'])  ; z_levels = np.append(z_levels, low_col['z(km)'])
    h2o         = np.append(h2o, low_col['h2o(cm-3)'].values/low_col['air(cm-3)'].values*18.016/28.966*1e3)
    T           = np.append(high_col['T(K)'], T)          ; T        = np.append(T, low_col['T(K)'])
    T[np.where(pressures>sfp)]   = np.nan
    h2o[np.where(pressures>sfp)] = np.nan 
    z_levels[pressures>sfp]      = np.nan
    fT          = interp1d(np.log(pressures),T,fill_value='extrapolate')
    fwv         = interp1d(np.log(pressures),h2o,fill_value='extrapolate')
    fz          = interp1d(np.log(pressures),z_levels,fill_value='extrapolate')
    T           = fT(np.log(P_levels))
    wv = fwv(np.log(P_levels))
    geo = fz(np.log(P_levels))
    return T, geo, wv, sfp, skt

def amt2tau(amt):
    return 2.13 * -1.*np.log(1. - amt)


def cld_amt(y, m, h2o, T, pressures, vza, ll, delr, rclr):
    taup = modis.transmission.tran_modisd101(y,m,wvmr=h2o,temp=T,ozmr=modis.transmission.ostd,
                                             theta=vza,kban=31,jdet=0,nl=modis.transmission.nl)
    rwcld = modis.emission_rte.fm_modrad_emis(T[ll], 1, pressures[ll], pressures, T, taup, 31, ll)
    amt = delr / (rwcld - rclr)
    
    return amt
#-------------------------------------------------------------------------------------------------

#-------------------------------------GET BIAS CORRECTION ANNUAL NPZ FILE-------------------------
bias_dir        = '/data/keeling/a/arkam2/c/modis_rect/bias_correction/daily_bias_files/combined/'
bias_file       = bias_dir+str(year)+'/MMCTH_CSR_Bias_'+str(year)+'.npz'
bias            = np.load(bias_file)['bias']

f        = open('sample_output_2016_ct.txt','w')
f1       = open('sample_output_2016_op.txt','w')

geode    = np.loadtxt('/data/gdi/c/arkam2/Reanalysis/standard/OUTPUT.DAT')
mod_band = [31,36,35,34,33]

print(len(np.unique(A[:,0])))
for i in range(len(np.unique(A[:,0]))):#len(np.unique(A[:,0]))
    #218 independent scenes
    print(i)
    orbit   = int(A[np.where(A[:,0]==int(np.unique(A[:,0])[i]))[0][0],0])
    start   = int(A[np.where(A[:,0]==int(np.unique(A[:,0])[i]))[0][0],1])
    end     = int(A[np.where(A[:,0]==int(np.unique(A[:,0])[i]))[0][0],2])
    date    = Mtk.orbit_to_time_range(orbit)[0]; path = Mtk.orbit_to_path(orbit); path = str(path).zfill(3)
    hr      = int(date[11:13]);    minute=int(date[14:16]);   year = int(date[0:4]);
    month   = int(date[5:7]);      day   = int(date[8:10])
    mn      = str(month).zfill(2); dy    = str(day).zfill(2); yr   = str(year)
    try:
        misr_file = glob.glob('/data/keeling/a/arkam2/c/MISR/'+yr+'/'+yr+'.'+mn+'.'+dy+'/'+ \
                              'MISR_AM1_TC_CLOUD_P'+path+'_O0'+str(orbit)+'*.hdf')[0]
    except:
        try:
            dy = str(day+1).zfill(2)
            misr_file = glob.glob('/data/keeling/a/arkam2/c/MISR/'+yr+'/'+yr+'.'+mn+'.'+dy+'/'+ \
                                  'MISR_AM1_TC_CLOUD_P'+path+'_O0'+str(orbit)+'*.hdf')[0]
        except:
            pass
    mt             = Mtk.MtkFile(misr_file)
    block_num_time = np.array(mt.block_metadata_field_read('PerBlockMetadataTime', 'BlockCenterTime'))
    block_num_list = np.array(mt.block_metadata_field_read('PerBlockMetadataCommon', 'Block_number'))
    date           = block_num_time[start]
    hr             = int(date[11:13]); minute=int(date[14:16])
    h1, h2         = get_timelimits(hr)
    doy            = datetime.datetime.strptime(date,'%Y-%m-%dT%H:%M:%S.%fZ').timetuple().tm_yday
    month1,sst1,skin1,sfp1,temp1,shum1,geo1,data_2mT1,data_2mD1,sst2, \
        skin2,sfp2,temp2,shum2,geo2,data_2mT2,data_2mD2  =  find_reanalysis(hr,doy)
#     print(year,yr,month,day,doy,hr,minute,mn)
    if (minute<5):
        minute = 55
        hr     = hr-1
    else:
        minute=minute-(minute%5)
    MOD021_file = []; MOD03_file = []; data_33 = []; data_35 = []; data_36 = []; data_31 = []
    mod_lat     = []; mod_lon    = []; data_34 = []; mod_vza = []; data_29 = []; data_32 = []
    day         = str(doy).zfill(3); hour = str(hr).zfill(2); mint = str(minute).zfill(2)
    MOD021_file = np.append(MOD021_file, glob.glob('/data/keeling/a/arkam2/satellite/TerraDataArchive/MODIS/MOD021KM/'+yr+'/'+ \
                                                   day+'/MOD021KM.A'+yr+day+'.'+hour+mint+'*.hdf')[0])
    MOD03_file  = np.append(MOD03_file, glob.glob('/data/keeling/a/arkam2/satellite/TerraDataArchive/MODIS/MOD03/'+yr+'/'+ \
                                                 day+'/MOD03.A'+yr+day+'.'+hour+mint+'*.hdf')[0])
    if (minute == 55):
        if (hr<23):
            minute = 0
            hr     = hr+1
        else:
            minute = 0
            hr     = 0
            doy    = doy+1
    else:
        minute     = minute+5
    day            = str(doy).zfill(3); hour = str(hr).zfill(2); mint = str(minute).zfill(2)
    MOD021_file    = np.append(MOD021_file, 
                               glob.glob('/data/keeling/a/arkam2/satellite/TerraDataArchive/MODIS/MOD021KM/'+yr+'/'+ \
                                                   day+'/MOD021KM.A'+yr+day+'.'+hour+mint+'*.hdf')[0])
    MOD03_file     = np.append(MOD03_file, 
                               glob.glob('/data/keeling/a/arkam2/satellite/TerraDataArchive/MODIS/MOD03/'+yr+'/'+ \
                                                 day+'/MOD03.A'+yr+day+'.'+hour+mint+'*.hdf')[0])
    hdf_geo        = SD(MOD03_file[0], SDC.READ)
    lati           = hdf_geo.select('Latitude'); latit = lati[:,:]
    long           = hdf_geo.select('Longitude'); longit = long[:,:]
    mod_lat        = np.append(mod_lat,latit.ravel())
    mod_lon        = np.append(mod_lon,longit.ravel())
    hdf_geo        = SD(MOD03_file[1], SDC.READ)
    lati           = hdf_geo.select('Latitude'); latit = lati[:,:]
    long           = hdf_geo.select('Longitude'); longit = long[:,:]
    mod_lat        = np.append(mod_lat,latit.ravel())
    mod_lon        = np.append(mod_lon,longit.ravel())
    hdf            = SD(MOD021_file[0], SDC.READ)
    data_id = hdf.select('EV_1KM_Emissive'); attrs = data_id.attributes(full=1)
    offsets = attrs["radiance_offsets"][0]; scales = attrs["radiance_scales"][0]

    #####???????????????
    data_29 = np.append(data_29, np.ravel((data_id[8,:,:]  - offsets[8]) * scales[8])) 
    data_32 = np.append(data_32, np.ravel((data_id[11,:,:] - offsets[11]) * scales[11]))
######################## Guangyu added .... Warning !!!!!!!!! 
    data_31 = np.append(data_33, np.ravel((data_id[10,:,:] - offsets[10]) * scales[10])) 
    data_33 = np.append(data_33, np.ravel((data_id[12,:,:] - offsets[12]) * scales[12]))
    data_34 = np.append(data_34, np.ravel((data_id[13,:,:] - offsets[13]) * scales[13]))
    data_35 = np.append(data_35, np.ravel((data_id[14,:,:] - offsets[14]) * scales[14]))
    data_36 = np.append(data_36, np.ravel((data_id[15,:,:] - offsets[15]) * scales[15]))
    hdf     = SD(MOD021_file[1], SDC.READ)
    data_id = hdf.select('EV_1KM_Emissive'); attrs = data_id.attributes(full=1)
    offsets = attrs["radiance_offsets"][0]; scales = attrs["radiance_scales"][0]
    data_29 = np.append(data_29, np.ravel((data_id[8,:,:]  - offsets[8]) * scales[8])) 
    data_32 = np.append(data_32, np.ravel((data_id[11,:,:] - offsets[11]) * scales[11]))
    data_31 = np.append(data_33, np.ravel((data_id[10,:,:] - offsets[10]) * scales[10])) 
    data_33 = np.append(data_33, np.ravel((data_id[12,:,:] - offsets[12]) * scales[12]))
    data_34 = np.append(data_34, np.ravel((data_id[13,:,:] - offsets[13]) * scales[13]))
    data_35 = np.append(data_35, np.ravel((data_id[14,:,:] - offsets[14]) * scales[14]))
    data_36 = np.append(data_36, np.ravel((data_id[15,:,:] - offsets[15]) * scales[15]))
    modis_points = np.dstack([mod_lat,mod_lon])[0]
    for j in range(len(A[np.where(A[:,0]==orbit)])): #len(A[np.where(A[:,0]==orbit)])
        index = np.where(A[:,0]==orbit)[0][j]; cp=np.zeros(4); ctp=np.zeros(4)
        if (A[index,6]>-0.5 and A[index,5]>0. and A[index,9]<5 and (A[index,7]-A[index,6])>1.):
            ctp           = []; cp = []
            cats_lat      = A[index,3]; cats_lon = A[index,4]; vza=A[index,13]
            CATS_coord    = [cats_lat, cats_lon]
            dist, results = do_kdtree(modis_points, CATS_coord)
            rad           = np.zeros(len(mod_band));        bi        = np.zeros(len(mod_band))
            rad[0]        = data_36[results]       ;        rad[1]    = data_35[results]; 
            rad[2]        = data_34[results]       ;        rad[3]    = data_33[results]; rad[4] = data_31[results]
            nlat,nlon     = find_coord_index_reanalysis(sst1,cats_lat,cats_lon)
            n, m          = find_coord_index(sst1.latitude.data[nlat],sst1.longitude.data[nlon])
            if land_mask[n, m] == 1:
                ls        = 0
            else:
                ls        = 3
            bi            = bias[int(doy)-1, find_coord_index_reanalysis(sst1,cats_lat,cats_lon)[0],:,ls]
            bias[np.where(np.isnan(bias)==True)] = 0
            month   = pd.to_datetime(doy-1, unit='D', origin=str(year)).month
            tmp1    = interp3d(temp1,cats_lat,cats_lon+180); tmp2 = interp3d(temp2,cats_lat,cats_lon+180)
            shm1    = interp3d(shum1,cats_lat,cats_lon+180); shm2 = interp3d(shum2,cats_lat,cats_lon+180)
            gpt1    = interp3d(geo1,cats_lat,cats_lon+180) ; gpt2 = interp3d(geo2,cats_lat,cats_lon+180)
            stp1    = interp2d(sst1,cats_lat,cats_lon+180) ; stp2 = interp2d(sst2,cats_lat,cats_lon+180)
            skt1    = interp2d(skin1,cats_lat,cats_lon+180); skt2 = interp2d(skin2,cats_lat,cats_lon+180)
            spr1    = interp2d(sfp1,cats_lat,cats_lon+180) ; spr2 = interp2d(sfp2,cats_lat,cats_lon+180)
            t2m1    = interp2d(data_2mT1,cats_lat,cats_lon+180)
            t2m2    = interp2d(data_2mT2,cats_lat,cats_lon+180)
            d2m1    = interp2d(data_2mD1,cats_lat,cats_lon+180)
            d2m2    = interp2d(data_2mD2,cats_lat,cats_lon+180)
            T1, g1, h2o1, sp1, st1 = find_reanalysis_vars(doy, month, \
                     stp1, skt1, spr1, tmp1, shm1, gpt1, t2m1, d2m1)
            T2, g2, h2o2, sp2, st2 = find_reanalysis_vars(doy, month, \
                     stp2, skt2, spr2, tmp2, shm2, gpt2, t2m2, d2m2)            
            T             = interp_time(T1, T2, h1, h2, hr+float(mint)/60.)
            h2o           = interp_time(h2o1, h2o2, h1, h2, hr+float(mint)/60.)
            gz            = interp_time(g1, g2, h1, h2, hr+float(mint)/60.)
            sfp           = interp_time(sp1, sp2, h1, h2, hr+float(mint)/60.)
            skt           = interp_time(st1, st2, h1, h2, hr+float(mint)/60.)
            delz          = float(geode[np.where(np.logical_and(int(geo1.latitude.data[nlat])==geode[:,0],
                             int(geo1.longitude.data[nlon])==geode[:,1])),2])/1e3
            fz            = interp1d(geo1[:,nlat,nlon]/9.8e3, np.log(geo1.level.data), fill_value='extrapolate') 
            fP            = interp1d(np.log(geo1.level.data), geo1[:,nlat,nlon]/9.8e3, fill_value='extrapolate') 
            catsp         = np.exp(fz(A[index,5]+delz)); misrtcold = [modis.planck_bright.modis_bright(rad=rad[kban],band=mod_band[kban],units=1) p = np.exp(fz(A[index,6]+delz))
          
                                                  for kban in range(len(mod_band))]
            met_date      = np.array([int(year),int(month)+1,int(day),int(hr)])
            cp, ec1       = modis.co2cld_onepixel(h2o,T,sfp,sfp,skt,vza,tcold,met_date,
                                                  mod_lat[results],mod_lon[results],
                                                  land_mask[n,m],bi)
            
            ll0           = np.where(np.abs(modis.transmission.pstd - cp[0]) == 
                          np.abs(modis.transmission.pstd - cp[0]).min())[0][0]
            ll1           = np.where(np.abs(modis.transmission.pstd - cp[2]) == 
                          np.abs(modis.transmission.pstd - cp[2]).min())[0][0]
            amt11         = cld_amt(year, month, h2o, T, modis.transmission.pstd, vza, ll0, ec1[1], ec1[0])
            amt12         = cld_amt(year, month, h2o, T, modis.transmission.pstd, vza, ll1, ec1[1], ec1[0])
            
            imisr         = np.where(np.abs(P_levels-misrp) == np.min(np.abs(P_levels-misrp)))[0][0]
            ########## Pixel Level Processing ###########
            ctp, ec2      = modis.co2cld_onepixel_misr(h2o,T,sfp,sfp,skt,vza,tcold[:5],met_date,mod_lat[results],
                                             mod_lon[results],land_mask[n,m],bi,misrp)
            
            f.write(str(cp[0])+"  "+str(cp[2])+"  "+str(A[index,8])+"  "+str(A[index,9])+"  "+str(catsp)+" "+
                        str(misrp)+"  "+str(ctp[0])+"  "+str(ctp[2])+"  "+str(np.exp(fz(A[index,5]-A[index,12]+delz)))+"  "+
                        str(fP(np.log(cp[0]))+delz)+"    "+str(fP(np.log(cp[2]))+delz)
                        +"    "+str(fP(np.log(ctp[0]))+delz)+"    "+str(fP(np.log(ctp[2]))+delz)+"    "+
                        str(A[index,5])+"    "+str(A[index,7])+"    "+str(A[index,6])+"    "+"\n")

            print(index, str(cp[0])+"  "+str(cp[2])+"  "+str(A[index,8])+"  "+str(A[index,9])+"  "+str(catsp)+" "+
                        str(misrp)+"  "+str(ctp[0])+"  "+str(ctp[2]))

            ll0 = np.where(np.abs(modis.transmission.pstd - ctp[0]) == 
                          np.abs(modis.transmission.pstd - ctp[0]).min())[0][0]
            ll1 = np.where(np.abs(modis.transmission.pstd - ctp[2]) == 
                          np.abs(modis.transmission.pstd - ctp[2]).min())[0][0]
            amt21 = cld_amt(year, month, h2o, T, modis.transmission.pstd, vza, ll0, ec2[1], ec2[0])
            amt22 = cld_amt(year, month, h2o, T, modis.transmission.pstd, vza, ll1, ec2[1], ec2[0])
                
            f1.write(str(amt11)+"  "+str(amt12)+"  "+str(amt21)+"  "+str(amt22)+"  "+str(A[index,13]*0.00999)+"  "+
                    str(1-np.exp(-1*(3.70939087e+01*A[index,10] + 1.07857124e-02)/2.13))+"   "+
                    str(A[index,11])+"   "+str(3.70939087e+01*A[index,10] + 1.07857124e-02)+
                    "\n")
#             print(index, str(amt11)+"  "+str(amt12)+"  "+str(amt21)+"  "+str(amt22)+"  "+str(A[index,13]*0.00999)+"  "+
#                     str(1-np.exp(-1*(3.70939087e+01*A[index,10] + 1.07857124e-02)/2.13))+"   "+
#                     str(amt2tau(amt11))+"   "+str(amt2tau(amt12))+"   "+str(amt2tau(amt21))+"   "+
#                     str(amt2tau(amt22))+"   "+str(A[index,11])+"   "+str(3.70939087e+01*A[index,10] + 1.07857124e-02)+
#                     "\n")            
#     except:
#         continue
f.close(); f1.close()