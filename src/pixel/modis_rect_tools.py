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
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#------------------------------------------MODIS Parameters---------------------------------------

P_levels = modis.transmission.pstd
freq     = [7.055309E+02, 7.196677E+02, 7.317089E+02,7.483224E+02, 9.081998E+02, 1.173198E+03,8.315149E+02]
mbands   = [31,33,34,35,36]
#Parameters to be changed ---->
year = 2016; days = 366; sc_type = 'water'

#-------------------------------------------------------------------------------------------------

#-------------------------Local Directories and Reanalysis Files----------------------------------

main_directory   = '/data/keeling/a/arkam2/c/Reanalysis/era5_files/'+str(year)+'/'
main_directory2  = '/data/keeling/a/arkam2/c/Reanalysis/standard/'
geode            = np.loadtxt('/data/gdi/c/arkam2/Reanalysis/standard/OUTPUT.DAT')
land_mask        = np.loadtxt('/data/keeling/a/arkam2/LandOceanMask/ANCILLARY_1200/data/' + \
                                'land_ocean_masks_xdeg/land_ocean_mask2_hd.asc',skiprows=6)
bias_directory   = '/data/keeling/a/arkam2/c/modis_rect/bias_correction/daily_bias_files/'+str(sc_type)+'/'+str(year)+'/' 

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
def find_coord_index_reanalysis(lat,lon,d):
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
            'latitude' : (['point'],      lat),
            'longitude': (['point'],      lon)
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
            'latitude': (['npoints'] , lat),
            'longitude': (['npoints'], lon)
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
    data_2mT12 = xr.open_dataset(os.path.join(main_directory,'era5_12UTC_'+str(year)+'_2mT.nc')).t2m[doy,:,:]
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

#------------------------------------CREATE ATMOSPHERE AT 101 MOD06 LEVELS-----------------------------------------
def find_reanalysis_vars(i, day, month, sstemp, skint, sfpres, temp, shum, geo, data_2mT, data_2mD, zouts = 'all',
                        afglms=afglms, afglmw=afglmw, afglsw=afglsw, afglss=afglss, afglt=afglt):
    
    column_lat = sstemp.latitude[i]
    column_lon = sstemp.longitude[i]

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

    T = temp.data[:,i]
    skt = skint.data[i]
    q = shum.data[:,i]
    pressures = temp.level.data.astype(np.float32)
    sfp = sfpres.data[i]/100.0
    z_levels = geo.data[:,i]/9.8e3

    t2 = data_2mT.data[i]
    d2 = data_2mD.data[i]
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
    pressures = pressures[np.where(z_levels>0.)]
    z_levels = z_levels[np.where(z_levels>0.)]
    z_levels = np.append(z_levels, np.array([0.0]))
    pressures = np.append(pressures, sfp)
    if not np.isnan(sstemp.data[i]):
        T = np.append(T, sstemp.data[i])
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
    T           = fT(np.log(P_levels)); wv = fwv(np.log(P_levels)); geo = fz(np.log(P_levels))
    return T, geo, wv, sfp, skt
#-------------------------------------------------------------------------------------------------
