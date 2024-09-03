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

mod031  = '/data/keeling/a/arkam2/c/DATA_MODIS/MOD03/2016/074/MOD03.A2016074.0330.061.2017325025726.hdf'
mod061  = '/data/keeling/a/arkam2/c/DATA_MODIS/MOD06/2016/074/MOD06_L2.A2016074.0330.061.2017325190341.hdf'
mod021  = '/data/keeling/a/arkam2/c/DATA_MODIS/MOD021KM/MOD021KM.A2016074.0330.061.2017325033257.hdf'
misr    = '/data/keeling/a/arkam2/b/MISR/2016/2016.03.14/MISR_AM1_TC_CLOUD_P126_O086367_F01_0001.hdf'

#------------------------------------------MODIS Parameters---------------------------------------
#NOTE TO SELF: MODIFY YEAR

P_levels = modis.transmission.pstd
#Parameters to be changed ---->
year = 2016; mod_band = [31,36,35,34,33]

#-------------------------------------------------------------------------------------------------

#-------------------------Local Directories and Reanalysis Files----------------------------------

main_directory   = '/data/keeling/a/arkam2/c/Reanalysis/era5_files/'+str(year)+'/'
main_directory2  = '/data/keeling/a/arkam2/c/Reanalysis/standard/'
geode            = np.loadtxt('/data/gdi/c/arkam2/Reanalysis/standard/OUTPUT.DAT')
land_mask        = np.loadtxt('/data/keeling/a/arkam2/LandOceanMask/ANCILLARY_1200/data/' + \
                                'land_ocean_masks_xdeg/land_ocean_mask2_hd.asc',skiprows=6)

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
    dist, indexes = mytree.query(pt, distance_upper_bound=np.inf)
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
            'level':     (['level'], P_levels[np.where(P_levels>=2.)]),
            'latitude' : (['latitude'],      lat),
            'longitude': (['longitude'],      lon)
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
            'latitude':  (['latitude'] , lat),
            'longitude': (['longitude'], lon)
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
def find_reanalysis_vars(column_lat, column_lon, month, zouts = 'all',
                        afglms=afglms, afglmw=afglmw, afglsw=afglsw, afglss=afglss, afglt=afglt):

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


    high_col    = profile.iloc[np.logical_and(profile['p(mb)'].values<2., profile['p(mb)'].values>0.0005)]
    low_col     = profile.iloc[profile['p(mb)'].values>=900.]
    hpressures  = high_col['p(mb)']  ; hz_levels = high_col['z(km)']
    hh2o        = high_col['h2o(cm-3)'].values/high_col['air(cm-3)'].values*18.016/28.966*1e3
    lpressures  = low_col['p(mb)']   ; lz_levels = low_col['z(km)']
    lh2o        = low_col['h2o(cm-3)'].values/low_col['air(cm-3)'].values*18.016/28.966*1e3
    hT          = high_col['T(K)']   ; lT        = low_col['T(K)']

    fT          = interp1d(np.log(hpressures),hT,fill_value='extrapolate')
    fwv         = interp1d(np.log(hpressures),hh2o,fill_value='extrapolate')
    fz          = interp1d(np.log(hpressures),hz_levels,fill_value='extrapolate')
    Th          = fT(np.log(P_levels[np.where(P_levels<2.)]))
    wvh         = fwv(np.log(P_levels[np.where(P_levels<2.)]))
    geoh        = fz(np.log(P_levels[np.where(P_levels<2.)]))
    
    fT          = interp1d(np.log(lpressures),lT,fill_value='extrapolate')
    fwv         = interp1d(np.log(lpressures),lh2o,fill_value='extrapolate')
    fz          = interp1d(np.log(lpressures),lz_levels,fill_value='extrapolate')
    Tl          = fT(np.log(P_levels[np.where(P_levels>=1000.)]))
    wvl         = fwv(np.log(P_levels[np.where(P_levels>=1000.)]))
    geol        = fz(np.log(P_levels[np.where(P_levels>=1000.)]))
    
    return Th, geoh, wvh, Tl, geol, wvl
#-------------------------------------------------------------------------------------------------

                          
def do_kdtree(A, pt):
    mytree        = scipy.spatial.cKDTree(A,balanced_tree=False,compact_nodes=False)
    dist, indexes = mytree.query(pt,distance_upper_bound=np.inf)
    return indexes

def get_modis(hdf, hdf_geo):
    lati    = hdf_geo.select('Latitude');         latit      = lati[:,:]
    long    = hdf_geo.select('Longitude');         longit    = long[:,:]
    mod_lat = latit.ravel()
    mod_lon = longit.ravel()
    vza     = hdf_geo.select('SensorZenith');      modis_vza = vza[:,:]
    mod_vza = modis_vza.ravel();                     mod_vza = mod_vza*0.01
    ls      = hdf_geo.select('Land/SeaMask')[:,:];      ls   = ls.astype(float)
    ls[np.where(ls<5.)] = 1
    ls[np.where(ls>5.)] = 0; ls = ls.ravel()
    data_id = hdf.select('EV_1KM_Emissive');           attrs = data_id.attributes(full=1)
    offsets = attrs["radiance_offsets"][0]; scales = attrs["radiance_scales"][0]
#     data_29 = np.append(data_29, np.ravel((data_selected_id[8,:, :] - offsets[8])  * scales[8])) 
#     data_32 = np.append(data_32, np.ravel((data_selected_id[11,:,:] - offsets[11]) * scales[11]))
    data_31 = np.ravel((data_id[10,:,:] - offsets[10]) * scales[10])
    data_33 = np.ravel((data_id[12,:,:] - offsets[12]) * scales[12])
    data_34 = np.ravel((data_id[13,:,:] - offsets[13]) * scales[13])
    data_35 = np.ravel((data_id[14,:,:] - offsets[14]) * scales[14])
    data_36 = np.ravel((data_id[15,:,:] - offsets[15]) * scales[15])
    
    return mod_lat, mod_lon, mod_vza, ls, data_31, data_33, data_34, data_35, data_36

def get_mod06(hdf_06):
    height     = hdf_06.select('cloud_top_height_1km')[:,:]/1e3
    press      = hdf_06.select('cloud_top_pressure_1km')[:,:]/10
    method     = hdf_06.select('cloud_top_method_1km')[:,:]
    
    return height, press, method

def get_misr(hdf, start, end):
    mt          = Mtk.MtkFile(hdf)
    path_number = mt.attr_get('Path_number')
    r           = Mtk.MtkRegion(path_number, start, end) #60,74
    d           = mt.grid('Stereo_1.1_km').field('CloudTopHeight').read(r)
    d           = d.data() * 0.001

    # Read geolocation dataset
    map_info    = Mtk.MtkRegion(path_number, start, end).snap_to_grid(path_number, 1100)
    latlon      = np.asarray(map_info.create_latlon())
    misrlat     = latlon[0,:,:]
    misrlon     = latlon[1,:,:]
    
#     d[np.where(d<0)]                                            = np.nan
#     d[np.where(misrlat==-200.0) and np.where(misrlon==-200.0)]  = np.nan
    
    return misrlat, misrlon, d

hdf        = SD(mod021, SDC.READ)
hdf_geo    = SD(mod031, SDC.READ)
hdf_06     = SD(mod061, SDC.READ)
mod_lat, mod_lon, mod_vza, ls, data_31, \
data_33, data_34, data_35, data_36  = get_modis(hdf, hdf_geo)
mod_cth, mod_pres, mod_method       = get_mod06(hdf_06)
misrlat, misrlon, d                 = get_misr(misr, 64, 66)
hr                                  = float(mod031.split('/')[10].split('.')[2][0:2])+ \
                                            float(mod031.split('/')[10].split('.')[2][2:])/60
doy                                 = int(mod031.split('/')[10].split('.')[1][5:8])
month1,sst1,skin1,sfp1,temp1,shum1, \
    geo1,data_2mT1,data_2mD1,sst2, \
    skin2,sfp2,temp2,shum2,geo2,  \
               data_2mT2,data_2mD2  = find_reanalysis(hr,doy)
latitude                            = sst1.latitude.values 
longitude                           = sst1.longitude.values

start = time.time()

T1 = interp3d(temp1,latitude,longitude); T2 = interp3d(temp2,latitude,longitude); print(time.time()-start)
w1 = interp3d(shum1,latitude,longitude); w2 = interp3d(shum2,latitude,longitude); print(time.time()-start)
g1 = interp3d(geo1,latitude,longitude);  g2 = interp3d(geo2,latitude,longitude) ; print(time.time()-start)
s1 = interp2d(skin1,latitude,longitude); s2 = interp2d(skin2,latitude,longitude); print(time.time()-start)
p1 = interp2d(sfp1,latitude,longitude);  p2 = interp2d(sfp2,latitude,longitude) ; print(time.time()-start)

tpr = interp_time(T1, T2, 0, 6, hr)
wvr = interp_time(w1, w2, 0, 6, hr)
ps  = interp_time(p1, p2, 0, 6, hr)
geo = interp_time(g1, g2, 0, 6, hr)
sft = interp_time(s1, s2, 0, 6, hr)

print(time.time()-start)

clat                         = np.mean(misrlat[np.where(misrlat!=-200.)]) 
clon                         = np.mean(misrlon[np.where(misrlon!=-200.)])

Th, geoh, wvh, Tl, geol, wvl = find_reanalysis_vars(clat, clon, month1)
T_high                       = np.zeros((len(Th)  , np.shape(latitude)[0], np.shape(longitude)[0]))
z_high                       = np.zeros((len(geoh), np.shape(latitude)[0], np.shape(longitude)[0]))
w_high                       = np.zeros((len(wvh) , np.shape(latitude)[0], np.shape(longitude)[0]))
T_high[:,:,:]                = Th[:,np.newaxis,np.newaxis]
z_high[:,:,:]                = geoh[:,np.newaxis,np.newaxis]
w_high[:,:,:]                = wvh[:,np.newaxis,np.newaxis]

T_low                        = np.zeros((len(Tl)  , np.shape(latitude)[0], np.shape(longitude)[0]))
z_low                        = np.zeros((len(geol), np.shape(latitude)[0], np.shape(longitude)[0]))
w_low                        = np.zeros((len(wvl) , np.shape(latitude)[0], np.shape(longitude)[0]))
T_low[:,:,:]                 = Tl[:,np.newaxis,np.newaxis]
z_low[:,:,:]                 = geol[:,np.newaxis,np.newaxis]
w_low[:,:,:]                 = wvl[:,np.newaxis,np.newaxis]

tpr = interp_time(T1, T2, 0, 6, hr); tpr = np.append(np.zeros((12,721,1440)),tpr,axis=0)
wvr = interp_time(w1, w2, 0, 6, hr); wvr = np.append(np.zeros((12,721,1440)),wvr*1e3,axis=0)
ps  = interp_time(p1, p2, 0, 6, hr)
geo = interp_time(g1, g2, 0, 6, hr); geo = np.append(np.zeros((12,721,1440)),geo/9.8e3,axis=0)
sft = interp_time(s1, s2, 0, 6, hr)

tpr[0:12,:,:] = T_high       ; tpr[97:101,:,:] = T_low
wvr[0:12,:,:] = w_high       ; wvr[97:101,:,:] = w_low
geo[0:12,:,:] = z_high       ; geo[97:101,:,:] = z_low

#let's collocate MISR and MODIS

misr_points    = ([misrlat.ravel(), misrlon.ravel()]); dims = len(misr_points)
modis_points   = np.dstack([mod_lat,mod_lon])[0]
results        = do_kdtree(modis_points,np.transpose(misr_points))

rad            = np.zeros((np.shape(results)[0], 5))
rad[:,0]       = data_36[results] ; rad[:,1] = data_35[results]   ; rad[:,2] = data_34[results]
rad[:,3]       = data_33[results] ; rad[:,4] = data_31[results]
vza            = mod_vza[results] ; landsea  = ls[results]
mod_ctp        = mod_pres.ravel()[results];
method         = mod_method.ravel()[results];
mod_ht         = mod_cth.ravel()[results]

#let's collocate MISR and Reanalysis and find MISR CTP by interpolation

era5_points    = [[i,j] for i,j in zip(latitude,longitude)] 
results1       = do_kdtree(era5_points,np.transpose(misr_points))

tpr  =   tpr.reshape(tpr.shape[0],-1).T[results1]
wvr  =   wvr.reshape(wvr.shape[0],-1).T[results1]
geo  =   geo.reshape(geo.shape[0],-1).T[results1]
ps   =   np.ravel(ps)[results1]
sft  =   np.ravel(sft)[results1]

misr_ht                                              = d.ravel()
misr_ctp                                             = P_levels[np.argmin(np.abs(geo - misr_ht[:,np.newaxis]), axis=1)]
mm_flag                                              = np.ones(len(misr_ht))*-1
mm_flag[np.where(np.logical_and(misrlat.ravel!=-200.,
                             misrlon.ravel!=-200.))] = 0

mod_tmp                                              = mod_ht
mod_tmp[np.where(mod_tmp==-0.999)]                   = -9.999
ind                                                  = np.where(np.logical_and(np.logical_and((mod_tmp - misr_ht)>1,
                                                                                              method < 5), 
                                                                 np.logical_and(mod_tmp!=-9.999, misr_ht!=-9.999)))
mm_flag[ind[0]]                                      = 1
dims                                                 = len(mm_flag)
met_date                                             = np.array([2016,month1,14,hr])

tcld                                                 = np.zeros((dims,5))
for kban in range(len(mod_band)):
    tcld[:,kban]                                     = [modis.planck_bright.modis_bright(rad=rad[i,kban],band=mod_band[kban],units=1) 
                                                                                                     for i in range(dims)]

print('NOW BEGINNING LOOP')                                                                          
ctp, ec = modis.mm_cth(wvr,tpr,ps,ps,sft,vza,tcld,met_date,misrlat.ravel(),misrlon.ravel(),landsea,misr_ctp,mm_flag,dims=dims)
np.savetxt('loop_36_35.txt', ctp[:,0].reshape(misrlat.shape))
np.savetxt('loop_34_33.txt', ctp[:,2].reshape(misrlat.shape))