import sys
import os
import xarray as xr
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from src.pixel import modis
import numpy as np

from misr_cth_to_pressure import height_to_log_pressure
# from memory_profiler import profile


'''
 {
                'ctp': mod06.get_ctp(),
                'opt': mod06.get_opt(),
                'cphase': mod06.get_cphase(),
                'cth': mod06.get_cth(),
                'emissivity': mod06.get_emissivity()
            }
'''
# @profile
def main( ):
    int_press_o = np.array([4.9999999e-03, 1.6100001e-02, 3.8400002e-02, 7.6899998e-02, 1.3699999e-01,
                    2.2440000e-01, 3.4540001e-01, 5.0639999e-01, 7.1399999e-01, 9.7530001e-01,
                    1.2972000e+00, 1.6872000e+00, 2.1526000e+00, 2.7009001e+00, 3.3397999e+00,
                    4.0770001e+00, 4.9204001e+00, 5.8776002e+00, 6.9566998e+00, 8.1654997e+00,
                    9.5118999e+00, 1.1003800e+01, 1.2649200e+01, 1.4455900e+01, 1.6431801e+01,
                    1.8584700e+01, 2.0922400e+01, 2.3452600e+01, 2.6182899e+01, 2.9121000e+01,
                    3.2274399e+01, 3.5650501e+01, 3.9256599e+01, 4.3100101e+01, 4.7188202e+01,
                    5.1527802e+01, 5.6125999e+01, 6.0989498e+01, 6.6125298e+01, 7.1539803e+01,
                    7.7239601e+01, 8.3231003e+01, 8.9520401e+01, 9.6113800e+01, 1.0301720e+02,
                    1.1023660e+02, 1.1777750e+02, 1.2564560e+02, 1.3384621e+02, 1.4238480e+02,
                    1.5126640e+02, 1.6049590e+02, 1.7007840e+02, 1.8001830e+02, 1.9032030e+02,
                    2.0098869e+02, 2.1202769e+02, 2.2344150e+02, 2.3523380e+02, 2.4740849e+02,
                    2.5996909e+02, 2.7291910e+02, 2.8626169e+02, 3.0000000e+02, 3.1413690e+02,
                    3.2867529e+02, 3.4361761e+02, 3.5896649e+02, 3.7472409e+02, 3.9089261e+02,
                    4.0747379e+02, 4.2446979e+02, 4.4188190e+02, 4.5971179e+02, 4.7796069e+02,
                    4.9662979e+02, 5.1571997e+02, 5.3523218e+02, 5.5516687e+02, 5.7552478e+02,
                    5.9630621e+02, 6.1751123e+02, 6.3913977e+02, 6.6119202e+02, 6.8366730e+02,
                    7.0656543e+02, 7.2988568e+02, 7.5362750e+02, 7.7778967e+02, 8.0237140e+02,
                    8.2737128e+02, 8.5278802e+02, 8.7862012e+02, 9.0486591e+02, 9.3152362e+02,
                    9.5859113e+02, 9.8606659e+02, 1.0139476e+03, 1.0422319e+03, 1.0709170e+03,
                    1.1000000e+03])
    file =  '/data/keeling/a/gzhao1/f/mmcth/output/MODMISR_L2_CP_A2016074.0330_O086367_v01.nc'
    
    

    ds  = xr.open_dataset(file)
    
    misr_cth = ds['misr_cth'].values
    misr_cth[misr_cth<-888] = np.nan

    z_prof = ds['geopotential'].values
    
    # Initialize the output array
    pressure_at_cth = np.empty_like(misr_cth)

    # Call the Fortran subroutine
    misr_ctp = height_to_log_pressure(z_prof, misr_cth/1000., z_prof.shape[1], z_prof.shape[2])
    mod_bt = ds['bands_BT'].values  
     
    mod06 = ds['mod06_product'].values 
    mod_cth = mod06[3,:,:]
    mod_ctp = mod06[0,:,:]
    mod_cphase = mod06[2,:,:]
    
    # Get the indices of the selected pixels
    mod_bt_all_bands = np.all(mod_bt > 0, axis=0)

    height_diff = mod_cth - misr_cth
    
    condition = (height_diff > 1000) & (mod_cphase == 3) & mod_bt_all_bands
    
    if not np.any(condition) :
        pass

    selected_pixels = np.where(condition)
    
    selected_sst = ds['sst'].values[selected_pixels[0],selected_pixels[1]] 
    selected_dew2m  = ds['dew2m'].values[selected_pixels[0],selected_pixels[1]] 
    selected_landsea = ds['landsea'].values[selected_pixels[0],selected_pixels[1]]  # landsea = 11
    selected_skintmp = ds['skint'].values[selected_pixels[0],selected_pixels[1]]
     
    selected_sp_prof = ds['specific_humidity'].values[::-1, selected_pixels[0],selected_pixels[1]]/100 #Wprof - 1
    selected_z_prof = z_prof[:, selected_pixels[0],selected_pixels[1]] 
    selected_t_prof = ds['temperature'].values[::-1, selected_pixels[0],selected_pixels[1]] #tprof - 2
    selected_surfp  = ds['surface_pressure'].values[selected_pixels[0],selected_pixels[1]] /100 # psfc - 3
    selected_seasurfp  = ds['msp'].values[selected_pixels[0],selected_pixels[1]] /100 # pmsl - 4
    selected_surftemp = np.copy(selected_sst)
    selected_surftemp[selected_landsea == 1] =  selected_skintmp[selected_landsea == 1] #surftmp - 5
    selected_vza  = ds['vza'].values[selected_pixels[0],selected_pixels[1]] #view -6
    selected_bt = mod_bt[:,selected_pixels[0],selected_pixels[1]]   #tcold -7
    selected_misrcth = misr_cth[selected_pixels[0],selected_pixels[1]]   #MISR_CTH 
    selected_misrctp = misr_ctp[selected_pixels[0],selected_pixels[1]] -200  #MISR_CTP - 13
    selected_modcth = mod_cth[selected_pixels[0],selected_pixels[1]]   #MISR_CTH 
    selected_modctp = mod_ctp[selected_pixels[0],selected_pixels[1]]   #MISR_CTP
    selected_lat = ds['latitude'].values[selected_pixels[0],selected_pixels[1]] #latitude -6
    selected_lon = ds['longitude'].values[selected_pixels[0],selected_pixels[1]] #latitude -6
    
    npix = selected_sst.shape[0]  # npix should be 27944
    jday = 74  # Example jday

    
    # Before Call the Fortran subroutine
    print(f"Shape of selected_sp_prof (wprof): {selected_sp_prof.shape if selected_sp_prof is not None else 'None'}")
    print(f"Shape of selected_t_prof (tprof): {selected_t_prof.shape if selected_t_prof is not None else 'None'}")
    print(f"Shape of selected_surfp (psfc): {selected_surfp.shape if selected_surfp is not None else 'None'}")
    print(f"Shape of selected_seasurfp (pmsl): {selected_seasurfp.shape if selected_seasurfp is not None else 'None'}")
    print(f"Shape of selected_surftemp (surftmp): {selected_surftemp.shape if selected_surftemp is not None else 'None'}")
    print(f"Shape of selected_vza (view): {selected_vza.shape if selected_vza is not None else 'None'}")
    print(f"Shape of selected_bt (tcold): {selected_bt.shape if selected_bt is not None else 'None'}")
    print(f"Shape of selected_lat (rlat): {selected_lat.shape if selected_lat is not None else 'None'}")
    print(f"Shape of selected_lon (rlon): {selected_lon.shape if selected_lon is not None else 'None'}")
    print(f"Shape of selected_landsea (landsea): {selected_landsea.shape if selected_landsea is not None else 'None'}")
    print(f"Shape of selected_misrctp (misr_ctp): {selected_misrctp.shape if selected_misrctp is not None else 'None'}")    
    print(f"npix: {npix}")
    print(f"jday: {jday}")

    # pixel_idx = 10
    # selected_sp_prof = np.array([modis.transmission.wstd]).T
    # selected_t_prof = np.array([modis.transmission.tstd]).T
    # selected_surfp  = np.array([modis.transmission.pstd[-1]])
    # selected_seasurfp = selected_surfp
    # selected_surftemp = np.array([modis.transmission.tstd[-1]])
    # selected_vza  = np.array([selected_vza[pixel_idx]])
    # selected_bt = np.array([selected_bt[:,pixel_idx]]).T
    # selected_lat = np.array([selected_lat[pixel_idx]])
    # selected_lon = np.array([selected_lon[pixel_idx]])
    # selected_landsea = np.array([selected_landsea[pixel_idx]])
    # selected_misrctp = np.array([selected_misrctp[pixel_idx]])
    # npix = 1
    # Call the Fortran subroutine
    print(f"Shape of selected_sp_prof (wprof): {selected_sp_prof.shape if selected_sp_prof is not None else 'None'}")
    print(f"Shape of selected_t_prof (tprof): {selected_t_prof.shape if selected_t_prof is not None else 'None'}")
    print(f"Shape of selected_surfp (psfc): {selected_surfp.shape if selected_surfp is not None else 'None'}")
    print(f"Shape of selected_seasurfp (pmsl): {selected_seasurfp.shape if selected_seasurfp is not None else 'None'}")
    print(f"Shape of selected_surftemp (surftmp): {selected_surftemp.shape if selected_surftemp is not None else 'None'}")
    print(f"Shape of selected_vza (view): {selected_vza.shape if selected_vza is not None else 'None'}")
    print(f"Shape of selected_bt (tcold): {selected_bt.shape if selected_bt is not None else 'None'}")
    print(f"Shape of selected_lat (rlat): {selected_lat.shape if selected_lat is not None else 'None'}")
    print(f"Shape of selected_lon (rlon): {selected_lon.shape if selected_lon is not None else 'None'}")
    print(f"Shape of selected_landsea (landsea): {selected_landsea.shape if selected_landsea is not None else 'None'}")
    print(f"Shape of selected_misrctp (misr_ctp): {selected_misrctp.shape if selected_misrctp is not None else 'None'}")    
    print(f"npix: {npix}")
    print(f"jday: {jday}")
    # print(f'MODIS CTP: {mod_ctp[pixel_idx]}')

    print(f'MODIS CTH: {selected_modcth[0]}')
    print(f'MODIS CTP: {selected_modctp[0]}')
    ctp,_ = modis.process_selected_pixels(selected_sp_prof/100, selected_t_prof, selected_surfp, selected_seasurfp,
                                  selected_surftemp, selected_vza, selected_bt,selected_lat,
                                   selected_lon, selected_landsea, selected_surfp-50,jday,npix)
    # print(f'CTP Dimension: {ctp.shape}') 
    # print('before calling the funcitons')
    # index = 3  # Replace with the desired index  
    # print("SST:", selected_sst[index])
    # print("Dew 2m:", selected_dew2m[index])
    # print("Land/Sea:", selected_landsea[index])
    # print("Skin Temp:", selected_skintmp[index])
    # print("Specific Humidity Profile:", selected_sp_prof[:, index])
    # print("Z Profile:", selected_z_prof[:, index])
    # print("Temperature Profile:", selected_t_prof[:, index])
    # print("Surface Pressure (hPa):", selected_surfp[index])
    # print("Mean Sea Level Pressure (hPa):", selected_seasurfp[index])
    # print("Surface Temp:", selected_surftemp[index])
    # print("Viewing Zenith Angle:", selected_vza[index])
    # print("Brightness Temperature:", selected_bt[:, index])
    # print("MISR Cloud Top Height:", selected_misrcth[index])
    # print("MISR Cloud Top Pressure:", selected_misrctp[index])
    # print("MODIS Cloud Top Height:", selected_modcth[index])
    # print("MODIS Cloud Top Pressure:", selected_modctp[index])
    # print("Latitude", selected_lat[index])
    # print("Longitude:", selected_lon[index])


main()
 
# # # Assign the results from the selected pixels back to their original positions
# # output_ctp[selected_pixels] = ctp
# # output_eca[selected_pixels] = eca
