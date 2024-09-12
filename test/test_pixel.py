import sys
import os
import xarray as xr
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from src.pixel import modis
import numpy as np
from misr_cth_to_pressure import height_to_log_pressure
from memory_profiler import profile
'''
 {
                'ctp': mod06.get_ctp(),
                'opt': mod06.get_opt(),
                'cphase': mod06.get_cphase(),
                'cth': mod06.get_cth(),
                'emissivity': mod06.get_emissivity()
            }
'''
@profile
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
    
    file =  '/data/keeling/a/gzhao1/f/mmcth/output/MODMISR_L2_CP_A2002337.0215_OO015738_v01.nc'
    ds  = xr.open_dataset(file)
    misr_cth = ds['misr_cth'].values
    z_prof = ds['geopotential'].values

    # Initialize the output array
    pressure_at_cth = np.empty_like(misr_cth)

    # Call the Fortran subroutine
    misr_ctp = height_to_log_pressure(z_prof, misr_cth/1000., z_prof.shape[1], z_prof.shape[2])

    mod06 = ds['mod06_product'].values 
    mod_cth = mod06[3,:,:]
    mod_cphase = mod06[2,:,:]

    # Assuming modis_heights, misr_heights, and cloud_phase are 2D arrays

    height_diff = np.abs(mod_cth - misr_cth)
    condition = (height_diff > 1000) & (mod_cphase == 3)
    # Get the indices of the selected pixels
    selected_pixels = np.where(condition)

    # mod_emm = mod06[-1,:,:]
    # mod_ctp = mod06[0,:,:]

    selected_sst = ds['sst'].values[selected_pixels[0],selected_pixels[1]] 
    selected_dew2m  = ds['dew2m'].values[selected_pixels[0],selected_pixels[1]] 
    selected_sp_prof = ds['specific_humidity'].values[:, selected_pixels[0],selected_pixels[1]] 
    selected_z_prof = z_prof[:, selected_pixels[0],selected_pixels[1]] 
    selected_surftmp  = ds['surface_pressure'].values[selected_pixels[0],selected_pixels[1]] /100 # CL
    selected_landsea = ds['landsea'].values[selected_pixels[0],selected_pixels[1]]
    selected_skintmp = ds['skint'].values[selected_pixels[0],selected_pixels[1]]
    selected_vza  = ds['vza'].values[selected_pixels[0],selected_pixels[1]]
    print('before calling tehf funcitons')
main()
# # # Initialize empty 2D arrays for the results
# # output_ctp = np.full((rows, cols), np.nan)  # Cloud Top Pressure
# # output_eca = np.full((rows, cols), np.nan)  # Effective Cloud Amount

# # # Assign the results from the selected pixels back to their original positions
# # output_ctp[selected_pixels] = ctp
# # output_eca[selected_pixels] = eca
