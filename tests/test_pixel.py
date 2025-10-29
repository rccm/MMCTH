import sys
import os
import xarray as xr
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from src.pixel import modis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

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
PSTD = np.array([
     0.0050,   0.0161,   0.0384,   0.0769,   0.1370,
     0.2244,   0.3454,   0.5064,   0.7140,   0.9753,   1.2972,
     1.6872,   2.1526,   2.7009,   3.3398,   4.0770,   4.9204,
     5.8776,   6.9567,   8.1655,   9.5119,  11.0038,  12.6492,
    14.4559,  16.4318,  18.5847,  20.9224,  23.4526,  26.1829,
    29.1210,  32.2744,  35.6505,  39.2566,  43.1001,  47.1882,
    51.5278,  56.1260,  60.9895,  66.1253,  71.5398,  77.2396,
    83.2310,  89.5204,  96.1138, 103.0172, 110.2366, 117.7775,
   125.6456, 133.8462, 142.3848, 151.2664, 160.4959, 170.0784,
   180.0183, 190.3203, 200.9887, 212.0277, 223.4415, 235.2338,
   247.4085, 259.9691, 272.9191, 286.2617, 300.0000, 314.1369,
   328.6753, 343.6176, 358.9665, 374.7241, 390.8926, 407.4738,
   424.4698, 441.8819, 459.7118, 477.9607, 496.6298, 515.7200,
   535.2322, 555.1669, 575.5248, 596.3062, 617.5112, 639.1398,
   661.1920, 683.6673, 706.5654, 729.8857, 753.6275, 777.7897,
   802.3714, 827.3713, 852.7880, 878.6201, 904.8659, 931.5236,
   958.5911, 986.0666,1013.9476,1042.2319,1070.9170,1100.0000
], dtype=np.float64)

# 2. Standard temperature profile (K)
TSTD = np.array([
    190.19,203.65,215.30,226.87,237.83,247.50,256.03,263.48,267.09,270.37,
    266.42,261.56,256.40,251.69,247.32,243.27,239.56,236.07,232.76,230.67,
    228.71,227.35,226.29,225.28,224.41,223.61,222.85,222.12,221.42,220.73,
    220.07,219.44,218.82,218.23,217.65,217.18,216.91,216.70,216.70,216.70,
    216.70,216.70,216.70,216.70,216.70,216.70,216.70,216.70,216.70,216.70,
    216.70,216.70,216.71,216.71,216.72,216.81,217.80,218.77,219.72,220.66,
    222.51,224.57,226.59,228.58,230.61,232.61,234.57,236.53,238.48,240.40,
    242.31,244.21,246.09,247.94,249.78,251.62,253.45,255.26,257.04,258.80,
    260.55,262.28,264.02,265.73,267.42,269.09,270.77,272.43,274.06,275.70,
    277.32,278.92,280.51,282.08,283.64,285.20,286.74,288.25,289.75,291.22,
    292.68
], dtype=np.float32)

# 3. Standard specific-humidity profile (g kg-1)
WSTD = np.array([
    0.001,0.001,0.002,0.003,0.003,
    0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.03,
    0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,
    0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,
    0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.003,0.004,0.004,
    0.005,0.005,0.007,0.009,0.011,0.012,0.014,0.020,0.025,0.030,0.035,0.047,
    0.061,0.075,0.089,0.126,0.162,0.197,0.235,0.273,0.310,0.356,0.410,0.471,
    0.535,0.601,0.684,0.784,0.886,0.987,1.094,1.225,1.353,1.519,1.686,1.852,
    2.036,2.267,2.496,2.721,2.947,3.170,3.391,3.621,3.848,4.084,4.333,4.579,
    4.822,5.061,5.296,5.528
], dtype=np.float32)

# 4. Standard ozone-mixing-ratio profile (ppmv)
OSTD = np.array([
    0.47330,0.27695,0.28678,0.51816,0.83229,1.18466,1.69647,2.16633,3.00338,
    3.76287,4.75054,5.61330,6.33914,7.03675,7.50525,7.75612,7.81607,7.69626,
    7.56605,7.28440,7.01002,6.72722,6.44629,6.17714,5.92914,5.69481,5.47387,
    5.26813,5.01252,4.68941,4.35141,4.01425,3.68771,3.37116,3.06407,2.77294,
    2.50321,2.24098,1.98592,1.74840,1.54451,1.34582,1.17824,1.02513,0.89358,
    0.78844,0.69683,0.62654,0.55781,0.50380,0.45515,0.42037,0.38632,0.35297,
    0.32029,0.28832,0.25756,0.22739,0.19780,0.16877,0.14901,0.13190,0.11511,
    0.09861,0.08818,0.07793,0.06786,0.06146,0.05768,0.05396,0.05071,0.04803,
    0.04548,0.04301,0.04081,0.03983,0.03883,0.03783,0.03685,0.03588,0.03491,
    0.03395,0.03368,0.03349,0.03331,0.03313,0.03292,0.03271,0.03251,0.03190,
    0.03126,0.03062,0.02990,0.02918,0.02850,0.02785,0.02721,0.02658,0.02596,
    0.02579,0.02579
], dtype=np.float32)
PLEV = np.array([4.9999999e-03, 1.6100001e-02, 3.8400002e-02, 7.6899998e-02, 1.3699999e-01,
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
# @profile


def plot_profile(wprof, tprof,
                 psfc, surftmp, misr_ctp, modis_ctp, view,
                 fname='input_profile.png',
                 pstd=PSTD, wstd=WSTD, tstd=TSTD,plev=PLEV):
    """
    Draw ERA-5 & standard humidity / temperature profiles
    and mark surface pressure, MISR CTP and MODIS CTP levels.
    """

    def _add_level_lines(ax):
        """Horizontal lines + text labels at the three key pressure levels."""
        x_max = ax.get_xlim()[1]
        levels = [(psfc,   'Surface P'),
                  (misr_ctp,  'MISR CTP'),
                  (modis_ctp, 'MODIS CTP')]
        for y, lab in levels:
            ax.axhline(y, ls='--', lw=0.9, color='k')
            ax.text(x_max*0.98, y*1.03, f'{lab}: {y:.0f} hPa',
                    va='bottom', ha='right', fontsize=8)
        ax.text(x_max*0.98, 0.1, f'VZA = {str(view)}',ha='right', fontsize=12)

    # ------------------------------------------------------------------
    # 1. Water vapour mixing-ratio profile
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(wprof[wprof>0], plev[wprof>0], label='ERA5 + Standard Profile ')
    ax.plot(wstd, plev, label='MODIS Standard Profile')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_xlabel('Water mixing ratio [g kg$^{-1}$]')
    ax.set_ylabel('Pressure [hPa]')
    ax.grid(True, which='both', ls=':')
    _add_level_lines(ax)
    ax.legend()
    fig.tight_layout()
    fig.savefig('w_' + fname, dpi=300)
    plt.close(fig)

    # ------------------------------------------------------------------
    # 2. Temperature profile
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(tprof[tprof>0], pstd[tprof>0], label='ERA5')
    ax.plot(tstd,  pstd, label='Standard')
    ax.set_xscale('linear')
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Pressure [hPa]')
    ax.grid(True, which='both', ls=':')
    _add_level_lines(ax)
    ax.legend()
    fig.tight_layout()
    fig.savefig('t_' + fname, dpi=300)
    plt.close(fig)

def main( ):
    
    # file =  '/data/keeling/a/gzhao1/f/mmcth/output/MODMISR_L2_CP_A2016074.0330_O086367_v01.nc'
    file =  '/data/keeling/a/gzhao1/f/mmcth/output/MODMISR_L2_CP_A2016074.0330_O086367_debug.nc'    

    ds  = xr.open_dataset(file)
    
    misr_cth = ds['misr_cth'].values
    misr_cth[misr_cth<-888] = np.nan
    z_prof = ds['geopotential_exp'].values
    
    # Initialize the output array
    pressure_at_cth = np.empty_like(misr_cth)

    # Call the Fortran subroutine
    misr_ctp = height_to_log_pressure(z_prof, misr_cth/1000., z_prof.shape[1], z_prof.shape[2])
    mod_bt_rad = ds['bands_BT'].values  
     
    mod06 = ds['mod06_product'].values 
    mod_cth = mod06[3,:,:]
    mod_ctp = mod06[0,:,:]
    mod_cphase = mod06[2,:,:]
    
    # Get the indices of the selected pixels

    mod_bt = mod_bt_rad[0:5,:,:]
    mod_rad = mod_bt_rad[5:,:,:]

    mod_bt_all_bands = np.all(mod_bt > 0, axis=0)

    height_diff = mod_cth - misr_cth
    
    condition = (height_diff > 1000) & (mod_cphase == 3) & mod_bt_all_bands & (misr_cth >0)
    
    if not np.any(condition) :
        pass

    selected_pixels = np.where(condition)
    
    selected_sst = ds['sst'].values[selected_pixels[0],selected_pixels[1]] 
    selected_dew2m  = ds['swmr'].values[selected_pixels[0],selected_pixels[1]] 
    selected_landsea = ds['landsea'].values[selected_pixels[0],selected_pixels[1]]  # landsea = 11
    selected_skintmp = ds['skint'].values[selected_pixels[0],selected_pixels[1]]
    selected_sp_prof = ds['mixingratio_exp'].values[:, selected_pixels[0],selected_pixels[1]] #Wprof - 1
    selected_z_prof = z_prof[:, selected_pixels[0],selected_pixels[1]] 
    selected_t_prof = ds['temperature_exp'].values[:, selected_pixels[0],selected_pixels[1]] #tprof - 2
    selected_surfp  = ds['surface_pressure'].values[selected_pixels[0],selected_pixels[1]] /100 # psfc - 3
    selected_seasurfp  = ds['msp'].values[selected_pixels[0],selected_pixels[1]] /100 # pmsl - 4
    selected_surftemp = np.copy(selected_sst)
    selected_surftemp[selected_landsea == 1] =  selected_skintmp[selected_landsea == 1] #surftmp - 5
    selected_vza  = ds['vza'].values[selected_pixels[0],selected_pixels[1]] #view -6
    selected_bt = mod_bt[:,selected_pixels[0],selected_pixels[1]]   #tcold -7
    selected_misrcth = misr_cth[selected_pixels[0],selected_pixels[1]]   #MISR_CTH 
    selected_misrctp = misr_ctp[selected_pixels[0],selected_pixels[1]] #MISR_CTP - 13
    selected_modcth = mod_cth[selected_pixels[0],selected_pixels[1]]   #MISR_CTH 
    selected_modctp = mod_ctp[selected_pixels[0],selected_pixels[1]]   #MISR_CTP
    selected_lat = ds['latitude'].values[selected_pixels[0],selected_pixels[1]] #latitude -6
    selected_lon = ds['longitude'].values[selected_pixels[0],selected_pixels[1]] #latitude -6
    selected_rad = mod_rad[:,selected_pixels[0],selected_pixels[1]]  

    npix = selected_sst.shape[0]  # npix should be 27944
    jday = 74  # Example jday

     
  
  

    pix = 3 # this one works
    wprof_pix = selected_sp_prof[:, pix]
    tprof_pix = selected_t_prof[:, pix]
    psfc_pix = selected_surfp[pix]
    pmsl_pix = selected_seasurfp[pix]
    surftmp_pix = selected_surftemp[pix]
    view_pix = selected_vza[pix]
    tcold_pix = selected_bt[:, pix]
    trad_pix = selected_rad[1:, pix]
    rlat_pix = selected_lat[pix]
    rlon_pix = selected_lon[pix]
    landsea_pix = selected_landsea[pix]
    misr_ctp_pix = selected_misrctp[pix]    
    met_date = np.array([2014,12,12,0],dtype='int32').reshape(4) 


    print("MISR Cloud Top Pressure:", selected_misrctp[pix])
    print("MODIS Cloud Top Pressure:", selected_modctp[pix])
    print("Brightness Temperature:", selected_bt[:, pix])
    print("Radiance :", selected_rad[1:, pix]) 

    # plot_profile(wprof_pix,tprof_pix, psfc_pix, surftmp_pix, misr_ctp_pix,selected_modctp[pix],selected_vza[pix])
    print(trad_pix)
    ctp,cee,copt, pm = modis.co2cld_onepixel_misr(
        (wprof_pix),
        (tprof_pix), 
        psfc_pix, 
        psfc_pix, 
        surftmp_pix, 
        view_pix,
        trad_pix,  
        met_date, 
        rlat_pix, 
        rlon_pix,
        landsea_pix,
        misr_ctp_pix
        )
     
    ctp, cee,opt,pmask  = modis.process_selected_pixels(
        selected_sp_prof,
        selected_t_prof,
        selected_surfp, 
        selected_seasurfp,
        selected_surftemp, 
        selected_vza, 
        selected_rad[1:,:], 
        selected_lat,
        selected_lon,
        selected_landsea,
        selected_misrctp,
        met_date,
        selected_vza.shape[0],
    )
    print(ctp[ctp>0])
main()
 
# # # Assign the results from the selected pixels back to their original positions
# # output_ctp[selected_pixels] = ctp
# # output_eca[selected_pixels] = eca
