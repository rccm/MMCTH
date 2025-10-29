import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('Agg')


int_press = np.array([4.9999999e-03, 1.6100001e-02, 3.8400002e-02, 7.6899998e-02, 1.3699999e-01,
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

or_press = np.array([1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000], dtype=int)



plev  = 101        # vertical levels
nbct  = 5          # your band count
fname = "taup_2014_346.bin"
bands =[36,35,34,33,31]

taup = np.fromfile(fname, dtype=np.float32).reshape((plev, nbct), order='F')              # Fortran order
print(np.shape(taup))
plt.figure(figsize=(8, 6))
for k in range(nbct) :
    plt.plot(taup[:,k], int_press,label=f'Band {bands[k]}')
plt.yscale('log') 
plt.gca().invert_yaxis()  # Invert y-axis to show higher pressure levels at the bottom
plt.xlabel('Temperature')
plt.ylabel('Transmittance')
plt.legend()


# # Load the NetCDF file
# file_path = '/data/keeling/a/gzhao1/f/mmcth/output/check_org.nc'
# ds = xr.open_dataset(file_path)

# # Choose a specific latitude and longitude for the profile (e.g., central grid point)

# latitude_idx = ds.latitude.size // 2  # central latitude index
# longitude_idx = ds.longitude.size // 2  # central longitude index

# print(latitude_idx,longitude_idx)
# # Extract the 't_org' data at the specified latitude and longitude across all pressure levels
# t_org_profile = ds['t_org'][:, latitude_idx, longitude_idx]
# w_org_profile = ds['q_org'][:, latitude_idx, longitude_idx]

# file_path = '/data/keeling/a/gzhao1/f/mmcth/output/check_int.nc'

# ds_int = xr.open_dataset(file_path)
# t_int_profile = ds_int['temperature'][:, latitude_idx, longitude_idx]
# w_int_profile = ds_int['mixingratio'][:, latitude_idx, longitude_idx]


# org_file_path = '/data/keeling/a/gzhao1/f/mmcth/output/check_single.nc'
# dsor = xr.open_dataset(org_file_path)
# sp = dsor['surface_pressure'][latitude_idx, longitude_idx]
# print(f'surface pressure: {sp}')

# # Plotting the profile of 't_org' across the 37 pressure levels
# plt.figure(figsize=(8, 6))
# plt.plot(t_org_profile, or_press, marker='o', color='blue',label='37 ERA5 Pressure Levels')
# plt.plot(t_int_profile[t_int_profile>0], int_press[t_int_profile>0], marker='x', color='red', label='101 Pressure Levels')
# plt.yscale('log') 
# plt.gca().invert_yaxis()  # Invert y-axis to show higher pressure levels at the bottom
# plt.xlabel('Temperature')
# plt.ylabel('Pressure Level')
# plt.title(f'Temperature Profile at Central Latitude and Longitude')
# plt.legend()
# plt.grid(True)
plt.savefig('transmittanceprofile.png')
 