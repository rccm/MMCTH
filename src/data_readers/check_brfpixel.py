import numpy as np
import glob
import misr_read
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')

file_dir = '/data/gdi/d/MISRCM/cases/74901/'
cams = ['AN','CF','BF','AF','AN','AA','BA','CA','DA']
bands = ['blue', 'green', 'red','nir']

results = np.zeros((len(bands), len(cams)))
for i, band in enumerate(bands):
    for j, cam in enumerate(cams):
        file = f'{file_dir}MISR_AM1_GRP_TERRAIN_GM_P041_O074901_{cam}_F03_0024.hdf'
        a = misr_read.MISRGranule(file)
        if (band == 'red' or cam == 'AN'):
            b = a.get_brf(band, [62, 62])
            rad = a.get_rad(band, [62, 62])
            print(rad[0,381-1, 1384-1])
            results[i, j] = b[0, 381, 1384]
        else:
            b = a.get_brf(band, [62, 62])
            results[i, j] = b[0, 95, 346]
        break
    break

# import MisrToolkit as Mtk
# import numpy as np

# file_dir = '/data/gdi/d/MISRCM/cases/74901/'
# cams = ['AN','CF','BF','AF','AN','AA','BA','CA','DA']
# bands = ['blue', 'green', 'red','nir']

# results = np.zeros((len(bands), len(cams)))
# for i, band in enumerate(bands):
#     for j, cam in enumerate(cams):
#         file = f'{file_dir}MISR_AM1_GRP_TERRAIN_GM_P041_O074901_{cam}_F03_0024.hdf'
#         f = Mtk.MtkFile(file).grid('BlueBand').field('Blue BRF')
#         r = Mtk.MtkRegion(41, 62, 62)
#         c = f.read(r)
#         c = c.data()
#         print(c[381-1, 1384-1])
#         break
# results = np.zeros((len(bands), len(cams)))
# for i, band in enumerate(bands):
#     for j, cam in enumerate(cams):
#         file = f'{file_dir}MISR_AM1_GRP_TERRAIN_GM_P041_O074901_{cam}_F03_0024.hdf'
#         a = misr_read.MISRGranule(file)
#         if (band == 'red' or cam == 'AN'):
#             b = a.get_brf(band, [62, 62])
#             results[i, j] = b[0, 381-1, 1384-1]
#         else:
#             b = a.get_brf(band, [62, 62])
#             results[i, j] = b[0, 95-1, 346-1]

# np.savetxt('misr_brf_results.txt', results, fmt='%.6f')
# fig, ax = plt.subplots()
# data = np.loadtxt('misr_brf_results.txt')

# # Plot each row of data for each band with specified colors and line styles
# colors = ['blue', 'green', 'red', 'black']
# wavelengths = [446, 558, 672, 866] 

# linestyles = ['-', '-', '-', '-']
# for i, (band, color, linestyle) in enumerate(zip(wavelengths, colors, linestyles)):
#     ax.plot(cams[::-1], data[i][::-1], label=f'{band}', color=color, linestyle=linestyle)

# # Add labels and legend inside the plot
# ax.set_xlabel('Cameras')
# ax.set_ylabel('BRF')
# ax.set_ylim(0, 0.7) 
# ax.legend(loc='upper right')
# fig.savefig('brf.png')


import MisrToolkit as Mtk
import numpy as np

file_dir = '/data/gdi/d/MISRCM/cases/74901/'
cams = ['AN','CF','BF','AF','AN','AA','BA','CA','DA']
bands = ['blue', 'green', 'red','nir']

results = np.zeros((len(bands), len(cams)))
for i, band in enumerate(bands):
    for j, cam in enumerate(cams):
        file = f'{file_dir}MISR_AM1_GRP_TERRAIN_GM_P041_O074901_{cam}_F03_0024.hdf'
        f = Mtk.MtkFile(file).grid('BlueBand').field('Blue Radiance')
        r = Mtk.MtkRegion(41, 62, 62)
        c = f.read(r)
        c = c.data()
        print(c[381-1, 1384-1])