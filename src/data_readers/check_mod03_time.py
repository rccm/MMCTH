from datetime import datetime, timezone, timedelta
import numpy as np
from mod_read import MOD03Granule

mod03 = MOD03Granule(
    "/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2015/176/MOD03.A2015176.1650.061.2017321001513.hdf"
)

from datetime import datetime, timezone, timedelta
import numpy as np

mod03 = MOD03Granule(
    "/data/gdi/satellite/TerraDataArchive/MODIS/MOD03/2015/176/MOD03.A2015176.1650.061.2017321001513.hdf"
)

time2d = mod03.get_time_2d_since_2000(use_center_time=True, fast=True)
time1d = mod03.get_time_1d_since_2000(use_center_time=True)
ev_frames = mod03.get_ev_frames()

epoch2000 = datetime(2000, 1, 1, tzinfo=timezone.utc)

nscans = ev_frames.size
last_scan = nscans - 1
last_row0 = last_scan * 10
last_row9 = last_row0 + 9
last_valid_col = int(ev_frames[last_scan]) - 1

print("nscans =", nscans)
print("last scan index =", last_scan)
print("last scan rows =", last_row0, "to", last_row9)
print("last valid column in last scan =", last_valid_col)

print("\nLast scan, first 5 frame times from row 0 of that scan:")
for x in time2d[last_row0, 0:5]:
    print(epoch2000 + timedelta(seconds=float(x)))

print("\nLast scan, last 5 valid frame times from row 0 of that scan:")
for x in time2d[last_row0, last_valid_col-4:last_valid_col+1]:
    print(epoch2000 + timedelta(seconds=float(x)))

print("\nLast scan, last 5 valid frame times from row 9 of that scan:")
for x in time2d[last_row9, last_valid_col-4:last_valid_col+1]:
    print(epoch2000 + timedelta(seconds=float(x)))

print("\nrows last_row0 and last_row9 identical for valid frames?",
      np.allclose(
          time2d[last_row0, :last_valid_col+1],
          time2d[last_row9, :last_valid_col+1],
          equal_nan=True
      ))

print("\n1D representative time for last scan rows:")
for i in range(last_row0, last_row0 + 3):
    print(i, epoch2000 + timedelta(seconds=float(time1d[i])))

print("\nLast scan representative 1D row time equals scan-center time for all 10 rows?",
      np.allclose(time1d[last_row0:last_row0+10], time1d[last_row0], equal_nan=True))