
from pyhdf.SD import SD, SDC
from pyhdf.HDF import HDF, HC
from pyhdf.VS import *
import numpy as np
import re
import os
from datetime import datetime,timezone, timedelta 


class MODL1Granule():
    def __init__(self, mod_l1b_fullpath,destripe_flag=True):
        self.fullpath = mod_l1b_fullpath  # Use the passed full path
        self.destripe_flag = destripe_flag
        if os.path.isfile(self.fullpath):
            self.hdf_file = SD(self.fullpath)
        else:
            print(f'Corrupted file: {self.fullpath}')    
            raise FileNotFoundError("File not found: " + self.fullpath)
    def get_id(self):
        match = re.search(r"\.A(\d{7})\.(\d{4})\.", self.fullpath)  # Use self.fullpath to extract date/time
        if match:
            return match.group(0)
    def get_datetime(self):
        match = re.search(r"\.A(\d{7})\.(\d{4})\.", self.fullpath)  # Use self.fullpath to extract date/time
        if match:
            julian_date = match.group(1)
            time_str = match.group(2)
            date = datetime.strptime(julian_date, "%Y%j")
            time = datetime.strptime(time_str, "%H%M").time()
            return julian_date, date, time
        else:
            print("No match found in the filename.")
            return None, None, None
    def get_yr(self):
        julian_date, date, time = self.get_datetime()
        return julian_date[0:4]
    def get_doy(self):
        julian_date, date, time = self.get_datetime()
        return julian_date[4:]
    def get_month(self):
        julian_date, date, time = self.get_datetime()
        return date[4:6]
    def get_time(self):
        julian_date, date, time = self.get_datetime()
        return time
    def get_coremetadata(self):
        hdf = HDF(self.fullpath, HC.READ)  # Open the HDF file
        try:
            vs = hdf.vstart()  # Initialize the VS (Vdata) interface
            vdata = vs.attach('CoreMetadata.0')  # Access the CoreMetadata.0 Vdata by name
            metadata_records = vdata.read()  # Read the metadata fields
            metadata = ''.join(sum(metadata_records, []))  # Process the metadata
        finally:
            vdata.detach()  # Detach the Vdata
            vs.end()  # End the VS interface
            hdf.close()  # Close the HDF file
        return metadata
    def get_misrtime(self):
        # Parse the input datetime string
        _,dt,_ = self.get_datetime()
        if dt:
            formatted_str = dt.strftime("%Y.%m.%d")
            return formatted_str
        else:
            return None

    def parse_metadata_for_day_night_flag(self, metadata):
        match = re.search(r'OBJECT\s*=\s*DAYNIGHTFLAG\s*NUM_VAL\s*=\s*1\s*VALUE\s*=\s*"([^"]+)"', metadata)
        return match.group(1) if match else "DayNightFlag not found"

    def parse_metadata_for_orbit_number(self, metadata):
        match = re.search(r'OBJECT\s*=\s*ORBITNUMBER\s*CLASS\s*=\s*"1"\s*NUM_VAL\s*=\s*1\s*VALUE\s*=\s*(\d+)', metadata)
        return int(match.group(1)) if match else "Orbit Number not found"

    def get_orbit(self):
        metadata = self.get_coremetadata()
        return self.parse_metadata_for_orbit_number(metadata)

    def get_daynight(self):
        metadata = self.get_coremetadata()
        return self.parse_metadata_for_day_night_flag(metadata)
    
    
    def get_earth_sun_dist(self):
        file_ = SD(self.fullpath)
        earthsundist = getattr(file_, 'Earth-Sun Distance')
        file_.end()

        return earthsundist

    def get_data(self, fieldname, SD_field_rawData):
        hdf_file = SD(self.fullpath)
        if SD_field_rawData==0:
            data = hdf_file #science data set
        elif SD_field_rawData==1:
            data = hdf_file.select(fieldname) #field
        else:
            data = hdf_file.select(fieldname).get() #raw data
        return data, hdf_file
    
    def get_scale_and_offset(self,data_field, rad_or_ref=True):
        '''
        INPUT
            data:       numpy float array - get_data(filename,fieldname,SD_or_rawData=1)
            rad_or_ref: boolean           - True if radaince or False if reflectance offsets/scale desired
        RETURN
            2 numpy float arrays, scale factor & offset of size=number of bands
            for the chosen field
        '''
        if rad_or_ref:
            offset_name = 'radiance_offsets'
            scale_name  = 'radiance_scales'
        else:
            offset_name = 'reflectance_offsets'
            scale_name  = 'reflectance_scales'

        #grab scale and offset
     
        scale_factor = data_field.attributes()[scale_name]
        offset = data_field.attributes()[offset_name]
        return scale_factor, offset

    # def get_radiance_or_reflectance(self,data_raw, data_field, rad_or_ref = True, scale_factor_flag=True):
    #     '''
    #     INPUT
    #         data_raw:   get_data(filename, fieldname, SD_field_rawData=2)
    #         data_field: get_data(filename, fieldname, SD_field_rawData=1)
    #         rad_or_ref: boolean - True if radiance, False if reflectance
    #         scale_factor: boolean - return this as well if True
    #     RETURN
    #         radiance: numpy float array - shape=(number of bands, horizontal, vertical)
    #     '''
    #     #get dimensions of raw data
    #     num_bands = np.ma.size(data_raw, axis=0)
    #     num_horizontal = np.ma.size(data_raw, axis=1)
    #     num_vertical = np.ma.size(data_raw, axis=2)

    #     #reshape raw data to perform scale and offset correction
    
    #     data_raw_temp = np.reshape(data_raw,(num_bands, num_horizontal * num_vertical))
    #     scale_factor, offset = self.get_scale_and_offset(data_field, rad_or_ref=rad_or_ref)

    #     #fill values found in data
    #     detector_saturated = 65533
    #     detector_dead      = 65531
    #     missing_data       = 65535
    #     max_DN             = 32767
    #     min_DN             = 0

    #     #to replace DN outside of range and not any of the other fill vals
    #     fill_val_bad_data  = -999

    #     #save indices of where bad values occured occured
    #     detector_saturated_idx = np.where(data_raw_temp == detector_saturated)
    #     detector_dead_idx      = np.where(data_raw_temp == detector_dead)
    #     missing_data_idx       = np.where(data_raw_temp == missing_data)
    #     over_DN_max_idx        = np.where(data_raw_temp >  max_DN)
    #     below_min_DN_idx       = np.where(data_raw_temp <  min_DN)

    #     #mark all invalid data as nan
    #     data_raw_temp = data_raw_temp.astype(np.float)
    #     data_raw_temp[over_DN_max_idx]  = np.nan
    #     data_raw_temp[below_min_DN_idx] = np.nan
    #     #correct raw data to get radiance/reflectance values
    #     #correct first band manually
    #     if not scale_factor_flag:
    #         return data_raw_temp.reshape((num_bands, num_horizontal, num_vertical))
    #     data_corrected_total = (data_raw_temp[0,:] - offset[0]) * scale_factor[0]
    #     #for loop to put all the bands together
    #     for i in range(1,num_bands):
    #         #corrected band
    #         data_corrected = (data_raw_temp[i,:] - offset[i]) * scale_factor[i]

    #         #aggregate bands
    #         data_corrected_total = np.vstack((data_corrected_total, data_corrected))

    #     #add fill values back in
    #     data_corrected_total[over_DN_max_idx]        = fill_val_bad_data
    #     data_corrected_total[below_min_DN_idx]       = fill_val_bad_data
    #     data_corrected_total[detector_saturated_idx] = detector_saturated
    #     data_corrected_total[detector_dead_idx]      = detector_dead
    #     data_corrected_total[missing_data_idx]       = missing_data

    #     #get original shape and return radiance/reflectance
    #     return data_corrected_total.reshape((num_bands, num_horizontal, num_vertical)) 
    def destripe(self, dn, band_name):
        """
        EDF destriping via f2py-wrapped MODIS_EDF_DESTRIPE.
        Returns array in the same orientation as input: (nlines, npixels).
        """
        import numpy as np
        from .modis_edf_destripe import modis_edf_destripe  # f2py module function
        from pyhdf.SD import SD, SDC
        from pyhdf.HDF import HDF, HC

        INT32_MISSING = np.int32(-2147483648)

        # --- Config copied from your original code (zero-based ref det) ---
        DESTRIPE_CFG = {
            31: (4, [0,0,0,0,0,0,0,0,0,0]),
            32: (4, [0,0,0,0,0,0,0,0,0,0]),
            33: (2, [1,0,0,0,0,0,0,0,0,0]),
            34: (3, [0,0,0,0,0,1,1,1,0,0]),
            35: (4, [0,0,0,0,0,0,0,0,0,0]),
            36: (4, [0,0,0,0,0,0,1,0,0,0]),
        }

        def read_modis_mirror_side(hdf_path):
            h = HDF(hdf_path, HC.READ)
            vs = h.vstart()
            ref = vs.find("Level 1B Swath Metadata")
            v = vs.attach(ref)
            v.setfields("Mirror Side")
            data = v.read(v._nrecs)
            mir = np.asarray(data, dtype=np.int32).ravel()
            v.detach()
            vs.end()
            h.close()
            return mir

        def replace_bad_detectors_gumley(destriped_px_line, bad_flags, nscan):
            """
            destriped_px_line: shape (npixel, 10*nscan)
            bad_flags: length-10 list, 1=bad, 0=good
            """
            out = destriped_px_line.copy()
            bad_flags = list(bad_flags)

            for det in range(10):
                if bad_flags[det] != 1:
                    continue

                # nearest good neighbor
                rep = None
                best = 999
                for i in range(10):
                    if i != det and bad_flags[i] == 0:
                        d = abs(i - det)
                        if d < best:
                            best = d
                            rep = i
                if rep is None:
                    continue

                # replace detector line within each scan
                for s in range(nscan):
                    out[:, s*10 + det] = out[:, s*10 + rep]

            return out

        # -------------------------------
        # 1) Validate and normalize input
        # -------------------------------
        dn_arr = np.asarray(dn)

        # Keep track of missing sentinel used by your pipeline (if any)
        miss = (dn_arr.astype(np.int32, copy=False) == INT32_MISSING)

        # We must feed EDF with integer "image" equivalent to ibits(buffer,0,16)
        # If dn is already integer counts (int16/uint16/int32), do uint16->int32 conversion.
    
        if not np.issubdtype(dn_arr.dtype, np.integer):
            raise TypeError(
                f"EDF destriper expects integer scaled values (like EV_* SDS), got dtype={dn_arr.dtype}. "
                "If you currently pass radiance floats, destripe the raw EV_* integer SDS instead."
            )

        # Expect your pipeline orientation: (nlines, npixels) = (10*nscan, npixel)
        if dn_arr.ndim != 2:
            raise ValueError(f"Expected 2D array, got shape {dn_arr.shape}")

        nlines_in, npix_in = dn_arr.shape
        if nlines_in % 10 == 0:
            # shape is (nlines, npixel) -> transpose to (npixel, nlines)
            image_line_px = dn_arr
            image_px_line = image_line_px.T
            miss_px_line = miss.T
        elif npix_in % 10 == 0:
            # shape is (npixel, nlines) already
            image_px_line = dn_arr
            miss_px_line = miss
            nlines_in = image_px_line.shape[1]
            npix_in = image_px_line.shape[0]
        else:
            raise ValueError(
                f"Neither dimension is divisible by 10; cannot interpret scans. shape={dn_arr.shape}"
            )

        npixel = image_px_line.shape[0]
        nlines = image_px_line.shape[1]
        if nlines % 10 != 0:
            raise ValueError(f"Along-track dim must be multiple of 10, got {nlines}")
        nscan = nlines // 10

        # ---------------------------------------------
        # 2) Mimic Fortran driver: image = ibits(...,16)
        # ---------------------------------------------
        # Interpret raw 16 bits as unsigned then widen to int32.
        # This matches ibits(buffer,0,16) even if underlying is int16.
        image_i32 = np.asarray(image_px_line, dtype=np.uint16).astype(np.int32, copy=False)
        image_i32 = np.asfortranarray(image_i32)  # (npixel, 10*nscan), Fortran contiguous

        # ---------------------------------------
        # 3) Mirror side: enforce scan consistency
        # ---------------------------------------
        mir = read_modis_mirror_side(self.fullpath).astype(np.int32).ravel()

        # Many files store mirror side per scan (len == nscan). f2py expects nscan+1.
        # For exact matching with your subset behavior, we require mir length == nscan here.
        if mir.size != nscan and mir.size != (nscan + 1):
            raise ValueError(
                f"Mirror side length {mir.size} does not match nscan {nscan} (or nscan+1). "
                "This usually means your dn has been subset but mirror side was not subset the same way."
            )
        if mir.size == nscan:
            mir_use = np.concatenate([mir, mir[-1:]])  # length nscan+1 for f2py
        else:
            mir_use = mir  # already nscan+1
        # -----------------------------
        # 4) Call EDF destriping routine
        # -----------------------------
        ref_det, bad_flags = DESTRIPE_CFG[int(band_name)]

        destriped_i32, errflag = modis_edf_destripe(
            np.int32(ref_det),
            mir_use,
            image_i32,
            npixel=np.int32(npixel),
            nscan=np.int32(nscan),
        )
        errflag = int(errflag)

        destriped_i32 = np.asarray(destriped_i32, dtype=np.int32)

        # ---------------------------------------------------
        # 5) Replace bad detectors ONLY if errflag == 0 (ref)
        # ---------------------------------------------------
        if errflag == 0:
            destriped_i32 = replace_bad_detectors_gumley(destriped_i32, bad_flags, nscan)
        else:
            # Reference driver: if errflag != 0, it typically returns original image
            destriped_i32 = image_i32

        # -----------------------------
        # 6) Restore missing sentinel(s)
        # -----------------------------
        destriped_i32[miss_px_line] = INT32_MISSING

        if (nlines_in % 10) == 0 and dn_arr.shape[0] == nlines_in:
            return destriped_i32.T  # back to (nlines, npixel)
        else:
            return destriped_i32  # already (npixel, nlines)
    

    def get_band_data(self, band_name, data_type,scale_flag = True):
        """
        Retrieve radiance or reflectance data for a specified band.

        Parameters:
        - band_name (str): Name of the band to retrieve.
        - data_type (str): Type of data to retrieve ('radiance' or 'reflectance').

        Returns:
        - numpy.ndarray: Scaled data for the specified band.
        
        Raises:
        - ValueError: If the band name is not found or if an invalid data type is specified.
        """
        hdf_file = SD(self.fullpath, SDC.READ)        
        emissive_dataset = hdf_file.select('EV_1KM_Emissive')
        refsb_dataset = hdf_file.select('EV_1KM_RefSB')
        refsb_dataset_250m = hdf_file.select('EV_250_Aggr1km_RefSB')
        refsb_dataset_500m = hdf_file.select('EV_500_Aggr1km_RefSB')

        emissive_band_names = emissive_dataset.attributes()['band_names'].split(',')
        refsb_band_names = refsb_dataset.attributes()['band_names'].split(',')
        refsb_band_250m_names = refsb_dataset_250m.attributes()['band_names'].split(',')
        refsb_band_500m_names = refsb_dataset_500m.attributes()['band_names'].split(',')

        if band_name in emissive_band_names:
            dataset  = emissive_dataset
            band_names = emissive_band_names
            band_index = band_names.index(band_name)
            band_data_raw  = dataset[band_index, :, :]
            if self.destripe_flag and int(band_name) in [31,32,33,34,35,36]:
                band_data = self.destripe(band_data_raw,band_name)
            else:
                band_data = band_data_raw
            scales = dataset.attributes()['radiance_scales']
            offsets = dataset.attributes()['radiance_offsets']
            valid_range = dataset.attributes()['valid_range']
            band_data = band_data.astype(np.float32)
            band_data = np.where((band_data >= valid_range[0]) & (band_data <= valid_range[1]), band_data, np.nan)    
            # Scale the data using the appropriate scales and offsets
            if scale_flag:
                band_data = (band_data - offsets[band_index]) * scales[band_index]
                return band_data
            else:
                return band_data
        elif band_name in refsb_band_names:
            dataset = refsb_dataset
            band_names = refsb_band_names
            if data_type == 'radiance':
                scales = dataset.attributes()['radiance_scales']
                offsets = dataset.attributes()['radiance_offsets']
            elif data_type == 'reflectance':
                scales = dataset.attributes()['reflectance_scales']
                offsets = dataset.attributes()['reflectance_offsets']
            else:
                raise ValueError(f"Invalid data_type '{data_type}' specified. Use 'radiance' or 'reflectance'.")
        elif  band_name in  refsb_band_250m_names: 
            dataset = refsb_dataset_250m
            band_names =refsb_band_250m_names
            if data_type == 'radiance':
                scales = dataset.attributes()['radiance_scales']
                offsets = dataset.attributes()['radiance_offsets']
            elif data_type == 'reflectance':
                scales = dataset.attributes()['reflectance_scales']
                offsets = dataset.attributes()['reflectance_offsets']
            else:
                raise ValueError(f"Invalid data_type '{data_type}' specified. Use 'radiance' or 'reflectance'.")
        elif  band_name in  refsb_band_500m_names: 
            dataset = refsb_dataset_500m
            band_names =refsb_band_500m_names
            if data_type == 'radiance':
                scales = dataset.attributes()['radiance_scales']
                offsets = dataset.attributes()['radiance_offsets']
            elif data_type == 'reflectance':
                scales = dataset.attributes()['reflectance_scales']
                offsets = dataset.attributes()['reflectance_offsets']
            else:
                raise ValueError(f"Invalid data_type '{data_type}' specified. Use 'radiance' or 'reflectance'.")
 
        else:
            raise ValueError(f"Band name '{band_name}' not found in either dataset")

        # Find the index of the requested band name
        band_index = band_names.index(band_name)
        # Extract the data for the specific band
        band_data = dataset[band_index, :, :]
        # Set values outside the valid range to np.nan
        valid_range = dataset.attributes()['valid_range']
        # dn_band_data = np.copy(band_data)
        band_data = band_data.astype(np.float32)
        band_data = np.where((band_data >= valid_range[0]) & (band_data <= valid_range[1]), band_data, np.nan)    
        # Scale the data using the appropriate scales and offsets
        if scale_flag:
            band_data = (band_data - offsets[band_index]) * scales[band_index]
            return band_data
        else:
            return band_data
    
    def get_BT(self, band_name, units=1, satellite='TERRA'):
        """
        Compute brightness temperature for a MODIS infrared band.
        
        Parameters:
        rad       : 2D array of Planck radiance (units determined by UNITS)
        band      : MODIS IR band number (20-25, 27-36)
        units     : Flag defining radiance units
                    0 => milliWatts per square meter per steradian per inverse centimeter
                    1 => Watts per square meter per steradian per micron
        satellite : "TERRA" or "AQUA"
        
        Returns:
        Brightness temperature (Kelvin) as a 2D array
        """
        rad = self.get_band_data(str(band_name),'radiance')
        band = int(band_name)
        if band < 20 or band > 36 or band == 26:
            raise ValueError("Argument BAND must be in the range [20-25, 27-36]")
        
        if satellite.upper() == 'TERRA':
            # Effective central wavenumber (inverse centimeters)
            cwn = [
                2641.775, 2505.277, 2518.028, 2465.428, 2235.815, 2200.346, 0.0, 1477.967,
                1362.737, 1173.190, 1027.715, 908.0884, 831.5399, 748.3394, 730.8963, 718.8681, 704.5367
            ]
            # Temperature correction slope (no units)
            tcs = [
                0.9993411, 0.9998646, 0.9998584, 0.9998682, 0.9998819, 0.9998845, 0.0, 0.9994877,
                0.9994918, 0.9995495, 0.9997398, 0.9995608, 0.9997256, 0.9999160, 0.9999167, 0.9999191, 0.9999281
            ]
            # Temperature correction intercept (Kelvin)
            tci = [
                0.4770532, 0.09262664, 0.09757996, 0.08929242, 0.07310901, 0.07060415, 0.0, 0.2204921,
                0.2046087, 0.1599191, 0.08253401, 0.1302699, 0.07181833, 0.01972608, 0.01913568, 0.01817817, 0.01583042
            ]
        elif satellite.upper() == 'AQUA':
            # Effective central wavenumbers (inverse centimeters)
            cwn = [
                2647.409, 2511.760, 2517.908, 2462.442, 2248.296, 2209.547, 0.0, 1474.262,
                1361.626, 1169.626, 1028.740, 907.6813, 830.8411, 748.2978, 730.7766, 718.2094, 703.5007
            ]
            # Temperature correction slopes (no units)
            tcs = [
                0.9993363, 0.9998626, 0.9998627, 0.9998707, 0.9998737, 0.9998770, 0.0, 0.9995694,
                0.9994867, 0.9995270, 0.9997382, 0.9995270, 0.9997271, 0.9999173, 0.9999070, 0.9999198, 0.9999233
            ]
            # Temperature correction intercepts (Kelvin)
            tci = [
                0.4818401, 0.09426663, 0.09458604, 0.08736613, 0.07873285, 0.07550804, 0.0, 0.1848769,
                0.2064384, 0.1674982, 0.08304364, 0.1343433, 0.07135051, 0.01948513, 0.02131043, 0.01804156, 0.01683156
            ]
        else:
            raise ValueError('You must select either "TERRA" or "AQUA" as your satellite')

        def bright_m(w, r):
            """
            Compute brightness temperature given monochromatic Planck radiance.
            """
            # Constants
            h = 6.6260755e-34  # Planck constant (Joule second)
            c = 2.9979246e+8  # Speed of light in vacuum (meters per second)
            k = 1.380658e-23  # Boltzmann constant (Joules per Kelvin)

            c1 = 2.0 * h * c * c
            c2 = (h * c) / k

            # Convert wavelength to meters
            ws = 1.0E-6 * w

            # Compute brightness temperature
            return c2 / (ws * np.log(c1 / (1.0e+6 * r * ws**5) + 1.0))

        def brite_m(v, r):
            """
            Compute brightness temperature given monochromatic Planck radiance.
            """
            # Constants
            h = 6.6260755e-34  # Planck constant (Joule second)
            c = 2.9979246e+8  # Speed of light in vacuum (meters per second)
            k = 1.380658e-23  # Boltzmann constant (Joules per Kelvin)

            c1 = 2.0 * h * c * c
            c2 = (h * c) / k

            # Convert wavenumber to inverse meters
            vs = 1.0e+2 * v

            # Compute brightness temperature
            return c2 * vs / np.log(c1 * vs**3 / (1.0e-5 * r) + 1.0)

        # Compute brightness temperature
        if units == 1:
            # Radiance units are Watts per square meter per steradian per micron
            BT = (bright_m(1.0e+4 / cwn[band - 20], rad) - tci[band - 20]) / tcs[band - 20]
        else:
            # Radiance units are milliWatts per square meter per steradian per wavenumber
            BT = (brite_m(cwn[band - 20], rad) - tci[band - 20]) / tcs[band - 20]

        return BT


class MOD06Granule():
    def __init__(self, mod_l1b_fullpath):
        self.fullpath = mod_l1b_fullpath  # Use the passed full path
    @staticmethod
    def _fit_to_target(arr: np.ndarray, target_shape: tuple[int, int], pad_mode: str = "edge") -> np.ndarray:
        """Center-crop or pad (edge) to match target_shape = (rows, cols)."""
        tr, tc = target_shape
        r, c = arr.shape

        # rows
        if r > tr:
            off = (r - tr) // 2
            arr = arr[off:off + tr, :]
        elif r < tr:
            top = (tr - r) // 2
            bot = tr - r - top
            arr = np.pad(arr, ((top, bot), (0, 0)), mode=pad_mode)

        # cols
        r, c = arr.shape
        if c > tc:
            off = (c - tc) // 2
            arr = arr[:, off:off + tc]
        elif c < tc:
            left = (tc - c) // 2
            right = tc - c - left
            arr = np.pad(arr, ((0, 0), (left, right)), mode=pad_mode)

        return arr
    def get_data(
        self,
        fieldname: str,
        upscale_1km: bool = False,
        scale_factor_flag: bool = True,
        *,
        target_shape: tuple[int, int] | None = None,
        target_like: np.ndarray | None = None,
        pad_mode: str = "edge",
    ):
        """
        Read a MOD06 SDS, optionally upscale 5 km → 1 km using the SDS's
        Cell_*_Swath_Sampling, crop to its exact 1 km window, and (optionally)
        fit to a caller-provided target grid.
        """
        hdf_file = SD(self.fullpath)
        sds = hdf_file.select(fieldname)
        attrs = sds.attributes()

        # ---- read ----
        data_value = sds.get()

        # ---- scale / mask (if requested) ----
        if scale_factor_flag:
            sf = float(attrs.get('scale_factor', 1.0))
            off = float(attrs.get('add_offset', 0.0))
            vr = attrs.get('valid_range', None)
            fill = attrs.get('_FillValue', None)

            data_value = data_value.astype(np.float32, copy=False)
            if fill is not None:
                data_value = np.where(data_value == fill, np.nan, data_value)
            if vr is not None:
                lo, hi = float(vr[0]), float(vr[1])
                data_value = np.where((data_value >= lo) & (data_value <= hi), data_value, np.nan)
            data_value = (data_value - off) * sf
        else:
            # ensure float so NaNs (from padding/cropping) are valid
            if np.issubdtype(data_value.dtype, np.integer):
                data_value = data_value.astype(np.float32, copy=False)

        # ---- optional upscale using SDS sampling triplets ----
        # triplets are (start, length, step) on the 1 km grid
        sA = attrs.get('Cell_Along_Swath_Sampling', (0, data_value.shape[0], 1))
        sC = attrs.get('Cell_Across_Swath_Sampling', (0, data_value.shape[1], 1))
        startA, lenA, incA = [int(x) for x in sA]
        startC, lenC, incC = [int(x) for x in sC]

        if upscale_1km and (incA != 1 or incC != 1):
            # repeat 5 km → 1 km
            rep = np.repeat(np.repeat(data_value, incA, axis=0), incC, axis=1)

            # guard occasional off-by-one in starts provided by some granules
            if startA + lenA > rep.shape[0] and startA > 0:
                startA -= 1
            if startC + lenC > rep.shape[1] and startC > 0:
                startC -= 1

            # crop the exact 1 km window advertised by the SDS
            data_value = rep[startA:startA + lenA, startC:startC + lenC]

        # ---- fit to a target grid if provided ----
        if target_like is not None and target_shape is None:
            target_shape = target_like.shape
        if target_shape is not None:
            data_value = self._fit_to_target(data_value, target_shape, pad_mode=pad_mode)

        return data_value
        
    def get_ctp(self):
        return self.get_data('cloud_top_pressure_1km')    

    def get_cttemp(self):
        return self.get_data('cloud_top_temperature_1km') 
    
    def get_opt(self):
        return self.get_data('Cloud_Optical_Thickness')
   
    def get_cer(self):
        return self.get_data('Cloud_Effective_Radius') 
    
    def get_ctm(self):
        return self.get_data('cloud_top_method_1km') 

    def get_multi(self):
        return self.get_data('Cloud_Multi_Layer_Flag')
    
    def get_qa(self):
        """
        Read MOD06 Quality_Assurance_1km as raw packed bytes.

        Returns
        -------
        np.ndarray
            Shape (ny, nx, 9), dtype int8.
        """
        hdf_file = SD(self.fullpath, SDC.READ)
        try:
            sds = hdf_file.select("Quality_Assurance_1km")
            qa = sds.get()  # raw packed bytes

            if np.ma.isMaskedArray(qa):
                qa = qa.filled(0)

            qa = np.asarray(qa)

            # Keep raw byte values; do NOT scale or apply valid_range masking
            if qa.dtype != np.int8:
                qa = qa.astype(np.int8, copy=False)

            return qa
        finally:
            hdf_file.end()
                  
    def get_cphase(self):
        '''
        0: Cloud Mask unavaible, missing data
        1: Cloud Mask Clear or Probably Clear, or Pixel Restored to clear sky
        2: Liquid Water
        3: Ice
        4: Undertermined
        '''
        return self.get_data('Cloud_Phase_Optical_Properties')    
    
    def get_multilayer_withPH(self):
        '''
        0 indicates no cloud was detected, 
        1 indicates a monolayer cloud and values
        2–10 indicate the cumulative weight of the positive multilayer tests. 
        '''

        return self.get_data('Cloud_Multi_Layer_Flag')
    
    def get_multilayer_withoutPH(self, sum_flag = False):

        '''
        Use quality-assurance at 1km  to turn off the Pavolonis–Heidinger (PH) test
        0 indicates no cloud was detected, 
        1 indicates a monolayer cloud and values
        2 indicates multi-layer clouds without the PH test

        '''

        cloud_multi_layer_flag = self.get_data('Cloud_Multi_Layer_Flag')
        quality_assurance_1km  = self.get_data('Quality_Assurance_1km')
        byte5 = quality_assurance_1km[:, :, 5]


        # Create a bitmask to isolate Bits 0-3 (other multilayer tests)
         
        other_tests_mask = 0b00001111

        other_tests_result = byte5 & other_tests_mask
 
        
        multilayer_detected = other_tests_result != 0

        new_cloud_multi_layer_flag = np.copy(cloud_multi_layer_flag)

        # Set pixels with multilayer detection (excluding PH test) to 2 (or another appropriate value)
        if  sum_flag:
            positive_tests_count = (
                (other_tests_result & 0b0001) +
                ((other_tests_result & 0b0010) >> 1) +
                ((other_tests_result & 0b0100) >> 2) +
                ((other_tests_result & 0b1000) >> 3)
            )
            new_cloud_multi_layer_flag = positive_tests_count  # Multilayer cloud detected
        else:
            new_cloud_multi_layer_flag[~multilayer_detected] = 2  # Single-layer cloud
        # Set pixels without multilayer detection to 1 (single-layer cloud)
        new_cloud_multi_layer_flag[~multilayer_detected] = 1  # Single-layer cloud
        # Handle pixels where no cloud is detected using Cloud_Multi_Layer_Flag
        no_cloud = cloud_multi_layer_flag == 0
        new_cloud_multi_layer_flag[no_cloud] = 0  # No cloud detected

        return new_cloud_multi_layer_flag

    def get_cth(self):
        return self.get_data('cloud_top_height_1km') 
    
    def get_emissivity(self,upscale_1km = True):
        # cth = self.get_cth()
        # emiss_1km = self.get_data('Cloud_Effective_Emissivity',upscale_1km = upscale_1km,target_like=cth) 
        return self.get_data('cloud_emissivity_1km') 

class MOD03Granule:
    def __init__(self, mod_l1b_fullpath):
        self.fullpath = mod_l1b_fullpath

        # simple caches so repeated calls do not reread the same SDSs
        self._cache = {}

    def get_data(self, fieldname, scale_factor_flag=True):
        hdf_file = SD(self.fullpath, SDC.READ)
        try:
            data = hdf_file.select(fieldname)
            data_value = data.get()

            if scale_factor_flag and "scale_factor" in data.attributes():
                scale_factor = data.attributes()["scale_factor"]
                data_value = data_value * scale_factor

            return data_value
        finally:
            hdf_file.end()

    def list_sds_names(self):
        hdf_file = SD(self.fullpath, SDC.READ)
        try:
            return list(hdf_file.datasets().keys())
        finally:
            hdf_file.end()

    def get_sza(self):
        return self.get_data("SolarZenith")

    def get_saa(self):
        return self.get_data("SolarAzimuth")

    def get_vza(self):
        return self.get_data("SensorZenith")

    def get_vaa(self):
        return self.get_data("SensorAzimuth")

    def get_lat(self):
        return self.get_data("Latitude", scale_factor_flag=False)

    def get_lon(self):
        return self.get_data("Longitude", scale_factor_flag=False)

    def get_landsea_mask(self, binary_output=True):
        """
        0: Shallow ocean
        1: Land
        2: Coastline
        3: Shallow inland water
        4: Ephemeral water
        5: Deep inland water
        6: Moderate or continental ocean
        7: Deep ocean
        """
        landsea_mask = self.get_data("Land/SeaMask", scale_factor_flag=False)

        if binary_output:
            binary_mask = np.zeros_like(landsea_mask, dtype=int)
            binary_mask[landsea_mask == 1] = 1
            binary_mask[landsea_mask == 2] = 2
            return binary_mask

        return landsea_mask

    # -----------------------------
    # Cached MOD03 timing SDS readers
    # -----------------------------
    def get_scan_number(self):
        if "Scan number" not in self._cache:
            self._cache["Scan number"] = self.get_data(
                "Scan number", scale_factor_flag=False
            ).astype(np.int32)
        return self._cache["Scan number"]

    def get_ev_frames(self):
        if "EV frames" not in self._cache:
            self._cache["EV frames"] = self.get_data(
                "EV frames", scale_factor_flag=False
            ).astype(np.int32)
        return self._cache["EV frames"]

    def get_ev_start_time(self):
        if "EV start time" not in self._cache:
            arr = self.get_data("EV start time", scale_factor_flag=False).astype(np.float64)
            arr[arr <= -1.0e9] = np.nan
            self._cache["EV start time"] = arr
        return self._cache["EV start time"]

    def get_ev_center_time(self):
        if "EV center time" not in self._cache:
            arr = self.get_data("EV center time", scale_factor_flag=False).astype(np.float64)
            arr[arr <= -1.0e9] = np.nan
            self._cache["EV center time"] = arr
        return self._cache["EV center time"]

    def get_mirror_side(self):
        if "Mirror side" not in self._cache:
            self._cache["Mirror side"] = self.get_data(
                "Mirror side", scale_factor_flag=False
            ).astype(np.int16)
        return self._cache["Mirror side"]

    # -----------------------------
    # Time conversion helpers
    # -----------------------------
    @staticmethod
    def _tai_minus_utc_for_utc_date(dt_utc):
        """
        Return TAI-UTC (seconds) for a UTC datetime.
        """
        leap_table = [
            (datetime(1993, 1, 1, tzinfo=timezone.utc), 27),
            (datetime(1993, 7, 1, tzinfo=timezone.utc), 28),
            (datetime(1994, 7, 1, tzinfo=timezone.utc), 29),
            (datetime(1996, 1, 1, tzinfo=timezone.utc), 30),
            (datetime(1997, 7, 1, tzinfo=timezone.utc), 31),
            (datetime(1999, 1, 1, tzinfo=timezone.utc), 32),
            (datetime(2006, 1, 1, tzinfo=timezone.utc), 33),
            (datetime(2009, 1, 1, tzinfo=timezone.utc), 34),
            (datetime(2012, 7, 1, tzinfo=timezone.utc), 35),
            (datetime(2015, 7, 1, tzinfo=timezone.utc), 36),
            (datetime(2017, 1, 1, tzinfo=timezone.utc), 37),
        ]

        tai_minus_utc = leap_table[0][1]
        for eff_date, val in leap_table:
            if dt_utc >= eff_date:
                tai_minus_utc = val
            else:
                break
        return tai_minus_utc

    @staticmethod
    def _native_mod03_to_seconds_since_2000(native_seconds):
        """
        Convert native MOD03 EV times to UTC-based
        seconds since 2000-01-01T00:00:00Z.

        Assumes MOD03 EV times are TAI93-like seconds:
        seconds since 1993-01-01 on the MODIS TAI-related timescale.
        """
        arr = np.asarray(native_seconds, dtype=np.float64)
        out = np.full(arr.shape, np.nan, dtype=np.float64)

        epoch_1993_utc = datetime(1993, 1, 1, tzinfo=timezone.utc)
        epoch_2000_utc = datetime(2000, 1, 1, tzinfo=timezone.utc)

        # At 1993-01-01 UTC, TAI-UTC = 27 s
        tai_minus_utc_at_1993_epoch = 27.0

        flat_in = arr.ravel()
        flat_out = out.ravel()

        for i, x in enumerate(flat_in):
            if not np.isfinite(x):
                continue

            # first guess using epoch offset only
            utc_guess = epoch_1993_utc.timestamp() + (x - tai_minus_utc_at_1993_epoch)
            dt_guess = datetime.fromtimestamp(utc_guess, tz=timezone.utc)

            tai_minus_utc_now = MOD03Granule._tai_minus_utc_for_utc_date(dt_guess)

            # convert native TAI93-like seconds to UTC seconds
            utc_seconds = epoch_1993_utc.timestamp() + (
                x - (tai_minus_utc_now - tai_minus_utc_at_1993_epoch)
            )

            flat_out[i] = utc_seconds - epoch_2000_utc.timestamp()

        return out

    # -----------------------------
    # Time array builders
    # -----------------------------
    def get_time_2d_native(self, use_center_time=True, fps=3000.0):
        """
        Original loop-based version.
        Returns shape (nscans*10, mframes) in native MOD03 time units.
        """
        ev_frames = self.get_ev_frames()
        nscans = ev_frames.size
        nrows = nscans * 10
        max_frames = int(ev_frames.max())

        if use_center_time:
            t_scan = self.get_ev_center_time()
        else:
            t_scan = self.get_ev_start_time()

        time2d = np.full((nrows, max_frames), np.nan, dtype=np.float64)

        for iscan in range(nscans):
            nframes = int(ev_frames[iscan])

            if nframes <= 0 or np.isnan(t_scan[iscan]):
                continue

            frame_idx = np.arange(nframes, dtype=np.float64)

            if use_center_time:
                frame_times = t_scan[iscan] + (frame_idx - (nframes - 1) / 2.0) / fps
            else:
                frame_times = t_scan[iscan] + (frame_idx + 0.5) / fps

            row0 = iscan * 10
            row1 = row0 + 10
            time2d[row0:row1, :nframes] = frame_times[None, :]

        return time2d

    def get_time_2d_since_2000_fast(self, use_center_time=True, fps=3000.0):
        """
        Pixel observation time was derived from the matched MOD03 geolocation granule. MOD03 provides the Earth-view center time for each scan, 
        and the MODIS scan geometry provides the number of Earth-view frames within each scan together with the frame sampling rate.
        Using the scan center time as the reference, a time offset was calculated for each cross-track frame,
        yielding a time for each frame within the scan. 
        That frame time was then assigned to the corresponding pixels and replicated across the 10 along-track detector rows associated with that scan as a scan-level approximation. 
        The resulting time field was then converted to UTC-based seconds since 2000-01-01 00:00:00 for storage in the MM product.
        """
        ev_frames = self.get_ev_frames().astype(np.int32)
        nscans = ev_frames.size
        nrows = nscans * 10
        max_frames = int(ev_frames.max())

        if use_center_time:
            t_scan_native = self.get_ev_center_time()
        else:
            t_scan_native = self.get_ev_start_time()

        t_scan_2000 = self._native_mod03_to_seconds_since_2000(t_scan_native)

        # frame index array: shape (1, max_frames)
        frame_idx = np.arange(max_frames, dtype=np.float64)[None, :]

        # number of valid frames per scan: shape (nscans, 1)
        nframes2d = ev_frames[:, None].astype(np.float64)

        if use_center_time:
            offsets = (frame_idx - (nframes2d - 1) / 2.0) / fps
        else:
            offsets = (frame_idx + 0.5) / fps
        scan_frame_time = t_scan_2000[:, None] + offsets

        # mask columns beyond EV frames for each scan
        valid = frame_idx < nframes2d
        scan_frame_time = np.where(valid, scan_frame_time, np.nan)


        time2d = np.repeat(scan_frame_time, 10, axis=0)

        if time2d.shape != (nrows, max_frames):
            raise ValueError(
                f"Unexpected time2d shape {time2d.shape}, expected {(nrows, max_frames)}"
            )

        return time2d

    def get_time_2d_since_2000(self, use_center_time=True, fps=3000.0, fast=True):
        """
        Main public 2D method.
        """
        if fast:
            return self.get_time_2d_since_2000_fast(
                use_center_time=use_center_time, fps=fps
            )

        time2d_native = self.get_time_2d_native(
            use_center_time=use_center_time, fps=fps
        )
        return self._native_mod03_to_seconds_since_2000(time2d_native)

    def get_time_1d_since_2000(self, use_center_time=True):
        """
        Fast 1D along-track time array of length nscans*10.

        One representative scan time is assigned to each set of 10 rows.
        """
        if use_center_time:
            t_scan_native = self.get_ev_center_time()
        else:
            t_scan_native = self.get_ev_start_time()

        t_scan_2000 = self._native_mod03_to_seconds_since_2000(t_scan_native)
        return np.repeat(t_scan_2000, 10) 