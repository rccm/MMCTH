PROGRAM HDF_DESTRIPE_NEW

!F90
!
!DESCRIPTION:
!    Destripe a MODIS L1B 1KM HDF file.
!
!INPUT PARAMETERS:
!    MODIS L1B 1KM file (via Process Control File)
!    MODIS destriping configuration file (via Process Control File)
!
!OUTPUT PARAMETERS:
!    MODIS L1B 1KM file with destriped scaled integers
!    (NOTE: The input L1B 1KM file is ireversibly modified)
!
!REVISION HISTORY:
!    Created by Liam.Gumley@ssec.wisc.edu
!
!TEAM-UNIQUE HEADER:
!    Copyright(C) 2004, University of Wisconsin-Madison MODIS Group
!
!END

IMPLICIT NONE

include 'hdf.f90'
include 'MOD_PRDS.inc'

!- Parameters
integer :: lines_1km
parameter (lines_1km = 10)

!- Local variables
character :: rcsid*100
integer :: band_data(12, 36)
character :: header*132
character :: in_file*256
integer :: hdfid, varid, attid, rtn
integer :: rank, dimsizes(3), data_type, num_attrs
character :: sds_name*100
integer :: num_pixels, num_lines, num_scans
integer :: band_loop, band
integer :: start(3), stride(3), edge(3)
integer :: ref_det
integer :: errflag
character :: pgename*8, errtext*512
integer :: pcf_ver
integer :: det_index, min_diff, i, det_diff, rep_index

!- Allocatable arrays
integer, dimension(:), allocatable :: mirror_side
integer, dimension(:, :), allocatable :: image, destripe
integer(selected_int_kind(4)), dimension(:, :), allocatable :: buffer

!- External functions
integer :: sfstart, sfselect, sfn2index, sfginfo, sfendacc, sfend, &
           sfrdata, sfwdata, sffattr, sfscatt
external   sfstart, sfselect, sfn2index, sfginfo, sfendacc, sfend, &
           sfrdata, sfwdata, sffattr, sfscatt
integer :: get_file_name, get_modis_mirror, get_band_config, get_band_index, strlen
external   get_file_name, get_modis_mirror, get_band_config, get_band_index, strlen

integer :: num_args
integer :: FlagRA
character :: FlagBuff*10

!- Data statements
data rcsid /'$Id: hdf_destripe_new.f90,v 1.8 2004/06/24 14:25:27 gumley Exp $'/
data pgename /'MOD_PRDS'/

!------------------------------------------------------------------------------
! INITIALIZATION
!------------------------------------------------------------------------------

num_args = command_argument_count()

if(num_args == 1) then
   call get_command_argument(1,FlagBuff)
   read (FlagBuff,*) FlagRA
else
   !This is the default value
   FlagRA = 0
endif

!- Print progress message to logfile
write(errtext, '(''Started MODIS L1B 1KM destriping version: '', a)') rcsid
call message(pgename, errtext, 0, 3)

if( FlagRA == 1) then
   !- Get the MODIS L1B 1KM input file name
   !- (MOD_PCF_NUM is defined in MOD_PRDS.inc)
   pcf_ver = 1
   rtn = get_file_name(MOD_PCF_NUM_RA, pcf_ver, in_file)
   
   if (rtn /= 0) then
      errtext = 'Error getting MODIS L1B 1KM filename from PCF [OPERATOR ACTION: Check PCF]'
      call message('get_file_name', errtext, 0, 2)
   endif
else
    !- Get the MODIS L1B 1KM input file name
   !- (MOD_PCF_NUM is defined in MOD_PRDS.inc)
   pcf_ver = 1
   rtn = get_file_name(MOD_PCF_NUM, pcf_ver, in_file)
   
   if (rtn /= 0) then
      errtext = 'Error getting MODIS L1B 1KM filename from PCF [OPERATOR ACTION: Check PCF]'
      call message('get_file_name', errtext, 0, 2)
   endif

endif

!- Open the input file in SDS mode for read/write
!- (DFACC_RDWR and FAIL are defined in hdf.f90)
hdfid = sfstart(in_file, DFACC_RDWR)
if (hdfid == FAIL) then
   write(errtext, '(''Error opening MODIS L1B 1KM file: '', a)') in_file
  call message(pgename, errtext, 0, 2)
endif

!- Select the 1KM emissive band SDS
varid = sfselect(hdfid, sfn2index(hdfid, 'EV_1KM_Emissive'))
if (varid == FAIL) then
  errtext = 'MODIS L1B 1KM emissive bands were not found in input file [OPERATOR ACTION: Contact SDST]'
  call message(pgename, errtext, 0, 2)
endif

!- Check if the file is already destriped
attid = sffattr(hdfid, 'UW_DESTRIPE_LWIR')
if (attid /= FAIL) then
  errtext = 'MODIS L1B 1KM input file is already destriped [OPERATOR ACTION: Contact SDST]'
  call message(pgename, errtext, 0, 2)
endif

!- Print progress message to logfile
write(errtext, '(''Opened MODIS L1B 1KM input file: '', a)') in_file
call message(pgename, errtext, 0, 3)

!------------------------------------------------------------------------------
! GET INPUT FILE INFORMATION
!------------------------------------------------------------------------------

!- Get the number of pixels, lines, and scans
rtn = sfginfo(varid, sds_name, rank, dimsizes, data_type, num_attrs)
rtn = sfendacc(varid)
num_pixels = dimsizes(1)
num_lines  = dimsizes(2)
num_scans  = num_lines / lines_1km
if (num_scans <= 1) then
  errtext = 'Number of MODIS L1B 1KM scans is <= 1 [OPERATOR ACTION: Contact SDST]'
  call message(pgename, errtext, 0, 2)
endif

!- Allocate runtime arrays
allocate(image(num_pixels, num_lines))
allocate(destripe(num_pixels, num_lines))
allocate(buffer(num_pixels, num_lines))
allocate(mirror_side(num_scans))

!- Get mirror side data
rtn = get_modis_mirror(in_file, mirror_side)
if (rtn /= 0) then
  errtext = 'Error reading MODIS L1B 1KM mirror side data [OPERATOR ACTION: Contact SDST]'
  call message(pgename, errtext, 0, 2)
endif

!- Set start, stride, edge arrays
start(1:3) = 0
stride(1:3) = 1
edge(1) = num_pixels
edge(2) = num_lines
edge(3) = 1

!- Get MODIS destriping band configuration
!- (CFG_PCF_NUM is defined in MOD_PRDS.inc)
pcf_ver = 1
rtn = get_band_config(CFG_PCF_NUM, pcf_ver, band_data, header)
if (rtn /= 0) then
  errtext = 'Error reading MODIS destriping configuration file [OPERATOR ACTION: Contact SDST]'
  call message(pgename, errtext, 0, 2)
endif

!- Report destriping configuration to logfile
call message(pgename, header, 0, 3)

!------------------------------------------------------------------------------
! DESTRIPE EACH BAND
!------------------------------------------------------------------------------

!- Loop over each band to be destriped
do band_loop  = 1, 36

  !- Get band number
  band = band_data(1, band_loop)
  if (band == 0) goto 20

  !- Start access to SDS for this band
  if (band >=  1 .and. band <=  2) &
    varid = sfselect(hdfid, sfn2index(hdfid, 'EV_250_Aggr1km_RefSB'))
  if (band >=  3 .and. band <=  7) &
    varid = sfselect(hdfid, sfn2index(hdfid, 'EV_500_Aggr1km_RefSB'))
  if (band >=  8 .and. band <= 19) &
    varid = sfselect(hdfid, sfn2index(hdfid, 'EV_1KM_RefSB'))
  if (band >= 20 .and. band <= 25) &
    varid = sfselect(hdfid, sfn2index(hdfid, 'EV_1KM_Emissive'))
  if (band == 26) &
    varid = sfselect(hdfid, sfn2index(hdfid, 'EV_1KM_RefSB'))
  if (band >= 27 .and. band <= 36) &
    varid = sfselect(hdfid, sfn2index(hdfid, 'EV_1KM_Emissive'))

  !- Get band index in SDS array for this band
  start(3) = get_band_index(band)

  !- Read the input image for this band
  rtn = sfrdata(varid, start, stride, edge, buffer)

  !- Convert input image from short to long
  image = ibits(buffer, 0, 16)

  !- Get reference detector for this band
  ref_det = band_data(2, band_loop)

  !- Compute destriped image
  call modis_edf_destripe(num_pixels, num_scans, ref_det, mirror_side, &
    image, destripe, errflag)

  !- If destriping was successful, write the destriped image
  if (errflag == 0) then

    !- Replace bad detectors with nearest good neighbor
    do det_index = 1, 10

      !- Check if this detector should be replaced
      if (band_data(det_index + 2, band_loop) == 1) then

        !- Get the index of the nearest good detector
        min_diff = 9
        do i = 1, 10
          det_diff = abs(i - det_index)
          if (i /= det_index .and. &
              det_diff < min_diff .and. &
              band_data(i + 2, band_loop) == 0) then
            rep_index = i
            min_diff = det_diff
          endif
        end do

        !- Replace the data for this detector
        do i = 1, num_scans
          destripe(:, (i - 1) * 10 + det_index) = destripe(:, (i - 1) * 10 + rep_index)
        end do

      endif

    end do

    ! Convert destriped image from long to short
    buffer = ibits(destripe, 0, 16)

    !- Write the destriped image for this band
    rtn = sfwdata(varid, start, stride, edge, buffer)

    !- Print success message to logfile
    write(errtext, '(''Successfully destriped band '', 1x, i2)') band
    call message(pgename, errtext, 0, 3)

  else

    !- Print failure message to logfile
    write(errtext, '(''Could not destripe band '', 1x, i2)') band
    call message(pgename, errtext, 0, 1)

  endif

  !- End access to SDS for this band
  rtn = sfendacc(varid)

  20 continue

end do

!------------------------------------------------------------------------------
! CLEANUP AND EXIT
!------------------------------------------------------------------------------

!- Write a new global attribute to show this file is destriped
rtn = sfscatt(hdfid, 'UW_DESTRIPE_LWIR', DFNT_CHAR8, strlen(rcsid), rcsid)

!- Write a new global attribute to record destriping configuration
rtn = sfscatt(hdfid, 'UW_DESTRIPE_CONFIG', DFNT_CHAR8, strlen(header), header)

!- Close the input file
rtn = sfend(hdfid)

!- Deallocate the input and destriped image arrays
deallocate(image)
deallocate(destripe)
deallocate(buffer)
deallocate(mirror_side)

!- Print progress message to logfile
write(errtext, '(''Successfully destriped MODIS L1B 1KM file: '', a)') in_file
call message(pgename, errtext, 0, 3)

!- Return exit code zero to shell
call exit(0)

END PROGRAM HDF_DESTRIPE_NEW

!------------------------------------------------------------------------------

FUNCTION GET_FILE_NAME(PCF_NUM, PCF_VER, PCF_NAME)

!F90
!
!DESCRIPTION:
!    Get the name of a file defined in the Process Control File (PCF).
!
!INPUT PARAMETERS:
!    PCF_NUM    Logical number of the file in the PCF
!    PCF_VER    Version number of the file in the PCF
!
!OUTPUT PARAMETERS:
!    PCF_NAME   Name of the file
!
!REVISION HISTORY:
!    Created by Liam.Gumley@ssec.wisc.edu
!
!TEAM-UNIQUE HEADER:
!    Copyright(C) 2004, University of Wisconsin-Madison MODIS Group
!
!END

IMPLICIT NONE

include 'PGS_SMF.f'

!- Arguments
integer, intent(in) :: pcf_num, pcf_ver
character(*), intent(out) :: pcf_name
integer :: get_file_name

!- Local variables
integer :: result

!- External functions
integer :: pgs_pc_getreference
external   pgs_pc_getreference

!- Call PGS toolkit routine to get file name given PCF number
!- (PGS_S_SUCCESS is defined in PGS_SMF.f)
result = pgs_pc_getreference(pcf_num, pcf_ver, pcf_name)
if (result /= PGS_S_SUCCESS) then
  get_file_name = -1
else
  get_file_name = 0
endif

END FUNCTION GET_FILE_NAME

!------------------------------------------------------------------------------

FUNCTION GET_MODIS_MIRROR(IN_FILE, MIRROR_SIDE)

!F90
!
!DESCRIPTION:
!    Get scan mirror data from a MODIS 1KM radiance HDF file.
!
!INPUT PARAMETERS:
!    IN_FILE    Name of the MODIS 1KM radiance HDF file
!
!OUTPUT PARAMETERS:
!    MIRROR_SIDE    Array of mirror side data (values are 0 or 1)
!
!REVISION HISTORY:
!    Created by Liam.Gumley@ssec.wisc.edu
!
!TEAM-UNIQUE HEADER:
!    Copyright(C) 2004, University of Wisconsin-Madison MODIS Group
!
!END

IMPLICIT NONE

include 'hdf.f90'

!- Arguments
character(*), intent(in) :: in_file
integer, intent(out) :: mirror_side(*)
integer :: get_modis_mirror

!- Local variables
integer :: file_id, num_dds_block, rtn, vdata_ref, vdata_id, &
  num_scans, record_index

!- External functions
integer :: hopen, vfstart, vsffnd, vsfatch, vsfsfld, vsfelts, &
  vsfseek, vsfrd, vsfrdc, vsfdtch, vfend, hclose
external   hopen, vfstart, vsffnd, vsfatch, vsfsfld, vsfelts, &
  vsfseek, vsfrd, vsfrdc, vsfdtch, vfend, hclose

!- Open the input HDF file for read only (DFACC_READ is defined in hdf.inc)
num_dds_block = 0
file_id = hopen(in_file, DFACC_READ, num_dds_block)
if (file_id == -1) then
  get_modis_mirror = -1
  return
endif

!- Start vdata interface
rtn = vfstart(file_id)

!- Get reference number of vdata containing scan type and mirror side
vdata_ref = vsffnd(file_id, 'Level 1B Swath Metadata')
if (vdata_ref == 0) then
  get_modis_mirror = -2
  return
endif

!- Attach to vdata
vdata_id = vsfatch(file_id, vdata_ref, 'r')

!- Get number of records (scans)
num_scans = vsfelts(vdata_id)

!- Read mirror side data
record_index = 0
rtn = vsfseek(vdata_id, record_index)
rtn = vsfsfld(vdata_id, 'Mirror Side')
rtn = vsfrd(vdata_id, mirror_side, num_scans, FULL_INTERLACE)

!- Detach from vdata
rtn = vsfdtch(vdata_id)

!- End vdata interface
rtn = vfend(file_id)

!- Close the input file
rtn = hclose(file_id)

!- Set return flag
get_modis_mirror = 0

END FUNCTION GET_MODIS_MIRROR

!------------------------------------------------------------------------------

FUNCTION GET_BAND_INDEX(BAND)

!F90
!
!DESCRIPTION:
!    Get index within a SDS for a given band in a MODIS 1KM radiance HDF file
!
!    EV_250_Aggr1km_RefSB contains bands 1,2
!    EV_500_Aggr1km_RefSB contains bands 3,4,5,6,7
!    EV_1KM_RefSB contains bands 8,9,10,11,12,13lo,13hi,14lo,14hi,15,16,17,18,19,26
!    EV_1KM_Emissive contains bands 20,21,22,23,24,25,27,28,29,30,31,32,33,34,35,36
!
!INPUT PARAMETERS:
!    BAND    MODIS band number
!
!OUTPUT PARAMETERS:
!    GET_BAND_INDEX    Index of the requested MODIS band number within the SDS
!
!REVISION HISTORY:
!    Created by Liam.Gumley@ssec.wisc.edu
!
!TEAM-UNIQUE HEADER:
!    Copyright(C) 2004, University of Wisconsin-Madison MODIS Group
!
!END

IMPLICIT NONE

!- Arguments
integer, intent(in) :: band
integer :: get_band_index

!- Local variables
integer :: band_index(36)

!- Data statements
data band_index / &
  0, 1, &
  0, 1, 2, 3, 4, &
  0, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13, &
  0, 1, 2, 3, 4, 5, 14, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 /

!- Get the band index with the SDS for this band
get_band_index = band_index(band)

END FUNCTION GET_BAND_INDEX

!------------------------------------------------------------------------------

FUNCTION GET_BAND_CONFIG(PCF_NUM, PCF_VER, BAND_DATA, HEADER)

!F90
!
!DESCRIPTION:
!    Read MODIS destriping band configuration file
!
!INPUT PARAMETERS:
!    PCF_NUM    Logical number of the file in the PCF
!    PCF_VER    Version number of the file in the PCF
!
!OUTPUT PARAMETERS:
!    BAND_DATA    Array of destriping configuration information
!    HEADER       Text string containing configuration file header
!
!REVISION HISTORY:
!    Created by Liam.Gumley@ssec.wisc.edu
!
!TEAM-UNIQUE HEADER:
!    Copyright(C) 2004, University of Wisconsin-Madison MODIS Group
!
!END

IMPLICIT NONE

include 'PGS_IO.f'
include 'PGS_SMF.f'

!- Arguments
integer, intent(in)  :: pcf_num, pcf_ver
integer, intent(out) :: band_data(12, 36)
character(*), intent(out) :: header
character errtext*200
integer :: get_band_config

!- Local variables
integer :: record_length, lun, rtn, iostat, row, i
character :: string*100

!- External functions
integer :: pgs_io_gen_openf, pgs_io_gen_closef
external   pgs_io_gen_openf, pgs_io_gen_closef

!- Open destriping configuration file successfully or exit fatal error
!- (PGSD_IO_GEN_RSEQFRM is defined in PGS_SMF.f)
!- (PGS_S_SUCCESS is defined in PGS_IO.f)

record_length = 1
rtn = pgs_io_gen_openf(pcf_num, PGSD_IO_GEN_RSEQFRM, record_length, lun, pcf_ver)

if (rtn /= PGS_S_SUCCESS) then
! rhucek 11/22/05: added custom error message and exit fatal error
  errtext = 'Error opening MODIS destriping configuration file: '   &
            // '[OPERATOR ACTION: Contact SDST]'
  call message('GET_BAND_CONFIG', errtext, 0, 2)
endif

!- Set the band_data array to zeroes
iostat    = 0
band_data = 0

!- Read the header line successfully or exit fatal error
read (lun, '(a)', iostat=iostat) header

! rhucek 11/22/05: added error check and if error, exit immediately
if (iostat /= 0) then
   errtext = 'Error reading MODIS destriping configuration file header: '  &
             // '[OPERATOR ACTION: Contact SDST]'
   call message('GET_BAND_CONFIG', errtext, 0, 2)
endif

!- Read all remaining lines from the input file
iostat = 0
do row = 1, 36

  !- Read the next line
  read (lun, '(a)', iostat=iostat) string
  if (iostat /= 0)  goto 10

  !- Decode the band data for this line
  read (string, *, iostat=iostat) (band_data(i, row), i = 1, 12)
  if (iostat /= 0) then
    get_band_config = -2
    return
  endif

end do

10 continue

! rhucek 11/22/05: distinguish read errors (iostat>0) from end of file (iostat<0)
if (iostat < 0) then               ! end of file
   rtn = pgs_io_gen_closef(lun)    ! Close input file

   if (rtn /= PGS_S_SUCCESS) then
      get_band_config = -2
      return
   endif
else                               ! read error
   errtext = 'Error reading MODIS destriping configuration file header: ' &
             // '[OPERATOR ACTION: Contact SDST]'
   call message('GET_BAND_CONFIG', errtext, 0, 2)
endif

!- Set return flag
get_band_config = 0

END FUNCTION GET_BAND_CONFIG
