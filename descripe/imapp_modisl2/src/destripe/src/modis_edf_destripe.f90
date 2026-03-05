!-------------------------------------------------------------------------------

SUBROUTINE MODIS_EDF_DESTRIPE(NPIXEL, NSCAN, REF, MIR, IMAGE, DESTRIPE, ERRFLAG)

!F90
!
!DESCRIPTION:
!    Destripe one band of MODIS 1KM scaled integer data using the EDF
!    algorithm of Weinreb et. al.
!
!INPUT PARAMETERS:
!    NPIXEL      Number of 1KM pixels across scan (usually 1354)
!    NSCAN       Number of earth scans
!    REF         Reference detector for this band (zero-based; product order)
!    MIR         Array of mirror size values for each earth scan
!    IMAGE       Array of MODIS 1KM scaled integers (converted to 32-bit)
!
!OUTPUT PARAMETERS:
!    DESTRIPE    Array of destriped MODIS 1KM scaled integers (converted to 32-bit)
!    ERRFLAG     Error flag (zero: no errors)
!
!REVISION HISTORY:
!    Created by Liam.Gumley@ssec.wisc.edu
!
!TEAM-UNIQUE HEADER:
!    Copyright(C) 2004, University of Wisconsin-Madison MODIS Group
!
!END

IMPLICIT NONE

!- NOTE: ZERO-BASED ARRAY INDEXING IS USED THROUGHOUT!
      
!- Parameters
integer :: lines_1km, max_scan
parameter (lines_1km = 10, max_scan = 1000)

!- Input arguments
integer, intent(in)  :: npixel, nscan, ref, mir(0:nscan)
integer, intent(in)  :: image(0:npixel - 1, 0:nscan * lines_1km - 1)
integer, intent(out) :: destripe(0:npixel - 1, 0:nscan * lines_1km - 1)
integer, intent(out) :: errflag

!- Local variables
integer :: nvalid, ref_ind, npair, i, j, det, row(0:max_scan - 1), hist(0:32767)
integer :: count, ngood, sum, loc, lut(0:32767, 0:19), n
integer :: median_old, median_new, median_del
real :: edf(0:32767, 0:19)
real :: xtab(0:32767), ytab(0:32767), xint(0:32767), yint(0:32767)

!- External functions
integer :: modis_edf_median, modis_edf_valid
external   modis_edf_median, modis_edf_valid

!- Get the number of valid values in the input image
nvalid = modis_edf_valid(size(image), image)

!- If the number of valid values is less than one scan, return
if (nvalid <= (lines_1km * npixel)) then
  errflag = 1
  return
endif

!- Initialize the EDF array
edf = 0.0

!- Get index for reference detector on mirror side zero
ref_ind = ref
if (mir(ref) .eq. 1) ref_ind = ref + 10

!- Create the EDF for each detector
do det = 0, 19

  !- Compute row indices for this detector
  loc = nscan / 2 - 1
  row(0 : loc) = (/((i * 20) + det, i = 0, loc)/)
  
  !- Handle an odd number of scans
  if (mod(nscan, 2) .eq. 1) then
    loc = nscan / 2
    if (det <= 9) then
      row(loc) = row(loc - 1) + 20
    else
      row(loc) = row(loc - 1) - 20
    endif
  endif

  !- Compute the histogram for this detector
  hist = 0
  ngood = 0
  do j = 0, loc
    do i = 0, 1353
      count = image(i, row(j))
      if (count >= 0 .and. count <= 32767) then
        hist(count) = hist(count) + 1
        ngood = ngood + 1
      endif
    enddo
  enddo

  !- Compute the EDF for this detector
  sum = 0
  do i = 0, 32767
    sum = sum + hist(i)
    edf(i, det) = real(sum) / real(ngood)
  enddo

enddo

!- Create the lookup table (LUT) which maps each detector to the reference detector
xtab = edf(:, ref_ind)
ytab = (/(float(i), i = 0, 32767)/)
n = 32768
do det = 0, 19
  if (det .ne. ref_ind) then
    xint = edf(:, det)
    call interp(n, xtab, ytab, n, xint, yint)
    lut(:, det) = nint(yint)
  endif
enddo

!- Set the output image equal to the input image
destripe = image

!- Apply the LUT to all detectors except the reference detector
do det = 0, 19

  if (det /= ref_ind) then
  
    !- Compute row indices for this detector
    loc = nscan / 2 - 1
    row(0 : loc) = (/((i * 20) + det, i = 0, loc)/)
    if (mod(nscan, 2) == 1 .and. det <= 9) then
      loc = nscan / 2
      row(loc) = row(loc - 1) + 20
    endif
    
    !- Apply the LUT to valid input values
    do j = 0, loc
      where (image(:, row(j)) >= 0 .and. image(:, row(j)) <= 32767)
        destripe(:, row(j)) = lut(image(:, row(j)), det)
      endwhere
    enddo
    
  endif
  
enddo

!- Replace bad destriped values
where (destripe < 0 .or. destripe > 32767) destripe = image

!- Set median value of destriped image to median value of original image
median_old = modis_edf_median(size(image), image)
median_new = modis_edf_median(size(image), destripe)
median_del = median_new - median_old
if (median_del /= 0) then
  where (destripe >= 0 .or. destripe <= 32767) destripe = destripe - median_del
endif

! Set error flag
errflag = 0

END SUBROUTINE MODIS_EDF_DESTRIPE

!-------------------------------------------------------------------------------

FUNCTION MODIS_EDF_MEDIAN(N, IMAGE)

!F90
!
!DESCRIPTION:
!    Compute the median of an array of MODIS 1KM scaled integers.
!
!INPUT PARAMETERS:
!    N        Number of elements in input array
!    IMAGE    Array of MODIS 1KM scaled integers (converted to 32-bit)
!
!OUTPUT PARAMETERS:
!    MODIS_EDF_MEDIAN    Median value
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
integer, intent(in)  :: n
integer, intent(in)  :: image(n)
integer :: modis_edf_median

!- Local variables
integer :: hist(0:32767), i, count, sum

!- Compute histogram
hist = 0
do i = 1, n
  count = image(i)
  if (count >= 0 .and. count <= 32767) hist(count) = hist(count) + 1
enddo

!- Compute cumulative sum, and return when median is found
sum = 0
do i = 0, 32767
  sum = sum + hist(i)
  modis_edf_median = i
  if (sum >= (n / 2)) return
enddo

END FUNCTION MODIS_EDF_MEDIAN
!-------------------------------------------------------------------------------

FUNCTION MODIS_EDF_VALID(N, IMAGE)

!F90
!
!DESCRIPTION:
!    Count the number of valid values in an array of MODIS 1KM scaled integers.
!
!INPUT PARAMETERS:
!    N        Number of elements in input array
!    IMAGE    Array of MODIS 1KM scaled integers (converted to 32-bit)
!
!OUTPUT PARAMETERS:
!    MODIS_EDF_VALID    Number of valid values in the range 0 to 32767 inclusive
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
integer, intent(in)  :: n
integer, intent(in)  :: image(n)
integer :: modis_edf_valid

!- Local variables
integer :: sum, i

!- Compute number of valid values
sum = 0
do i = 1, n
  if (image(i) >= 0 .and. image(i) <= 32767) sum = sum + 1
enddo

!- Set return value
modis_edf_valid = sum

END FUNCTION MODIS_EDF_VALID
