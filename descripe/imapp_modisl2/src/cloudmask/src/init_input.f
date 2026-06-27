      subroutine init_input(modfil_L1B_1km,modfil_L1B_250m,
     +                      modfil_mod35,modfil_Geo,
     +                      buf_size,buf_geo_size,buf_anc_size,
     +                      qa_bits,bitarray,qa_bitarray)
      implicit none
      save
c----------------------------------------------------------------------
C!F77
c
c!Description:
C     Routine for initializing variables used in the cloud mask
c     MODIS software.
C
c!Input parameters:
c modfil_L1b_1km   Level 1B hdf file array - 1km data
c modfil_L1b_250m  Level 1B hdf file array - 250m data
c modfil_mod35     Output cloud mask file hdf file array
c modfil_Geo       Geolocation file hdf array
c buf_size         Line and element granule size
c buf_geo_size     Line and element geolocation file size
c buf_anc_size     Line and element ancillary data file size
c qa_bits          10 byte array contining qa bit results
c bitarray         Array containing a line of cloud mask bit values
c qa_bitarray      Array containing line of 10 byte qa results
c
c!Output Parameters:
c None.
c
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c     MODFILLEN  - size of hdf file arrays for mapi - mapi.inc
c
c!END
c----------------------------------------------------------------------
      INCLUDE 'mapi.inc'
      INCLUDE 'global.inc'

c ... array arguments

      integer modfil_L1B_1km(MODFILLEN),modfil_L1B_250m(MODFILLEN),
     +        modfil_mod35(MODFILLEN),modfil_Geo(MODFILLEN),
     +        buf_size(2),buf_geo_size(2),
     +        buf_anc_size(2)
      byte bitarray(npixel,6),qa_bitarray(npixel,10),qa_bits(10)

c ... local scalars

      integer i,j

      do 10 i = 1 , MODFILLEN
         modfil_L1b_1km(i) = 0
         modfil_L1b_250m(i) = 0
         modfil_mod35(i) = 0
         modfil_Geo(i) = 0
  10  continue

c ... set max buf size of radiance data
      buf_size(1) = nx
      buf_size(2) = ny

c ... set max buf size of geolocation data
      buf_geo_size(1) = npixel
      buf_geo_size(2) = scans_cube
c ... set max buf size of ancillary data file
      buf_anc_size(1) = npixel
      buf_anc_size(2) = scans_cube

c ... initialize testbits and qa_bits output products
c     Initialize array which holds cloud mask bit results for each line
      do 20 i = 1 , npixel
        do 30 j = 1 , 10
           if (j .le. 6) bitarray(i,j) = 0
           qa_bitarray(i,j) = 0
 30     continue
 20   continue

c ... Initialize qa_bits right away. Will set some in file open routine
      do 40 i = 1 , 10
         qa_bits(i) = 0
 40   continue

      return 
      end

 


