      subroutine read_reynsst ( npoints_x, npoints_y, 
     &                          year, month, day, hour, sst, error,
     &                          ice, file_format, sst_success ) 
     
      implicit none

c ... Include files 
      include 'PGS_IO.f'
      include 'PGS_SMF.f'
      include 'Atmos_AncData.inc'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:  Retrieves an array of global sea surface temperature  
C                (sst) data from the NCEP Reynolds sst data product.
C                Both formatted (oi.mean.bias) and unformatted (oisst)
C                products are supported.  
C
C
C !INPUT PARAMETERS:
C  Integer       npoints_x   - number of points around a latitude circle
C  Integer       npoints_y   - number of grid points along a meridian  
C
C !OUTPUT PARAMETERS:
C  character*(*) file_format - one of either 'unk', 'fmt' or 'unf'
C                        
C  Integer       year        - data year 
C  Integer       month       - data month
C  Integer       day         - data day
C  Intger        hour        - set to 0 
C
C  Logical       sst_success - .true. if able to retrieve sst data; 
C                              otherwise .false.
C
C  Real          sst         - global array of REYNSST sea surface
C                              temperature data in degree Celsius 
C  Real          error       - global array of REYNSST sea surface
C                              temperature error estimates
C  Byte          ice         - global fraction of ice concentration
C
C !REVISION HISTORY:
C
C
C !TEAM-UNIQUE HEADER:
C
C    This software wasdeveloped by the MODIS Science Data Support Team
C    together with MODIS Group, CIMSS/SSEC, UW-Madison for the National
C    Aeronautics and Space Administration, Goddard Space Flight Center,
C    under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by Rich Frey, MODIS Group, CIMSS/SSEC, UW-Madison 
C        and
C    Rich Hucek, Decision Systems Technologies, INC 
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C    rhucek@ltpmail.gsfc.nasa.gov
C
C
C !DESIGN NOTES:  
C
C    The NCEP Reynolds sst file format (formatted or unformatted) is
C    identified by comparing the size of the product staged for input 
C    against known, unique and fixed values for the two format types.
C    If a match is made, processing continues.  Otherwise, read_reynsst
C    aborts with a fatal error.
C
C    If a read error occurs, read_reynsst returns control to the 
C    calling program and sets the subroutine logical status flag,
C    sst_success, to .false.  
C
C
C  Externals:
C
C    Functions and Subroutines:
C        modis_io_gen_closef       (MODIS_IO_GEN_CLOSE.f in src_L2)
C        pgs_io_gen_openf_BE          (libPGSTK.a)
C        pgs_pc_getfilesize        (libPGSTK.a)
C
C    Named Constant:
C        LUN_REYNSST               (Atmos_AncData.inc in src_L2)
C        SIZE_OF_REYNSST_FMT       (Atmos_AncData.inc in src_L2)
C        SIZE_OF_REYNSST_UNF       (Atmos_AncData.inc in src_L2)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C
C !END
C-----------------------------------------------------------------------

c-----Declaration of function arguments
      character*(*) file_format
      integer       npoints_x, npoints_y, year, month, day, hour
      logical       sst_success
      real          sst(0:npoints_x-1, 0:npoints_y-1),
     +              error(0:npoints_x-1, 0:npoints_y-1)
      byte          ice(0:npoints_x-1, 0:npoints_y-1)


c-----Declaration of local variables 
      character*160 errmsg
      integer       filesize, header(8), i, ios, j, level, lun, pcfnum,
     1              reclen, status, version 


c-----Declaration of functions 
      integer  pgs_io_gen_openf_BE, pgs_pc_getfilesize
      external pgs_io_gen_openf_BE, pgs_pc_getfilesize

c-----------------------------------------------------------------------
c     INITIALIZATION
c-----------------------------------------------------------------------
      sst_success = .false.
      file_format = 'unk'


c-----------------------------------------------------------------------
c     DETERMINE FILE FORMAT 
c-----------------------------------------------------------------------
      pcfnum  = LUN_REYNSST
      reclen  = 1
      version = 1

      status  = pgs_pc_getfilesize(pcfnum, version, filesize) 

c ... can't get file size.  Process w/o NCEP sst
      if ( status .ne. PGS_S_SUCCESS ) then
         level = 1
         write( errmsg,'(''Unable to get size of sst file on LUN '', 
     1                     i12)') pcfnum
         call message( 'read_reynsst', errmsg // 
     1      ' [OPERATOR ACTION: Contact SDST]', status, level )

c ... file is formatted NCEP sst
      elseif (filesize .eq. SIZE_OF_REYNSST_FMT) then
         file_format = 'fmt'
         status = pgs_io_gen_openf_BE( pcfnum, PGSd_IO_Gen_RSeqFrm, 
     1                              reclen, lun, version )
      
         if ( status .ne. PGS_S_SUCCESS ) then
            level = 1
            write( errmsg,'(''Unable to open formatted sst '',
     1         ''product on LUN '',i12)') pcfnum
            call message( 'read_reynsst', errmsg //
     1         ' [OPERATOR ACTION: Contact SDST]', status, level )
         endif

c ... file is unformatted NCEP sst
      elseif (filesize .eq. SIZE_OF_REYNSST_UNF) then
         file_format = 'unf'
         status = pgs_io_gen_openf_BE( pcfnum, PGSd_IO_Gen_RSeqUnf,
     1                              reclen, lun, version )

         if ( status .ne. PGS_S_SUCCESS ) then
            level = 1
            write( errmsg,'(''Unable to open unformatted sst '',
     1         ''product on LUN '',i12)') pcfnum
            call message( 'read_reynsst', errmsg //
     1         ' [OPERATOR ACTION: Contact SDST]', status, level )
         endif

c ... file is corrupted or not NCEP sst product. QUIT immediately 
      else
         level   =  2
         status  = -1

         write( errmsg, '(''Size of sst file ('', i12, '' bytes) '',
     1      ''on LUN '', i12, '' matches '', A1, ''neither  '', 
     2      ''formatted nor unformatted REYNSST sst product.'')' ) 
     3      filesize, pcfnum, char(10)

         call message( 'read_reynsst', errmsg // char(10) 
     1      // '[OPERATOR ACTION: Check for corrupted sst ' 
     2      // 'file; Contact SDST]', status, level )
      endif


c-----------------------------------------------------------------------
c     READ SST DATA 
c-----------------------------------------------------------------------
      if ( status .eq. PGS_S_SUCCESS ) then

c ...    Initialize ice output in case we are reading an ascii file
         do i = 0, 359 
            do j = 0, 179
               ice(i,j) = 0
            enddo
         enddo

         if (file_format .eq. 'fmt') then
            read( lun, '( 8i5 )', iostat = ios ) header
            write(errmsg, '('' Reading sst ASCII data PCF#'',i12)' ) 
     1      pcfnum
            call message( 'read_reynsst', errmsg //
     1          ' [OPERATOR ACTION: Contact SDST]', 0, 3 )
         else if(file_format .eq. 'unf') then
            read( lun, iostat = ios ) header
         end if

         if ( ios .ne. 0 ) then
            level = 1
            write( errmsg,
     1         '(''Error reading sst header PCF#'',i12)') pcfnum
            call message( 'read_reynsst', errmsg //
     1         ' [OPERATOR ACTION: Contact SDST]', ios, level )
         else
            year  = header( 1 )
            month = header( 2 )
            day   = header( 3 )
            hour  = 0
         endif

         if (file_format .eq. 'fmt') then
            read( lun, '( 20f4.2 )', iostat = ios ) sst
         else if (file_format .eq. 'unf') then
            read( lun, iostat = ios ) ( (sst(i,j), i=0,359), j=0,179 )
            read( lun, iostat = ios ) ( (error(i,j), i=0,359), j=0,179 )
            read( lun, iostat = ios ) ( (ice(i,j), i=0,359), j=0,179 )
         end if

         if ( ios .ne. 0 ) then
            level = 1
            write( errmsg, '(''Error reading sst data PCF#'',i12)' ) 
     1      pcfnum
            call message( 'read_reynsst', errmsg //
     1          ' [OPERATOR ACTION: Contact SDST]', ios, level )
         else
            sst_success = .true.
         endif

         call modis_io_gen_closef( pcfnum, lun )
      endif

      end
