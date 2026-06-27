      Integer Function SetGlobalAttr(sdid, lines10km, elements10km) 

C----------------------------------------------------------------------
C !F77
C
C !Description:  SetGlobalAttr sets two global attributes that 
C                specify the size characteristics of the current 
C                MODIS granule. Two more global attributes identify the
C                Aerosol product and define the unpacking scheme in 
C                terms of slope and offset. 
C
C                     
C
C !INPUT PARAMETERS:
C
C    integer sdid          the HDF SD Interface ID used for access to 
C                          data objects at file level. 
C
C    integer lines10km     number of the scans in L1B data granule 
C
C    integer elements10km  maximum number of the 10km elements in a
C                          L1B data granule
C
C
C !OUTPUT PARAMETERS:
C    
C    NONE
C
C
C !REVISION HISTORY:
c  Revision 1.2  1997/02/25  18:55:14  vlin
c  Changed from function to a subroutine
C
C !TEAM-UNIQUE HEADER:
C
C    This software was developed by the MODIS Science Data Support Team
C    (SDST) for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS
C
C    Written by Richard Hucek
C
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C
C    rhucek@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C    There is no check to identify the unpredictable value of an
C    undefined string.  Consequently, users must take care to initialize
C    all string variables before passing them to CUT_NAME.
C
C    Externals:
C
C       Functions and Subroutines:
C          message                   (src_UW)
C          SetGlobalAttr             (science code)
C          sfsattr                   (dffunc.inc)
C          sfscatt                   (dffunc.inc)
C          strlen                    (c )
C
C
C       Named Constant:
C          DFNT_CHAR                 (hdf.inc)
C          DFNT_INT32                (hdf.inc)
C          FAIL                      (?)
C          FALSE                     (?)
C          SUCCEED                   (?)
C
C
C !END
C-----------------------------------------------------------------------

      implicit none
c     include 'Atmos_ECSMET.inc'
      include 'hdf.inc'
      include 'dffunc.inc'

c Function argument declarations
      integer  elements10km, lines10km, sdid

c local variables
      character*256  attrname, title
      character*2048 descrip
      integer count, numbertype, rtn_flag, sl, srtn, strlen
      logical error_flag

c function declarations
c     integer  sfscatt, sfsattr 
c     external sfscatt, sfsattr


C------------------------
C Initialization
C------------------------

      rtn_flag      =  FAIL
      error_flag    = .FALSE.
      SetGlobalAttr = -1      


c-----------------------------------------------
c ... set attribute 'Number_of_Instrument_Scans'
      srtn=sfsattr(sdid,'Number_of_Instrument_Scans',DFNT_INT32,1,lines10km)

      if (srtn .ne. 0) then
         call message( 'setup_create_hdf',
     &                 'Error defining dimensions for output file: swdefdim ' //
     &                 '[OPERATOR ACTION: Notify SDST.]', 0, 2)
      endif


c------------------------------------------------
c     set attribute 'Maximum_Number_of_1km_Frames' 

      srtn=sfsattr(sdid,'Maximum_Number_of_1km_Frames',DFNT_INT32,1,elements10km)

      if (srtn .ne. 0) then
         call message( 'setup_create_hdf',
     1                 'Error defining dimensions for output file: swdefdim ' //
     2                 '[OPERATOR ACTION: Notify SDST.]', 0, 2)
      endif


c--------------------------
c     set attribute 'title'

      attrname   = 'title'
      numbertype = DFNT_CHAR

      descrip    = char(10)
     1             // ' MODIS HDF File Specification MOD04_L2: MODIS Level 2 Aerosol ' // char(10) 
     2             // ' Land and Ocean Product                                       ' // char(10)
     3             // ' '

      sl         = strlen(descrip)
      count      = sl
      srtn       = sfscatt(sdid, attrname, numbertype, count, descrip)


c------------------------------------------
c     set attribute 'Slope_and_Offset_Usage'

      attrname   = 'Slope_and_Offset_Usage'
      numbertype = DFNT_CHAR

      descrip    = char(10)
     1 // ' The local SDS scale_factor and add_offset attributes are used for the      ' // char(10)
     2 // ' conversion of stored integer data to geophysical floating point numbers.   ' // char(10)
     3 // ' The implementation follows conventional HDF usage (See HDF Users Guide):   ' // char(10)
     4 // '                                                                            ' // char(10)
     5 // '       float value = scale_factor*(stored integer - add_offset)             ' // char(10)
     6 // '                                                                            ' // char(10)
     8 // ' The unit of the derived floating point value is indicated in the ''units'' ' // char(10)
     9 // ' local attribute which is also provided.                                    ' // char(10)
     1 // ' '

      sl         = strlen(descrip)
      count      = sl
      srtn       = sfscatt(sdid, attrname, numbertype, count, descrip)

      If (.not.error_flag) SetGlobalAttr = SUCCEED

      return
      end
