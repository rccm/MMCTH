      INTEGER FUNCTION GET_COEFF(PCF_NUM, MAX_BANDS, MAX_COLS, 
     &  HEADER, PLATFORM, DATA_FLAG, NUM_BANDS, NUM_COLS,
     &  BAND_DATA, COEFF_DATA, ERROR_TEXT)

c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c     Read an ASCII text coefficient file. Currrently used to read
c     the MOD07 radiance bias and detector mask files.
c
c     The format of a coefficient file is as follows:
c
c     Line 1: Text header
c     Line 2: Satellite name (TERRA or AQUA)
c     Line 3: Number of bands, number of data columns, data flag
c     Subsequent lines: Band number, data column 1, data column 2, ...
c
c     The contents of a sample file are shown below:
c
c
c !INPUT PARAMETERS:
c     PCF_NUM        PCF number for coefficient file
c     MAX_BANDS      Maximum number of bands in coefficient file
c     MAX_COLS       Maximum number of data columns in coefficient file
c
c !OUTPUT PARAMETERS:
c     GET_COEFF      Return status code
c                    ( 0 = Success)
c                    (-1 = Error opening input file)
c                    (-2 = Error reading header text)
c                    (-3 = Error reading platform name)
c                    (-4 = Error reading number of bands and columns)
c                    (-5 = Error reading band and coefficient data)
c                    (-6 = Error closing input file)
c     HEADER         Header string
c     PLATFORM       Platform name string
c     DATA_FLAG      Flag (for bias, indicates if slope and intercept
c                    are for TPW: flag = 1, or BT: flag = 2.  Any other
c                    flags will mean no radiance bias is applied)
c     NUM_BANDS      Number of bands in coefficient file
c     NUM_COLS       Number of data columns in coefficient file
c     BAND_DATA      Array of band numbers
c     COEFF_DATA     Array of coefficient data
c     ERROR_TEXT     Error text string
c                    (Blank if successful; explanatory text if failure)
c
c !REVISION HISTORY:
c     Created 12 March 2003
c
c     SWS added DATA_FLAG 13 March 2003
c
c !TEAM-UNIQUE HEADER:
c     Developed by Liam Gumley, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------
      
      implicit none

c ... Include files for SDP toolkit I/O
      include 'PGS_IO.f'
      include 'PGS_SMF.f'

c ... Input arguments
      integer pcf_num, max_bands, max_cols
      
c ... Output arguments
      character*(*) header, platform
      integer num_bands, num_cols
      integer band_data(max_bands)
      real coeff_data(max_cols, max_bands)
      character*(*) error_text

c ... Local variables
      integer record_length, file_version, lun, status
      integer iostat
      integer row, col, data_flag

c ... External functions for SDP toolkit I/O
      integer pgs_io_gen_openf, pgs_io_gen_closef
      external pgs_io_gen_openf, pgs_io_gen_closef

c...  Initialize the output arguments
      do col = 1, len(header)
        header(col:col) = ' '
      enddo
      do col = 1, len(platform)
        platform(col:col) = ' '
      end do
      num_bands = -1
      num_cols = -1
      do row = 1, max_bands
        band_data(row) = -9999.0
      enddo
      do row = 1, max_bands
        do col = 1, max_cols
          coeff_data(col, row) = -9999.0
        enddo
      enddo
      
c ... Open input file
c ... (read-only sequential formatted, i.e., ASCII text)
      record_length = 1
      file_version = 1
      lun = -1
      status = pgs_io_gen_openf(pcf_num, PGSD_IO_GEN_RSEQFRM,
     &  record_length, lun, file_version)
      if (status .ne. PGS_S_SUCCESS) then
        get_coeff = -1
        error_text = 'Error opening input file'
        return
      endif

c ... Read header text
      read(lun, '(a)', iostat=iostat) header
      if (iostat .ne. 0) then
        get_coeff = -2
        error_text = 'Error reading header text'
        return
      endif

c ... Read platform name  
      read(lun, '(a)', iostat=iostat) platform
      if (iostat .ne. 0) then
        get_coeff = -3
        error_text = 'Error reading platform name'
        return
      endif
      
c ... Read number of bands and columns
      read(lun, *, iostat=iostat) num_bands, num_cols, 
     &                            data_flag
      if (iostat .ne. 0) then
        get_coeff = -4
        error_text = 'Error reading number of bands and columns'
        return
      endif

c ... Read band and coefficient data
      do row = 1, num_bands
        read(lun, *, iostat=iostat)
     &    band_data(row), (coeff_data(col, row), col = 1, num_cols)
        if (iostat .ne. 0) then
          get_coeff = -5
          error_text = 'Error reading band and coefficient data'
          return
        endif
      enddo           
     
c ... Close input file
      status = pgs_io_gen_closef(lun)
      if (status .ne. PGS_S_SUCCESS) then
        get_coeff = -6
        error_text = 'Error closing input file'
        return
      endif

c ... Set return value of function
      get_coeff = 0
      
      END
