      subroutine mod06_close_files( cube_handle, geo_handle, cld_handle,
     &                        out_handle, h_eco1 )

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Close mod06 output files, output files.
c
c!Input Parameters:
c    cube_handle MAPI hdf and SDS information on L1B file
c    geo_handle  MAPI hdf and SDS information on geolocation file
c    cld_handle  MAPI hdf and SDS information on cloud mask file
c    out_handle  MAPI hdf and SDS information on output  file
c    h_eco1      File handle for ecosystem file - 1km global Olson
c
c!Output Parameters:
c   None
c
c!Revision History:
c $Id: mod06_close_files.f,v 1.7 1999/05/24 15:40:44 kis Exp $
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

      implicit none
      save
      
      include 'PGS_SMF.f'
      include 'PGS_IO.f'
      include 'mapi.inc'
      include 'mod06uw_pcfnum.inc'
      include 'mod06uw_debug.inc'
      
c ... arguments

      integer cube_handle(MODFILLEN), geo_handle(MODFILLEN), 
     *        cld_handle(MODFILLEN), out_handle(MODFILLEN),
     *        h_eco1

c ... local variables

c     Local scalars
      integer rtn,ier

c ... set program name for error messaging
      character*32 routine_name
      parameter ( routine_name = 'mod06_close_files' )

c ... external functions

      integer  pgs_io_gen_closef
      external pgs_io_gen_closef, message
            
c ... Close the 1km L1B file

      rtn = clmfil(cube_handle)
      if( rtn .eq. 0 ) then
        call message( routine_name, 
     &  'Success closing scan cube 1km data file ', 0, 3 )
      else
        call message( routine_name,
     &  'Error closing scan cube 1km data file.' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
      endif

c ... Close the geolocation file

      rtn = clmfil(geo_handle)
      if( rtn .eq. 0 ) then
        call message( routine_name,
     &  'Success closing geolocation data file ', 0, 3 )
      else
        call message( routine_name,
     &  'Error closing geolocation data file.' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
      endif

c ... Close the cloud mask input file

      rtn = clmfil(cld_handle)
      if( rtn .eq. 0 ) then
        call message( routine_name,
     &  'Success closing cloud mask data file ', 0, 3 )
      else
        call message( routine_name,
     &  'Error closing cloud mask data file.' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
      endif

c ... Close ecosystem files

c ... First the 1km NA map

      ier = pgs_io_gen_closef(h_eco1)
      if( ier .eq. PGS_S_SUCCESS ) then
       call message( routine_name, 
     &  'Success closing 1km NA eco. file using pgs_io_gen_closef ',
     &  0, 3 )
      else
       call message( routine_name,
     &  'Error closing 1km NA eco. file using pgs_io_gen_closef' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
      endif

c ... Close the output file

      rtn = clmfil(out_handle)
      if( rtn .eq. 0 ) then
        call message( routine_name,
     &  'Success closing output MOD_PR06 file ', 0, 3 )
      else
        call message( routine_name,
     &  'Error closing output MOD_PR06 file.' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
      endif

      end
