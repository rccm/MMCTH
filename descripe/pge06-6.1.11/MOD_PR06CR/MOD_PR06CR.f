      program MOD_PR06Create

      implicit none

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Driver program for MODIS Cloud Product HDF File Create Routine.
c
c!Input Parameters:
c    You must have a L1B 1km file to get the correct number of
c    lines and elements for the output product file.  This is
c    done through the PCF file.
c
c!Output Parameters:
c    Output consists of the creation of the MOD06 cloud product
c    skeleton HDF file.  The file name is taken from the PCF
c    file.
c
c!Revision History:
c $Log: MOD_PR06CR.f,v $
c Revision 1.1.1.1  2005/02/22 17:15:54  gumley
c MODIS Operational Code
c
c Revision 1.1  1998/10/29 22:51:06  kis
c Created
c
c Revision 1.6  1998/10/12 18:36:53  gumley
c Latest baselined V2 source from Rich Hucek
c
c 10/12/98   Kathy Strabala  UW Madison
c Stripped out from Driver program of MOD06_CT.
c
c Revision 1.3  1998/05/15 19:42:58  gumley
c Removed calls to multiple versions of mod06_open_files and mod06_close_files
c
c!Team-unique Header:
c
c
c!End
c-----------------------------------------------------------------------

C#######################################################################
C     DECLARATIONS
C#######################################################################

c ... Include files

      include 'mod06_create.inc'
      include 'mod06_pcfnum.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'

c ... Local variables

      character*32 routine_name
      parameter ( routine_name = 'MOD_PR06CR' )

      character*2 out_number
      character*( PGSd_PC_FILE_PATH_MAX ) out_file

      integer file_version, rtn
      integer out_handle( MODFILLEN )

      integer nscans, npixels,
     &  out_lines_1km, out_elements_1km,
     &  out_lines_5km, out_elements_5km

c ... External functions
      external clmfil

      integer  pgs_pc_getreference
      external pgs_pc_getreference


C#######################################################################
C     MAIN PROGRAM
C#######################################################################

c ... Get number of scans and pixels for 1km bands from the L1B file

      call mod06_Create_Get_Metadata( nscans, npixels )

c ... Get output file name from process control file

      file_version = 1
      out_file = ' '
      rtn = PGS_S_SUCCESS - 1
      rtn = pgs_pc_getreference( out_pcfnum, file_version, out_file )
      if( rtn .ne. PGS_S_SUCCESS ) call message( routine_name,
     &  'Error getting output filename from PCF.' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )

c ... Create the MOD06 output HDF file, and open it using MAPI

      out_number = '06'
      out_lines_1km = nscans * max_line
      out_elements_1km = npixels
      out_lines_5km = out_lines_1km / isamp
      out_elements_5km = out_elements_1km / isamp
      call setup_create_hdf( out_number, out_file, out_lines_1km,
     &  out_elements_1km, out_lines_5km, out_elements_5km, out_handle )

c ... Close the output file using MAPI routine.

      rtn = MAPIOK - 1
      rtn = clmfil( out_handle )
      if( rtn .eq. MAPIOK ) then
         call message( routine_name,
     &  'Sucessfully closed output MOD06_L2 file.', 0, 3 )
      else
         call message( routine_name,
     &  'Error closing output MOD06_L2 file.' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 2 )
      endif

c ... Return exit code 0 to shell

      call exit( 0 )

      end
