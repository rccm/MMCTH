      SUBROUTINE MODIS_IO_GEN_TEMP_OPENF( FILE_PCF_INDEX, FILE_ACCESS, 
     &  RECORD_LENGTH, FILE_HANDLE )

C-----------------------------------------------------------------------
C
C !F77
C
C !DESCRIPTION:
C 	Open a temporary file using the SDPTK generic file open function.
C
C !INPUT PARAMETERS:
C 	FILE_PCF_INDEX  PCF index number assigned to the file.
C 	FILE_ACCESS     SDPTK file access type.
C 	RECORD_LENGTH   If direct access file, record length in bytes
C 			    (ignored for other file types).
C
C !OUTPUT PARAMETERS:
C 	FILE_HANDLE     FORTRAN logical unit number assigned to the file.
C
C !REVISION HISTORY:
c 01/14/98 fhliang
c filled in prolog.
c
C 	20-JUL-1998 Liam Gumley CIMSS/SSEC
C         Messaging is now handled by message.f
C         Added operation action message.
C
C 	13-NOV-1997 Liam Gumley CIMSS/SSEC
C 	    Changed first argument of call to PGS_IO_GEN_TEMP_OPENF from
C           PGSD_IO_GEN_ENDURANCE to PGSD_IO_GEN_NOENDURANCE to indicate
C           that this is a temporary file which only endures for the
C           duration of the calling PGE (as advised by Rich Hucek).
C
C 	04-AUG-1997 Liam Gumley CIMSS/SSEC
C 	    Created using MODIS_IO_GEN_OPENF.f as a template.
C
c!Team-Unique Header:
c    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
C !END
C
C-----------------------------------------------------------------------

      implicit none
      
c ... Included files
      
      include 'PGS_SMF.f'
      include 'PGS_PC.f'
      include 'PGS_PC_9.f'
      include 'PGS_IO.f'
      include 'PGS_IO_1.f'
      include 'PGS_MODIS_39500.f'
      
c ... Input arguments
                                             
      integer file_pcf_index, file_access, record_length, file_handle

c ... Local variables

      integer returnstatus
      character*10 pcf_string
      
c ... Declaration of functions 

      integer pgs_io_gen_temp_openf
      external pgs_io_gen_temp_openf
      
c ... Open file using sdptk generic file open function 

      returnstatus = 0
      returnstatus = pgs_io_gen_temp_openf( pgsd_io_gen_noendurance, 
     &  file_pcf_index, file_access, record_length, file_handle )

c ... Check return status of sdptk file open call and
c ... report success or error message to log file

      write( pcf_string, '(i10)' ) file_pcf_index
      
      if ( returnstatus .eq. pgs_s_success ) then
 
        call message( 'modis_io_gen_temp_openf',
     &  'Successfully opened temporary file PCF index ' // pcf_string,
     &  0, 3 )

      else
      
        call message( 'modis_io_gen_temp_openf',
     &  'Error opening temporary file PCF index ' // pcf_string //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 2 )
      
      endif

      end
