      SUBROUTINE OPENR_TEMP( PCFNUM, RECLEN, LUN )
      
C-----------------------------------------------------------------------
C
C !F77
C
C !DESCRIPTION:
C      Open a binary direct access temporary file for read-only
C      using the PGS toolkit.
C
C !INPUT PARAMETERS:
C      PCFNUM    File number in Process Control File.
C      RECLEN    Record length (bytes).
C
C !OUTPUT PARAMETERS:
C      LUN       Logical unit number assigned to the opened file.
C
C !REVISION HISTORY:
c 01/14/98 fhliang
c filled in prolog.
c
C      04-AUG-1997 Liam Gumley CIMSS/SSEC
C          Created using OPENR.f as a template.
C
c!Team-Unique Header:
c    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
C !END
C
C-----------------------------------------------------------------------

      implicit none

      include 'PGS_IO.f'
      
c ... input arguments

      integer pcfnum, reclen, lun

c ... open the file for read,direct,unformatted
      
      call modis_io_gen_temp_openf( pcfnum, pgsd_io_gen_rdirunf,
     &  reclen, lun )

      end
