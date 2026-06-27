      subroutine mod35_out(s_pixels,nlin,klin,lines_in_edge,
     +                     modfil_mod35,bitarray,qa_bitarray)

      IMPLICIT NONE
      SAVE

c     Common block used for writing debug output statements
      common / bug / debug, h_output

C!F77------------------------------------------------------------------
C
C!DESCRIPTION:  Write Scientific Data Set (SDS) for cloud
C               mask algorithm. The process is complicated
c               by the way the cloud mask is processed. The
c               method used is a nlcntx x necntx sliding box,
c               so more than one line of data is being used
c               at a time.  This program adds a line of data
c               to the holding array, then writes it out
c               once the tenth line has been filled.
c
C!INPUT PARAMETERS:
c s_pixels      Context of pixels in current scan
C nlin          Total number of lines to process
c klin          Counts number of lines output to bit file
c lines_in_edge Number of lines outside of processing region
C modfil_mod35 	Array used to reference MODIS HDF file
C bitarray      Byte array containing line of cloud mask output
C qa_bitarray   Byte array containing line of cloud mask qa output
C
C!OUTPUT PARAMETERS:
c None
C
C!REVISION HISTORY:
C
C!TEAM-UNIQUE HEADER:
C
C!REFERENCES and CREDITS:
C
C!DESIGN NOTES:
C       Functions:
C          CRMAR                           (mapi.inc)
C          PMAR                            (mapi.inc)
C
C       Named Constant:
C          MAPIOK                          (mapi.inc)
c          MODFILLEN                       (mapi.inc)
C
C!END------------------------------------------------------------------

      INCLUDE 'mapi.inc'
      INCLUDE 'global.inc'

c ... scalar arguments
      integer nlin,klin,lines_in_edge
c ... array arguments
      integer modfil_mod35(MODFILLEN),s_pixels(nlcntx)
      byte bitarray(npixel,6),qa_bitarray(npixel,10)


c ... Parameter statements - dimensions for product and qa arrays
      integer byte_dim,qa_dim
      parameter(byte_dim = 6, qa_dim = 10)

c ... number of lines to collect beforing writing to output file
      integer  lines_per_write
      parameter(lines_per_write = scans_cube)

c ... modis api toolkit declarations
      integer dims(3),start(3)

c ... local scalars
      character*100 info_string
      character*72 arrayname
      integer buf_line_counter,line_counter,ipxl,rtn,ibyte,h_output,
     +        debug,nele,i,iline,buf_offset,buf_index

c ... local arrays
      byte sstore1(npixel*scans_cube*byte_dim)
      byte sstore2(qa_dim*npixel*scans_cube)

c ... initialize line counter
      data line_counter /0/


c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Processing subroutine mod35_out '',/)')
      endif
c ................................................................

C     Initialization
c     Need to know how many pixels to store for this scan
      if (klin .le. lines_in_edge) then
         i = klin
      else if (klin .gt. nlin-lines_in_edge) then
         i = nlcntx
      else
         i = ((nlcntx-1) / 2) + 1
      endif

c     Enter hdf file array dimensions
      nele = s_pixels(i)

      line_counter = line_counter + 1
      write(info_string,'(I10)') line_counter

c     Determine current line in work buffer (buf), and transfer
c     data from bitarray to work buffer.
      buf_line_counter = mod(line_counter-1,lines_per_write) + 1

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''lc1 lc2 nele nlin'',4i10)') line_counter,
     +          buf_line_counter,nele,nlin
      endif
c ................................................................
      iline = buf_line_counter
      do ibyte = 1 , byte_dim
          buf_offset = (ibyte-1)*nele*lines_per_write
          do ipxl = 1 , nele
            buf_index = buf_offset + (iline-1) * nele + ipxl
            sstore1(buf_index) = bitarray(ipxl,ibyte)
          end do
      end do

c ... Now for the qa array
      do ipxl = 1 , nele
         buf_offset = (iline-1)*nele*qa_dim
         do ibyte = 1 , qa_dim
            buf_index = buf_offset + (ipxl-1) * qa_dim + ibyte
            sstore2(buf_index) = qa_bitarray(ipxl,ibyte)
         end do
      end do

c ... if the number of lines in buffer is 10, then write out
      if (buf_line_counter.eq.lines_per_write) then

         arrayname = 'Cloud_Mask'
         start(1) = 0
         start(2) = lines_per_write*(line_counter/lines_per_write-1)
         start(3) = 0
         dims(1) = nele
         dims(2) = lines_per_write
         dims(3) = byte_dim
         rtn = pmar(modfil_mod35,arrayname,' ',start,dims,sstore1)
         if (rtn .ne. MAPIOK) then
           call message( 'mod35_out',
     &     'Could not write product bits to output file using pmar' //
     &     ' [OPERATOR ACTION: Check system resources]',
     &     0, 2 )
         endif

         arrayname = 'Quality_Assurance'
         dims(1) = qa_dim
         dims(2) = nele
         dims(3) = lines_per_write
         start(1) = 0
         start(2) = 0
         start(3) = lines_per_write*(line_counter/lines_per_write-1)
         rtn = pmar(modfil_mod35,arrayname,' ',start,dims,sstore2)
         if (rtn .ne. MAPIOK) then
           call message( 'mod35_out',
     &     'Could not write qa_bits to output file using pmar' //
     &     ' [OPERATOR ACTION: Check system resources]',
     &     0, 2 )
         endif


c ... If we have come to the last line, then you'll want to
c ... write out anyway.
      else if (line_counter.eq.nlin) then

         arrayname = 'Cloud_Mask'
         dims(1) = nele
         dims(2) = buf_line_counter
         dims(3) = byte_dim
         start(1) = 0
         start(2) = lines_per_write * (line_counter/lines_per_write)
         start(3) = 0
         rtn = pmar(modfil_mod35,arrayname,' ',start,dims,sstore1)
         if (rtn .ne. MAPIOK) then
           call message( 'mod35_out',
     &     'Could not write product bits to output file using pmar' //
     &     ' [OPERATOR ACTION: Check system resources]',
     &     0, 2 )
         endif

         arrayname = 'Quality_Assurance'
         dims(1) = qa_dim
         dims(2) = nele
         dims(3) = buf_line_counter
         start(1) = 0
         start(2) = 0
         start(3) = lines_per_write*(line_counter/lines_per_write)
         rtn = pmar(modfil_mod35,arrayname,' ',start,dims,sstore2)
         if (rtn .ne. MAPIOK) then
           call message( 'mod35_out',
     &     'Could not write qa_bits to output file using pmar' //
     &     ' [OPERATOR ACTION: Check system resources]',
     &     0, 2 )
         endif
      end if

      return
      end
