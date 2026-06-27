      subroutine fill_bit_line(nc,nmtests,nbands,bad_value,bad_geo,
     +                         snglnt,desert,testbits,qa_bits,bitarray,
     +                         qa_bitarray)

      implicit none
      save

c----------------------------------------------------------------------
c!F77 
c
c!Description:
c     Routine for placing results of the cloud mask product and
c     qa array into a line of data values.
c
c!Input parameters:
c nc            Current processing element
c nmtests       Number of tests applied to this pixel
c nbands        Number of bands successfully read for this pixel
c bad_value     Logical value indicating band radiance or reflectance 
c               value
c bad_geo       Logical variable flagging bad lat/long data
c snglnt        Logical variables where true indicates sun glint
c desert        Logical varibable whre true indicates desert ecosystem
c testbits      four-byte integer containing bit results
c qa_bits       Byte array containing qa bit results
c
c!Output Parameters:
c bitarray      Array containing line of 48 bit test results
c qa_bitarray   Array containing line of 10 byte qa results
c
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c----------------------------------------------------------------------

      include 'global.inc'

c     scalar arguments
      integer nc,nmtests,nbands
      logical bad_value,snglnt,desert,bad_geo

c     scalar arrays
      byte testbits(6),bitarray(npixel,6),qa_bits(10),
     +     qa_bitarray(npixel,10) 

c     local scalars 
      integer i,debug,h_output

c     external routines
      external set_bit,set_qa_bit

c ... Common statement for debug purposes
      common / bug / debug, h_output

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Within fill_bit_line routine '',/)')
        write(h_output,'(10x/,''Bad_value,Bad_geo, desert, sunglint= '',
     +        4L5,/)') bad_value, bad_geo, desert, snglnt
      endif
c .............................................................

c ... Fill in final pixel values before putting into line
c     array.  This is where the quality of the cloud mask
c     product is determined.  

c ... Now decide quality of pixel confidence based upon
c ...  number of tests used and processing path

c ... Next, how many tests and bands were used? If none
c ... set all output product bits to zero.
      if (nmtests .eq. 0 .or. nbands .eq. 0 .or. bad_geo) then
        do 100 i = 1 , 6
          testbits(i) = 0
  100   continue
c ...   and set the usable qa_bit to 0 (not useful)
        qa_bits(1) = 0

c ... If there were still some bands that were useful,
c ... then scale the qaulity accordingly
c ... (< 3 set quality bit to 4)
      else if (nmtests .lt. 3) then
        call set_bit(testbits,0)  
        call set_qa_bit(qa_bits,0)
        call set_qa_bit(qa_bits,3)

c ... (< 7 set quality bit to 6)
      else if (nmtests .lt. 7) then
        call set_bit(testbits,0)  
        call set_qa_bit(qa_bits,0)
        call set_qa_bit(qa_bits,2)
        call set_qa_bit(qa_bits,3)

c     Else set qaulity to highest value of 7
      else
         call set_bit(testbits,0)
         do i = 0 , 3
           call set_qa_bit(qa_bits,i)
         enddo
      endif

c ... Now if area is in difficult processing path region then
c ...  reduce quality to 6
      if (snglnt) then
         if (qa_bits(1) .eq. 15) qa_bits(1) = 13
      endif

c ... save bit flags for the current element in the line array
      do 200 i = 1 , 10
        if (i .le. 6) bitarray(nc,i) = testbits(i)
        qa_bitarray(nc,i) = qa_bits(i)
  200 continue

      return
      end
