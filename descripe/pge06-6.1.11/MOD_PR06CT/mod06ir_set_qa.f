      subroutine mod06ir_set_qa( npixels,qual_flag_ir, conf_flag_ir,
     *  maxl, clus, gfit, Chisq, nmaxl, icode, ng_irp, 
     *  ni, nw, nx, nu ) 

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Set the appropriate MOD06 QA bits for a granule of
c    data - this one is for the IRPHASE portion of MOD06
c
c!Input Parameters:
c
c    QA flags:
c    npixels       Number of pixels in scan
c    qual_flag_ir  quality of product
c    conf_flag_ir  confidence in product
c    maxl          flag for implemtation of maximum
c                     likelihood estimator
c    clus          flag for implemtation of cluster
c                     analysis
c    gfit          goodness of fit of MLE result
c    Chisq         Chi_sqared result
c    nmaxl         Number of pixels used in max likelihood
c                     extimator
c    icode                       holds 1D array of phase product
c
c
c!Output Parameters:
c
c          granule based counters
c    ng_irp        sum of good irphase retrievals
c    ni            sum of ice clouds
c    nw            sum of water clouds
c    nx            sum of mixed phase clouds
c    nu            sum of uncertain phase clouds
c
c
c    The following arrays in COMMON /MOD06_DATA/ are initialized:
c    product_qa     Output product qa array 
c
c!Revision History:
c $Id: mod06ir_set_qa.f,v 1.6 1999/06/11 22:36:10 kis Exp $
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------
      
      implicit none
      save
      
      include 'mod06uw_data.inc'
      include 'mod06uw_debug.inc'

c ... scalar arguments
      real gfit, Chisq
      integer maxl, clus, nmaxl, ng_irp, ni, nw, nx, nu, npixels

c ... array arguments
      integer icode( max_samp_line*max_samp_pixel )
      byte qual_flag_ir ( max_samp_line*max_samp_pixel ),
     +     conf_flag_ir ( max_samp_line*max_samp_pixel )

c ... parameter
      real epsilon
      parameter (epsilon = 1.0e-10)

c ... local variables
      logical lval
      integer i, j, k, ii, jj, kk, kount, ipos, ival, ibyte
      integer qa(80)

c ... First setbits in single dimension holder and set meta counters
      do i = 1 ,  max_samp_line*max_samp_pixel
       
c ...   quality bits
        if (icode(i) .ne. -1 .and. icode(i) .ne. 100 .and. 
     +      icode(i) .ne. 0) then
           qual_flag_ir(i) = 1
c ...      confidence bits
           if (icode(i) .ne. 3 .and. icode(i) .ne. 6) then
              conf_flag_ir(i) = 1
           else
              conf_flag_ir(i) = 2
           endif
c ...      if good product add to successful. retrieval meta counter
           ng_irp = ng_irp + 1
c ...      if water cloud add to water cloud counter
           if (icode(i) .eq. 1 .or. icode(i) .eq. 5 )
     +          nw = nw + 1
c ...      if ice cloud add to ice cloud counter
           if (icode(i) .eq. 2 .or. icode(i) .eq. 4 )
     +          ni = ni + 1
c ...      if mixed phase cloud add to mixed cloud counter
           if (icode(i) .eq. 3) nx = nx + 1
c ...      if uncertain phase cloud add to uncertain cloud counter
           if (icode(i) .eq. 6) nu = nu + 1
        endif

c ...   debug statement ............................................
        if (debug .gt. 3) then
           write(h_output,'(10x,'' Metadata sums '',7i6/)')
     +         ng_irp, ni, nx, nw, nu, icode(i), maxl
        endif
c ................................................................
      end do

      kount = 0

c ... Now fill the output qa bit array
      do i = 1, ( max_line / isamp ) 
         do j = 1, ( npixels / isamp ) 
             kount = kount + 1
c ...        Set single fov quality flag bit
             if (qual_flag_ir(kount) .eq. 1) then
               call set_qa_bit(product_qa(1,j,i),16)

c ...          Set single fov confidence bits
               if (conf_flag_ir(kount) .eq. 1) then
                 call set_qa_bit(product_qa(1,j,i),17)
                 call set_qa_bit(product_qa(1,j,i),18)
               else if (conf_flag_ir(kount) .eq. 2) then
                 call set_qa_bit(product_qa(1,j,i),17)
               endif
  
c ...          Set bit 48 if maximum likelihood estimator was invoked
               if (maxl .eq. 1) 
     +           call set_qa_bit(product_qa(1,j,i),48)
  
c ...          Set bit 49 if cluster analysis was performed
               if (clus .eq. 1) 
     +           call set_qa_bit(product_qa(1,j,i),49)

c ...          Set bit 50 if fit was good
               if (abs(gfit-1.0) .lt. epsilon) 
     +           call set_qa_bit(product_qa(1,j,i),50)

c ...          Set bit 51 if Chisquare solution was good
               if (nint(Chisq) .lt. nmaxl) 
     +           call set_qa_bit(product_qa(1,j,i),51)
            endif

         end do
      end do

c ... debug statement ............................................
      if (debug .gt. 3) then
c        strip out and print bit values for current pixel

         kount = 0
         do i = 1, ( max_line / isamp ) 
           do j = 1, ( npixels / isamp )
              kount = kount + 1
              ipos = 0
              do ii = 1 , 80
                qa(ii) = 0
              enddo
              do ibyte = 1,10
                ival = product_qa(ibyte,j,i)
                do k = 1,8
                  ipos = ipos + 1
                  lval = btest(ival,k-1)
                  if(lval) then
                     qa(ipos) = 1
                  else
                     qa(ipos) = 0
                  end if
                enddo
              enddo

         write(h_output,'(1x,3x,'' QA Bits Results'')')
         write(h_output,'(1x,3x,'' Line Numuber '', i6)') i
         write(h_output,'(1x,3x,''FOV No'',3x,16(i2,2x))') (kk,kk=0,15)
         write(h_output,'(1x)')
         write(h_output,'(1x,i8,2x,16i4)') j,(qa(jj),jj=1,16)
         write(h_output,'(1x,i8,2x,16i4)') j,(qa(jj),jj=17,32)
         write(h_output,'(1x,i8,2x,16i4)') j,(qa(jj),jj=33,48)
         write(h_output,'(1x,i8,2x,16i4)') j,(qa(jj),jj=49,64)
         write(h_output,'(1x,i8,2x,16i4,/)') j,(qa(jj),jj=65,80)

             write(h_output,'(10x,'' Quality bits '',5i6,2f8.2/)')
     +        qual_flag_ir(kount),conf_flag_ir(kount),icode(kount),
     +        maxl,nmaxl,gfit,Chisq

             write(h_output,'(10x,'' qa_bits pixel results ''/10i5/)')
     +        product_qa(1,j,i),product_qa(2,j,i),product_qa(3,j,i),
     +        product_qa(4,j,i),product_qa(5,j,i),product_qa(6,j,i),
     +        product_qa(7,j,i),product_qa(8,j,i),product_qa(9,j,i),
     +        product_qa(10,j,i)
         end do
       end do

      endif
c ................................................................

      return
      end
