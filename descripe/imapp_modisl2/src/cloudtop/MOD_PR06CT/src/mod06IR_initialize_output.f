      subroutine mod06IR_initialize_output( qual_flag_ir, conf_flag_ir,
     *  maxl, clus, gfit, Chisq, nmaxl, icode, indat_bt, var)
      
c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Initialize output product arrays for MOD06IR 5x5 sampling,
c    as well as QA flag information
c
c!Input Parameters:
c
c    QA flags:
c    qual_flag_ir                quality of product
c    conf_flag_ir                confidence in product
c    maxl                        flag for implemtation of maximum
c                                  likelihood estimator
c    clus                        flag for implemtation of cluster
c                                  analysis
c    gfit                        goodness of fit of MLE result
c    Chisq                       Chi_sqared result
c
c    Processing Variables:
c
c    indat_bt                    holds bt values for processing
c    icode                       holds 1D array of phase product
c    var                         holds variance values for processing
c    nmaxl                       Number of pixels used in max likelihood
c                                  extimator
c!Output Parameters:
c
c    QA flags:
c    qual_flag_ir                quality of product
c    conf_flag_ir                confidence in product
c    maxl                        flag for implemtation of maximum
c                                  likelihood estimator
c    clus                        flag for implemtation of cluster
c                                  analysis
c    gfit                        goodness of fit of MLE result
c    Chisq                       Chi_sqared result
c
c    Processing Variables:
c
c    indat_bt                    holds bt values for processing
c    icode                       holds 1D array of phase product
c    var                         holds variance values for processing
c    nmaxl                       Number of pixels used in max likelihood
c                                  extimator

c
c    The following arrays in COMMON /MOD06_DATA/ are initialized:
c    radiance_variance           Variance of radiances over box
c                                 for all three ir phase bands
c    cloud_phase                 Tri-spectral irphase result
c    product_qa                  Output product qa array
c    sfc_temp                    Ouput product for surface temperature
c    sfc_pres                    Ouput product for surface pressure
c
c!Revision History:
c $Id: mod06IR_initialize_output.f,v 1.5.2.1 2001/04/30 13:58:38 raf Exp $
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

      implicit none
      save
          
      include 'mod06uw_data.inc'
      
c ... arguments 
      real gfit, Chisq
      integer maxl, clus, nmaxl

c ... array arguments
      integer indat_bt( max_samp_line*max_samp_pixel , irphase_band+5 ),
     +        icode( max_samp_line*max_samp_pixel )
      real var ( max_samp_line*max_samp_pixel , irphase_band )
      byte qual_flag_ir ( max_samp_line*max_samp_pixel ),
     +     conf_flag_ir ( max_samp_line*max_samp_pixel )

c ... local variables

      integer i, j, k

c ... Initialize all the output arrays

c ... radiance variance for cloud ir phase
      do k = 1, n_output
        do j = 1, max_samp_line
          do i = 1, max_samp_pixel
            radiance_variance( i, j, k ) = bad_value
          end do
        end do
      end do

c ... output product qa holder
      do k = 1, nproduct_qa
        do j = 1, max_samp_line
          do i = 1, max_samp_pixel
            product_qa( k, i, j ) = 0
          end do
        end do
      end do

c ... Phase product and surface type, surface p,and temp
      do j = 1, max_samp_line
        do i = 1, max_samp_pixel
          cloud_phase( i, j ) = bad_value
          day_night_flag( i, j ) = -1
          sfc_temp( i, j ) = bad_value
          sfc_pres( i, j ) = bad_value
        end do
      end do                      

c ... Brightness temperature
      do k = 1, n_output
        do j = 1, max_samp_line
          do i = 1, max_samp_pixel
            brightness_temp( i, j, k ) = bad_value
          end do
        end do
      end do

c ... initialize processing arrays

      do j = 1, max_samp_line * max_samp_pixel
        icode( j ) = -1
        qual_flag_ir( j ) = 0
        conf_flag_ir( j ) = 0
        do i = 1, irphase_band+5
          if (i .le. irphase_band) var( j, i ) = bad_value
          indat_bt( j, i ) = out_misg
        end do
      end do

c     Initialize flags
      maxl = 0
      clus = 0
      gfit = 0.0
      Chisq = 0.0
      nmaxl = 0


      return
      end
