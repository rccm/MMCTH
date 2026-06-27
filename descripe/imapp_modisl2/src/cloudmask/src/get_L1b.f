      subroutine get_L1b(modfil_L1B_1km,fname_1km,modfil_L1B_250m,
     +                   fname_250m,no_250,cube,buf_size,scan_type,
     +                   v250_band,v1km_band,vis250_band,rad_band)

      implicit none
      save

c ... Common statement for debug purposes
      common / bug / debug, h_output

c----------------------------------------------------------------------
c!F77
c
c!Description:
c     Routine for reading MODIS level 1B radiances and relfectances
c     and placing in appropriate arrays for processing.
C
c!Input Parameters:
c modfil_L1B_1km Level 1B hdf file array
c fname_1km     Level 1B hdf file name
c modfil_L1B_250m Level 1B hdf file array
c fname_250m    Level 1B hdf file name
c no_250        Logical flag - true if 250 m file valid
c cube          Scan cube counter
c buf_size      Line and element granule size
c scan_type     Day or night scan
c
c!Output Parameters:
c v250_band     Byte array containing valid flag for a scan of 250m
c               pixels
c v1km_band     Byte array containing valid flag for a scan of 1km
c               pixels
c vis250_band   Array containing the 2 250m resolution visible
c               channel reflectances for a scan cube
c rad_band      Array containing the km resolution channel
c               radiances and brightness temperations for a scan cube
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c----------------------------------------------------------------------

        include 'mapi.inc'
        include 'PGS_PC.f'
        include 'global.inc'

c       scalar arguments
        character*1   scan_type
        integer cube
        logical no_250

c       argument arrays
        character*(PGSd_PC_FILE_PATH_MAX) fname_1km, fname_250m
        byte v250_band(nx,ny,vis_band),v1km_band(nx,ny,inband)
        integer buf_size(2),modfil_L1B_1km(MODFILLEN),
     +          modfil_L1B_250m(MODFILLEN)
        real rad_band(nx,ny,inband),vis250_band(nx,ny,vis_band)

c       local scalars
        character*6 type
        character*6 band_no
        integer ll,mm,nn,nch,ii,jj,ElementSize,LineSize,debug,h_output,
     +          i,j,k,res
        logical rtn

c       local arrays
        byte v_flag(nx,ny),v1km_flag(nx,ny)
        byte buf_unc(nx,ny),buf_sam(nx,ny)
        integer data_size(2),kk,V_code(nx,ny,inband)
        real buf(nx,ny)

c       external subroutines
        external Read_L1B,message

c ...   debug statement ............................................
        if (debug .gt. 0) then
          write(h_output,'(10x/,''Subroutine get_L1b  '',/)')
        endif
c ..........................................................

c ...   Set array sizes needed to read data
        ElementSize = nx
        LineSize = ny

c ...   Initialize main radiance array holders
         do 101 ll = 1 , inband
           do 102 mm = 1 , buf_size(2)
             do 103 nn = 1 , buf_size(1)
               if (ll .le. 2) then
c                250 m visible channel array holder
                 vis250_band(nn,mm,ll) = -99.0
                 v250_band(nn,mm,ll) = -1
               end if
c              1km resolution array holder
               rad_band(nn,mm,ll) =  -99.0
               v1km_band(nn,mm,ll) = -1
               V_code(nn,mm,ll) = -1
  103        continue
  102      continue
  101    continue

c ...    Now extract the pixels if instrument is in Day mode

c ...    debug statement ............................................
         if (debug .gt. 3) then
          write(h_output,'(10x/,''Scan mode is: '',A1,/)') scan_type
         endif
c ..........................................................


c----------------------------------------------------------------------
         if (scan_type .eq. 'D') then

           Do 110 nch = 1 , inband

c ...        Initialize HDF call ouput variables
             data_size(1) = 0
             data_size(2) = 0
             rtn = .false.
             res = 1

c ...        get solar bands in units of relfectance...thermal
c ...        channels in terms of radiance
             if (nch .le. inband-ir_band .or. nch .eq. 26) then
                type = 'refl'
             else
                type = 'rad'
             end if

c .......... debug statement .......................................
             if (debug .gt. 3) then
               write(h_output,'(10x,''channel, cal type: '',i5,A7)')
     +         nch, type
             endif
c ..........................................................

c ...        Read a 1km scan cube-s worth of data out of the granule
             call Read_L1B(modfil_L1B_1km,fname_1km,cube,nch,'H',
     +          res,type,ElementSize,LineSize,buf,buf_unc,
     +          buf_sam,v1km_flag,V_code,data_size,rtn)

             if (rtn) then
                write(band_no,'(I6)') nch
                call message('get_L1b',
     +           'Failed to extract 1km L1B data for band '// band_no //
     +           ' [OPERATOR ACTION: Check input L1B, rerun PGE]',
     +            0, 2 )
             end if

c ...        Place the data into the 1km holding arrays
             do 120 ii = 1 , data_size(2)
               do 130 jj = 1 , data_size(1)
                  rad_band(jj,ii,nch) = buf(jj,ii)
                  v1km_band(jj,ii,nch) = v1km_flag(jj,ii)
 130           continue
 120         continue

 110       continue
c----------------------------------------------------------------------

c          Under certain conditions, change the earth view data validity
c          flag so that the data may be used in the cloud mask.

           do 300 ii = 1, data_size(2)
             do 310 jj = 1, data_size(1)

c              Check for saturated pixels in band 2. Assume saturation
c              if error code is 65528 and channel 1 contains valid data.

               if(v1km_band(jj,ii,2) .eq. -1) then
                 if(V_code(jj,ii,2) .eq. 65528 .and. v1km_band(jj,ii,1)
     +           .eq. 0) then
                   v1km_band(jj,ii,2) = 0
                 end if
               end if

 310         continue
 300       continue

c----------------------------------------------------------------------
c ...      Now for the 250 m pixels - only extract data if 250 m file
c ...        is valid - check no_250 logical flag
           if (.not. no_250) then
             Do 140 nch = 1 , 2

c ...          Initialize HDF call ouput variables
               data_size(1) = 0
               data_size(2) = 0
               rtn = .false.

c ...          get 250m solar bands in units of relfectance
               type = 'refl'
               res = 16

c ..........   debug statement .......................................
               if (debug .gt. 3) then
                 write(h_output,'(10x,''Band, cal, size:  '',i5,A7,2i5/)')
     +           nch, type, data_size(1), data_size(2)
               endif
c ..........................................................

c ...          Read a 250m scan cube-s worth of data out of the granule
               call Read_L1B(modfil_L1B_250m,fname_250m,cube,nch,
     +            ' ',res,type,ElementSize,LineSize,buf,buf_unc,
     +              buf_sam,v_flag,V_code,data_size,rtn)

               if (rtn) then
                  write(band_no,'(I6)') nch
                  call message('get_L1b',
     +          'Failed to extract 250m L1B data for band '// band_no //
     +          ' [OPERATOR ACTION: Check input L1B, rerun PGE]',
     +              0, 0 )
               end if

c ...          Place the data into the 250m  holding arrays
c ...          Make compatible with v1 array indexing
               do 150 ii = 1 , data_size(2) , 4
                 do 160 jj = 1 , data_size(1) , 4
                    do 170 kk = 0 , 3
                       do 180 ll = 0 , 3
                          vis250_band(jj+ll,ii+kk,nch) =
     +                            buf(jj+ll,ii+kk)
                          v250_band(jj+ll,ii+kk,nch) =
     +                       v_flag(jj+ll,ii+kk)
 180                   continue
 170                continue
 160              continue
 150            continue


 140         continue
c----------------------------------------------------------------------

c            Under certain conditions, change the earth view data validity
c            flag so that the data may be used in the cloud mask.

             do 320 ii = 1, data_size(2)
               do 330 jj = 1, data_size(1)

c                Check for saturated pixels in band 2. Assume saturation
c                if error code is 65533 and channel 1 contains valid data.

                 if(v250_band(jj,ii,2) .eq. -1) then
                   if(V_code(jj,ii,2) .eq. 65533 .and. v250_band(jj,ii,1)
     +             .eq. 0) then
                     v250_band(jj,ii,2) = 0
                   end if
                 end if

 330           continue
 320         continue

c----------------------------------------------------------------------
c ...   debug statement ............................................
        if (debug .gt. 3) then
          write(h_output,'(15x,'' CONTEXT FOR 250 CUBE '', 3I10)') cube,
     +                   data_size(1), data_size(2)
          write(h_output,'(2x,'' vis250_band 1(.66um)'',/,40(10f10.6/))')
     +              (((vis250_band(i,j,k),i=2780,2789),j=1,40),k=1,1)
          write(h_output,'(2x,'' vis250_band 2(.87um)'',/,40(10f10.6/))')
     +              (((vis250_band(i,j,k),i=2780,2789),j=1,40),k=2,2)
          write(h_output,'(2x,'' V_code 1(.66um)'',/,40(10i10/))')
     +              (((V_code(i,j,k),i=2780,2789),j=1,40),k=1,1)
          write(h_output,'(2x,'' V_code 2(.87um)'',/,40(10i10/))')
     +              (((V_code(i,j,k),i=2780,2789),j=1,40),k=2,2)
          write(h_output,'(2x,'' v250_band 1(.66um)'',/,40(10i10/))')
     +              (((v250_band(i,j,k),i=2780,2789),j=1,40),k=1,1)
          write(h_output,'(2x,'' v250_band 2(.87um)'',/,40(10i10/))')
     +              (((v250_band(i,j,k),i=2780,2789),j=1,40),k=2,2)
        endif
c ................................................................


         endif


c----------------------------------------------------------------------
         else
c ...      NOW proceed with the night only scans

           Do 190 nch = inband - ir_band + 1 , inband

c ...        Do not extract data if band is 26 (1.38 micron)
             if (nch .ne. 26) then

c ...          Initialize HDF call ouput variables
               data_size(1) = 0
               data_size(2) = 0
               rtn = .false.

c ...          get thermal channelsin terms of radiance
               type = 'rad'
               res = 1

c ..........   debug statement .......................................
               if (debug .gt. 3) then
                 write(h_output,'(10x,''channel, cal type: '',i5,A7)')
     +           nch, type
               endif
c ..........................................................

c ...           Read a 1km scan cube-s worth of data out of the granule
                call Read_L1B(modfil_L1B_1km,fname_1km,cube,nch,' ',
     +               res,type,ElementSize,LineSize,buf,buf_unc,
     +               buf_sam,v1km_flag,V_code,data_size,rtn)

                if (rtn) then
                  write(band_no,'(I6)') nch
                  call message('get_L1b',
     +             'Failed to extract 1km L1B data for band '
     +             // band_no // ' [OPERATOR ACTION: Check input
     +             L1B, rerun PGE]',0, 2 )
                end if

c ...           Place the data into the 1km holding arrays
                do 200 ii = 1 , data_size(2)
                   do 210 jj = 1 , data_size(1)
                    rad_band(jj,ii,nch) = buf(jj,ii)
                      v1km_band(jj,ii,nch) = v1km_flag(jj,ii)
 210               continue
 200            continue
              endif

 190     continue

         endif
c----------------------------------------------------------------------


c ... debug statement ............................................
      if (debug .gt. 3) then
        write(h_output,'(15x,'' CONTEXT FOR 1km CUBE '', 3I10)') cube,
     +                   data_size(1), data_size(2)
        write(h_output,'(15x,'' Pixels 696 - 705'')')
        write(h_output,'(2x,'' rad_band 20 (3.7um)'',/,20(5f10.4/))')
     +       (((rad_band(i,j,k),i=696,705),j=1,10),k=20,20)
        write(h_output,'(2x,'' rad_band 31 (11um)'',/,20(5f10.4/))')
     +       (((rad_band(i,j,k),i=696,705),j=1,10),k=31,31)
        write(h_output,'(2x,'' rad_band 1 (.66um)'',/,20(5f10.4/))')
     +       (((rad_band(i,j,k),i=696,705),j=1,10),k=1,1)
        write(h_output,'(2x,'' rad_band 2 (.87um)'',/,20(5f10.4/))')
     +       (((rad_band(i,j,k),i=696,705),j=1,10),k=2,2)
        write(h_output,'(2x,'' V_code 1 (.66um)'',/,20(5i10/))')
     +       (((V_code(i,j,k),i=696,705),j=1,10),k=1,1)
        write(h_output,'(2x,'' V_code 2 (.87um)'',/,20(5i10/))')
     +       (((V_code(i,j,k),i=696,705),j=1,10),k=2,2)
        write(h_output,'(2x,'' v1km_band 1 (.66um)'',/,20(5i10/))')
     +       (((v1km_band(i,j,k),i=696,705),j=1,10),k=1,1)
        write(h_output,'(2x,'' v1km_band 2 (.87um)'',/,20(5i10/))')
     +       (((v1km_band(i,j,k),i=696,705),j=1,10),k=2,2)
      endif
c ................................................................

         return
         end
