      subroutine read_zonbias(bias_file,zbias_oc,zbias_dl,
     *                        zbias_nl,csr_success)

c-----------------------------------------------------------------------
c !F77
c     
c !DESCRIPTION:
c      Read zonal mean clear-sky radiance biases from staged HDF file.
c
c !INPUT PARAMETERS:
c      bias_file Name of file to be read
c    
c !OUTPUT PARAMETERS:
c      znbias_oc   Biases computed for ocean surfaces
c      znbias_dl   Biases computed for daytime land surfaces
c      znbias_nl   Biases computed for nighttime land surfaces
c      csr_success Indicates successful read of file
c    
c !REVISION HISTORY: 
c
c 09/26/05:        Initial version of code
c
c !TEAM-UNIQUE HEADER:
c      Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

c     Include files.
      include 'mapi.inc'
      include 'mod06uw_debug.inc'

      integer bias_bands, bias_zones
      parameter (bias_bands = 5)
      parameter (bias_zones = 181)
      
c     Argument variables.
      character*255 bias_file
      integer modnum
      logical csr_success
      real zbias_oc(bias_zones,bias_bands),zbias_dl(bias_zones,bias_bands),
     *     zbias_nl(bias_zones,bias_bands),csr_data(bias_zones,bias_bands)

c     Local variables.
      character*80 arrnm, grpnm, data_type, attrib
      integer level, status, rank, nstats, nzones, fvalue,
     *        i, j, k, jj, iprd

      integer nprods
      parameter (nprods = 3)

c     Arrays.
      character*80 sds(nprods) / 'Ocean_Clear_Sky_Bias',
     *                           'Land_Day_Clear_Sky_Bias',
     *                           'Land_Night_Clear_Sky_Bias'/
      integer modfil(MODFILLEN), dims(2), start(2), edge(2)

      integer nms

      csr_success = .true.

c     Open zonal bias file.

      status = opmfil( bias_file, 'r', modfil)

      if( status .ne. 0 ) then

        level = 1
        csr_success = .false.
        call message( 'read_zonbias ',
     *    'Error opening input bias file ' //
     *    '[OPERATOR ACTION: Check if file exists]', status, level )
        return

      endif

c     Read data from bias file.  Get mean clear-sky biases for each   
c     1-degree latitude zone and 5 MODIS bands (31, 33-36).

      do iprd = 1, nprods

        arrnm = sds(iprd)
        grpnm=' '
        data_type = ' '
        rank = 2
        if(GMARDM(modfil, arrnm, grpnm, data_type, rank, dims)
     *            .ne. MAPIOK) then
          level = 1
          status = -1
          csr_success = .false.
          call message( 'read_zonbias ',
     *      'Failed to extract information for ' // arrnm //
     *      ' [OPERATOR ACTION: Contact SDST.]', status, level)
          go to 1000
        else
          nzones = dims(1)
          nstats = dims(2)
c         Extract fill value
          attrib = '_FillValue'
          nms = 1
          if (GMARIN (modfil,arrnm,grpnm,attrib,data_type,nms,fvalue)
     *                .ne. MAPIOK) then
            level = 1
            status = -1
            csr_success = .false.
            call message( 'read_zonbias ',
     *        'Failed to extract fill value for ' // arrnm //
     *        ' [OPERATOR ACTION: Contact SDST.]', status, level)
            go to 1000
          endif
        end if

c       Initialize array.
        do j = 1 , bias_bands
          do i = 1 , bias_zones
            csr_data(i,j) = fvalue*1.0
          enddo
        enddo

c       Read bias values.

c       Define start and end of HDF read.
        start(1) = 0
        start(2) = 0
        edge(1) = nzones
        edge(2) = nstats

c       Extract data
        if(GMAR(modfil,arrnm,grpnm,start,edge,csr_data)
     *          .ne. MAPIOK)then
          level = 1
          status = -1
          csr_success = .false.
          call message( 'read_zonbias ',
     *      'Failed to extract data from ' // arrnm //
     *      ' [OPERATOR ACTION: Contact SDST.]', status, level)
          go to 1000
        endif

        if(debug .gt. 0) then
          write(h_output,'(''Zonal clear-sky radiance bias data '',i5)') iprd
          do k = 1, 181
            write(h_output,'(i10,5f12.3)') k,(csr_data(k,jj),jj=1,5)
          enddo
        end if

        if(iprd .eq. 1) then
          do i = 1, bias_bands
            do j = 1, bias_zones
              zbias_oc(j,i) = csr_data(j,i)
            enddo
          enddo
        else if(iprd .eq. 2) then
          do i = 1, bias_bands
            do j = 1, bias_zones
              zbias_dl(j,i) = csr_data(j,i)
            enddo
          enddo
        else if(iprd .eq. 3) then
          do i = 1, bias_bands
            do j = 1, bias_zones
              zbias_nl(j,i) = csr_data(j,i)
            enddo
          enddo
        end if

      enddo

 1000 continue

c     Close bias file.
      status = clmfil(modfil)

      if( status .ne. 0 ) then

        level = 1
        csr_success = .false.
        call message( 'read_csr ',
     &    'Error closing input bias file. ' //
     &    '[OPERATOR ACTION: Check system resources]', status, level )

      endif

      return
      end
