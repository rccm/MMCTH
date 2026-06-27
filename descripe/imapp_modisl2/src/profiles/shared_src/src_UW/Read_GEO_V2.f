       INTEGER FUNCTION READ_GEO_V2(MODFIL,ICODE,SCANCUBE_NO,
     & NELES,NLINES,BUF,DATA_SIZE)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:   Retrieves one scan cube of MODIS geolocation data HDF
C                target array of typically 100 scan cubes (a granule).
C
C!INPUT PARAMETERS:
C    INTEGER     modfil       Array containing SD_ID, File_ID and
C                             File access type.
C    INTEGER     icode        Code representing geolocation variable
c                             to extract from file.  Code represents:
C                             1 - Latitude (degrees)
C                             2 - Longitude (degress)
C                             3 - height (meters)
C                             4 - sensor zenith angle (degrees)
C                             5 - sensor azimuth angle (degrees)
C                             6 - range (meters)
C                             7 - solar zenith angle (degrees)
C                             8 - solar azimuth angle (degress)
C                             9 - land/sea mask (degress)
C
C    INTEGER     ScanCube_No  The ScanCube Number.
C    INTEGER     neles        Number of elements in scan
C    INTEGER     lines        Number of lines in scan
C
C!OUTPUT PARAMETERS:
C    REAL        Buf          Buffer containing requested scan cube
C                             radiance data.
C    INTEGER     Data_Size(2) Array showing the real size of SDS EV data.
C
C!REVISION HISTORY:
C
C!TEAM-UNIQUE HEADER:
C This software is developed by the MODIS Science Data Support
C Team for the National Aeronautics and Space Administration,
C Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C!REFERENCES AND CREDITS:
C
C       WRITTEN BY:
C       Xiao-Yang Ding                08/15/95
C       Research and Data systems Corporation
C       SAIC/GSC MODIS Science Data Support Office
C       7501 Forbes Blvd, Seabrook MD 20706
C
C       ding@ltpmail.gsfc.nasa.gov
C
C       MODIFIED BY :  Kathy Strabala
C
C!DESIGN NOTES:
C    This subroutine checks the return status of all MODIS Application
C    Program Interface (M-API) function calls.  A successful M-API
C    function call is indicated by a return value 0. If a M-API
C    function call is not successful, an error message is written
C    to the LogStatus file and the process aborts. This is achieved
C    by passing an fatal mnemonic error message (i.e. .._F_..)
C    to the function MODIS_SMF_SETDYNAMICMSG.
C
C    Externals:
C
C      Functions:
C        GMAR                       (libmapi.a)
C        GMARDM                     (libmapi.a)
C
C      Named Constants:
C        DFACC_READ                 (hdf.inc)
C        MAPIOK                     (mapi.inc)
C        MODFILLEN                  (mapi.inc)
C        MODIS_F_GENERIC            (MODIS_39500.f)
C
C    Internals:
C
C      Subroutines:
C        MODIS_SMF_SETDYNAMICMSG    (mod35_fileoc.f)
C        CONCATENATE                (mod35_fileoc.f)
C
C      Variables:
C        arrnam       Name of the SDS array.
C        grpnm        Name of the data group containing the target
C        data_type    String describing the data type of the array.
C        msgbuf*      Message buffer.
C        Dim_Size(2)  Array specifying the size of the SDS data.
C        Edge(2)      Array specifying the number of data value to read.
C        Start(2)     Array specifying the starting location of data.
C        Fmax         Maximum frame number per scan line.
C        Lmax         Maximum line number per scan cube.
C        Rank         The number of dimensions in an array
C        Read_Geo     The function return value
C        j,k,L,temp   Temporary variable for holding integer value
C        TSCN         Total ScanCube Number.
C        angles(*),elev(*),count(*)
C                     A temporary buffer for data of the target array.
C        NofLine_PerCube The number of lines per scan cube
C        scale        Scale factor
C
C!END
C------------------------------------------------------------------------

       IMPLICIT NONE

       INCLUDE 'mapi.inc'
       INCLUDE 'hdf.inc'
       INCLUDE 'PGS_MODIS_39500.f'

C Declarations
       CHARACTER*80 arrnm,grpnm,data_type
       CHARACTER*80 msgbuf,msgbuf1,msgbuf2
       INTEGER ScanCube_No,NofLine_PerCube,Rank,j,k,L,
     &         Fmax,Lmax,TSCN,icode
       PARAMETER(Fmax=1500,Lmax=10)

       Real count(Fmax*Lmax)
       Integer*2  angles(Fmax*Lmax),elev(Fmax*Lmax),dist(Fmax*Lmax)
       INTEGER Start(2),Edge(2),Dim_Size(2),neles,nlines,Data_Size(*),
     &         modfil(MODFILLEN)
       REAL Buf(Fmax,Lmax),scale
       BYTE ls(Fmax*Lmax)

C Initialization

       grpnm=' '
       data_type = ' '
       Read_GEO_V2=-1
       NofLine_PerCube=10
       Rank = 2
C Checking for valid file and band numbers
       if(modfil(1).le.0.or.modfil(3).ne.DFACC_READ)then
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Invalid SD_ID or invalid file access type'//
     &   ' [OPERATOR ACTION: Run PGE again on non-currupt Geolocation' //
     &   ' File. If file is okay, contact SDST]','Read_Geo_V2')
       endif

       if(icode.eq.1)then
         arrnm='Latitude'
         scale = 1.0
       elseif(icode.eq.2)then
         arrnm='Longitude'
         scale = 1.0
       elseif(icode.eq.3)then
         arrnm='Height'
         scale=1.0
       elseif(icode.eq.4)then
         arrnm='SensorZenith'
         scale=0.01
       elseif(icode.eq.5)then
         arrnm='SensorAzimuth'
         scale=0.01
       elseif(icode.eq.6)then
         arrnm='Range'
         scale=25.
       elseif(icode.eq.7)then
         arrnm='SolarZenith'
         scale=0.01
       elseif(icode.eq.8)then
         arrnm='SolarAzimuth'
         scale=0.01
       else
         arrnm='Land/SeaMask'
         scale = 1.0
       endif

C Retrieving the rank, dimensions and data type of SDS data.
      if(GMARDM(modfil, arrnm, grpnm, data_type, Rank, Dim_Size)
     & .ne.MAPIOK) then
       CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     & 'GMARDM get dim_info failed, [OPERATOR ACTION: '//
     & 'Notify SDST]','Read_Geo_V2')
      endif

C  Additional input check of ScanCube_No.
       TSCN=Dim_Size(2)/NofLine_PerCube
       if(ScanCube_No.lt.1.or.ScanCube_No.gt.TSCN) then
         write(msgbuf,'(i4)') TSCN
         call Concatenate('ScanCube_No out of bounds; range 1 -',msgbuf,
     &   msgbuf1)
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   msgbuf1 // ' [OPERATOR ACTION: Check input Geolocation file]',
     &   'Read_Geo_V2')
       endif
C Check if Buf_Size is large enough.
       Data_Size(1) = Dim_Size(1)
       Data_Size(2) = NofLine_PerCube
       if(neles.lt.Dim_Size(1).or.
     &   nlines.lt.NofLine_PerCube) then
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'Data_Size array not large enough ','Read_Geo_V2')
         write(msgbuf,'(2i5)')neles,nlines
         call Concatenate('array size:',msgbuf,msgbuf1)
         write(msgbuf,'(2i5)')Dim_Size(1),NofLine_PerCube
         call Concatenate('...  Data_Size:',msgbuf,msgbuf2)
         call Concatenate(msgbuf1,msgbuf2,msgbuf)
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,msgbuf //
     &   ' [OPERATOR ACTION: Notify SDST]', 'Read_Geo_V2')
       endif

C Initialize the Buffer
       do j=1,nlines
        do k=1,neles
          Buf(k,j)=0.0
        enddo
       enddo

c Initialize the parameters needed from the file
      do j = 1 , Fmax*Lmax
        count(j) = 0.0
        angles(j) = 0
        elev(j) = 0
        dist(j) = 0
        ls(j) = 0
      enddo

C Get Geo data
      Start(1)=0
      Start(2)=(ScanCube_No-1)*NofLine_PerCube
      Edge(1)=Dim_Size(1)
      Edge(2)=NofLine_PerCube

c     Loop for getting Zenith and Azimuth angles
      if((icode.eq.4).or.(icode.eq.5).or.(icode.eq.7).or.
     &  (icode.eq.8))then
        if(GMAR(modfil,arrnm,grpnm,Start,Edge,angles).ne.MAPIOK)then
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &    'GMAR geolocation read for zenith and azumuth angles failed.'//
     &    ' [OPERATOR ACTION: Run PGE again on non-currupt Geolocation'//
     &    ' File. If file is okay, contact SDST]','Read_Geo_V2')
        endif
        L=0
        do k=1,Edge(2)
         do j=1,Edge(1)
          L=L+1
          Buf(j,k)= real(angles(L))*scale
         end do
        end do

c     Loop for getting Latitude/Longitude values
      elseif((icode.eq.1).or.(icode.eq.2))then
        if(GMAR(modfil, arrnm, grpnm, Start, Edge, count).ne.MAPIOK)then
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &    'GMAR geolocation read for latitude/longitude values failed.'//
     &    ' [OPERATOR ACTION: Run PGE again on non-currupt Geolocation'//
     &    ' File. If file is okay, contact SDST]','Read_Geo_V2')
        endif
        L=0
        do k=1,Edge(2)
         do j=1,Edge(1)
          L=L+1
          Buf(j,k)=count(L)*scale
         end do
        end do

c     Loop for getting the Height
      elseif(icode.eq.3)then
        if(GMAR(modfil, arrnm, grpnm, Start, Edge, elev).ne.MAPIOK)then
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &    'GMAR geolocation read for height failed.' //
     &    ' [OPERATOR ACTION: Run PGE again on non-currupt Geolocation'//
     &    ' File. If file is okay, contact SDST]','Read_Geo_V2')
        endif
        L=0
        do k=1,Edge(2)
          do j=1,Edge(1)
            L=L+1
            Buf(j,k)= real(elev(L))*scale
          end do
        end do

c     Loop for getting the Range
      elseif(icode.eq.6)then
        if(GMAR(modfil, arrnm, grpnm, Start, Edge, dist).ne.MAPIOK)then
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &    'GMAR geolocation read for range failed.' //
     &    ' [OPERATOR ACTION: Run PGE again on non-currupt Geolocation'//
     &    ' File. If file is okay, contact SDST]','Read_Geo_V2')
        endif
        L=0
        do k=1,Edge(2)
          do j=1,Edge(1)
            L=L+1
            Buf(j,k)= real(dist(L))*scale
          end do
        end do

c     Must be asking for the land/sea mask
      else
        if(GMAR(modfil, arrnm, grpnm, Start, Edge, ls).ne.MAPIOK)then
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &    'GMAR geolocation read for land/sea mask failed.' //
     &    ' [OPERATOR ACTION: Run PGE again on non-currupt Geolocation'//
     &    ' File. If file is okay, contact SDST]','Read_Geo_V2')
        endif
        L=0
        do k=1,Edge(2)
          do j=1,Edge(1)
            L=L+1
            Buf(j,k)= real(ls(L))
          end do
        end do
      endif
      Read_GEO_V2=0

      END
