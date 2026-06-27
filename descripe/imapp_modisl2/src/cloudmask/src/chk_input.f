      subroutine chk_input(max_pixel,nscans,lines_in_edge,no_250,
     +                     Begin_G_Date,Begin_G_Time,End_G_Date,
     +                     End_G_Time,beg_lin,nlins,beg_ele,
     +                     neles,pcf_satid,maxele,maxlin,max_count,
     +                     max5km_lin,max5km_ele,TAI_start,TAI_end)

      implicit none
      save

c ... Include files Needed for time conversion
      INCLUDE 'PGS_PC.f'
      include 'global.inc'
      include 'mod35.inc'

c     Common Block for debugging code
      common / bug / debug, h_output

c-------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C     Routine for checking consistency of input parameters.
C
c!Input Parameters:
C max_pixel     Maximum number of pixels per scan in granule
c nscans        Number of scans for this granule
c lines_in_edge Number of lines outside of processing region
c no_250        Logical flag true if we are to use 250m data
c Begin_G_Date  Beginning date of granule
c Begin_G_Time  Beginning time of granule
c End_G_Date    Ending date of granule
c End_G_Time    Ending time of granule
c beg_lin    Beginning line number to process
c nlins      Number of lines to process
c beg_ele    Beginning element number to process
c neles      Number of elements to process
c pcf_satid  Satellite identifier
c
c
c!Output Parameters:
c maxele        Maximum number of elements possible in granule
c maxlin        Number of lines in this granule
c max_count     Number of lines + lines_in_edge of granule
c max5km_lin    Max number of lines in 5km array
c max5km_ele    Max number of elements in 5km array
c TAI_start     Beginning TAI time of granule
c TAI_end       Ending TAI time of granule
c
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c Externals:  error messaging routine   -  message.f
c
c!END
c-------------------------------------------------------------------

c     scalar arguments
      character*70 Begin_G_Time,End_G_Time,Begin_G_Date,End_G_Date
      character*1 store3,store6
      character*4 pcf_satid
      double precision TAI_start,TAI_end
      integer maxlin,maxele,max_count,nscans,max_pixel,lines_in_edge,
     +        max5km_lin,max5km_ele,beg_lin,nlins,beg_ele,neles
      logical no_250

c     local scalars
      character*(PGSd_PC_VALUE_LENGTH_MAX) btime_string,etime_string,
     *           store1,store2,store4,store5,store7,store8,store9,
     *           store10
      integer debug,h_output,LUN_Start,LUN_Stop,rtn,chk_input_L2,
     *        version,itest,ibls,ibes
      double precision Start1,Start2,Stop1,Stop2

c ... external functions
      integer pgs_pc_getconfigdata,pgs_td_utctotai
      external pgs_pc_getconfigdata,pgs_td_utctotai

c     external subroutines
      external message

c     Satellite identifier from L1b file.
      include 'platform_name.inc'

      integer num_args
      integer FlagRA
      character FlagBuff*10

      num_args = command_argument_count()
      
      if(num_args == 1) then
         call get_command_argument(1,FlagBuff)
         read (FlagBuff,*) FlagRA
      else
      !This is the default value
         FlagRA = 0
      endif


c ... initialize variables
      itest = -1
      btime_string = '    '
      etime_string = '    '

      start1 = 0.0d0
      start2 = 0.0d0
      stop1 = 0.0d0
      stop2 = 0.0d0

c     Check satellite platform names for consistency.

      if(pcf_satid .eq. 'AM1M' .or. pcf_satid .eq. 'PM1M') then
c        Compare to satellite platform from L1b file.
         if(pcf_satid .eq. 'AM1M') then
           if(platform_name(1:5) .ne. 'Terra' .and.
     +        platform_name(1:5) .ne. 'terra' .and.
     +        platform_name(1:5) .ne. 'TERRA') then
              call message( 'chk_input',
     +        'Error: Inconsistent Satellite identifiers from L1b and PCF '
     +        //'[OPERATOR ACTION: Contact SDST]',0,2 )
           end if
         else if(pcf_satid .eq. 'PM1M') then
           if(platform_name(1:4) .ne. 'Aqua' .and.
     +        platform_name(1:4) .ne. 'aqua' .and.
     +        platform_name(1:4) .ne. 'AQUA') then
              call message( 'chk_input',
     +        'Error: Inconsistent Satellite identifiers from L1b and PCF '
     +        //'[OPERATOR ACTION: Contact SDST]',0,2 )
           end if
         end if
      else
        call message( 'chk_input',
     +  'Error: Invalid satellite instrument name from PCF file '
     +  // pcf_satid // '[OPERATOR ACTION: Contact SDST]',0,2 )
      end if

c ... maxlin is number of lines to process
c ..  max_count is number of lines + edge lines
c ... maxele is taken right from granule metadata
c ... 5 km variables are for reduced resolution swath files
      maxlin = nscans * scans_cube
      max5km_lin = maxlin / 5
      max5km_ele = max_pixel / 5
      if (beg_lin .eq. 0 .and. beg_ele .eq. 0) then
          beg_lin = 1
          beg_ele = 1
          max_count = maxlin + lines_in_edge
          maxele = max_pixel
      else
          ibls = beg_lin
          if(ibls .eq. 1) then
             beg_lin = 1
          else
             beg_lin = ibls - ((nlcntx-1) / 2)
          end if
          ibes = beg_ele
          if(ibes .eq. 1) then
             beg_ele = 1
          else
             beg_ele = ibes - ((necntx-1) / 2)
          end if
          max_count = beg_lin + ((nlcntx-1) + (nlins - 1))
          maxele = beg_ele + ((necntx-1) + (neles - 1))
          if (max_count .gt. maxlin .or. maxele .gt. max_pixel) then
          call message( 'chk_input',
     &    'Error: Number of lines or elements to process was .gt. max'
     &    // char(10) // '[OPERATOR ACTION:  Suspect system problem.    If a fault is '
     &    // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.]',
     &    0, 2 )
          endif
      endif

      if(maxlin .le. 0 .or. maxele .le. 0 .or. beg_lin .le. 0 .or.
     &   beg_ele .le. 0) then
        call message( 'chk_input',
     &  'Error: Number of lines or elements to process was .le. 0.'
     &  // char(10) // '[OPERATOR ACTION:  Suspect system problem. If a fault is '
     &  // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.]',
     &  0, 2 )
      endif

c ....................................................................
      if(debug .gt. 0) then
         write(h_output,fmt=
     + '(1x,''Debug  lines  '',i10, '' and elements '',i10,/)')
     +   maxlin, maxele
      endif
c ....................................................................

c ... Added the time input check
c ... Read the processing time from the pcf file.
c ... Must be in these LUN'S
      LUN_Start = 10258
      LUN_Stop = 10259

c ... Get the beginning processing time - in form of a string
      rtn = pgs_pc_getconfigdata(LUN_Start,btime_string)
      if (rtn .ne. 0) then
        call message( 'chk_input',
     &  'Error extracting process start time from pcf LUN 10258.'
     &  // char(10) // '[Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     &  // char(10) // 'entry.  If LUN is nonexistent or RP syntax is incorrect, stage '
     &  // char(10) // 'correct PCF and rerun PGE.  Otherwise, notify SDST.]',
     &  0, 2 )
      endif
c ... Get the end processing time - in form of a string
      rtn = pgs_pc_getconfigdata(LUN_Stop,etime_string)
      if (rtn .ne. 0) then
        call message( 'chk_input',
     &  'Error extracting process end time from pcf LUN 10259.'
     &  // char(10) // '[Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     &  // char(10) // 'entry.  If LUN is nonexistent or RP syntax is incorrect, stage '
     &  // char(10) // 'correct PCF and rerun PGE.  Otherwise, notify SDST.]',
     &  0, 2 )
      endif

c ....................................................................
      if (debug .gt. 0) then
        write(h_output,'(10x,'' Temporal Processing time from PCF '')')
        write(h_output,'(10x,'' Beginning Date/Time          Ending Date
     +/Time'',/,5x,A30,5x,A30/)') btime_string,etime_string
      endif
c ....................................................................

c ... Now add a T to the date stamp extracted from the geolocation
c ... file
      store3 = 'T'
      call concatenate_2(Begin_G_Date,store3,store1)
      call concatenate_2(End_G_Date,store3,store2)

c ... Now concatenate the T added date with the time
      call concatenate_2(store1,Begin_G_Time,store4)
      call concatenate_2(store2,End_G_Time,store5)

c ... Both the granule extracted and processing time/dates needs
c ... to have a Z added to the end for the PGS conversion tool
c ... to work
      store6 = 'Z'
      call concatenate_2(store4,store6,store7)
      call concatenate_2(store5,store6,store8)
      call concatenate_2(btime_string,store6,store9)
      call concatenate_2(etime_string,store6,store10)

c     Now convert our CCSDS ASCII time in Code A format (got that?)
c     and convert both to TAI time for comparison
      rtn = pgs_td_utctotai(store7,start1)
      if (rtn .ne. 0) then
        call message( 'chk_input',
     &  'Error converting granule start time to TAI.'
     &  // char(10) // ' [OPERATOR ACTION: Notify SDST]',
     &  0, 2 )
      endif
      rtn = pgs_td_utctotai(store8,stop1)
      if (rtn .ne. 0) then
        call message( 'chk_input',
     &  'Error converting granule stop time to TAI.'
     &  // char(10) // ' [OPERATOR ACTION: Notify SDST]',
     &  0, 2 )
      endif
      rtn = pgs_td_utctotai(store9,start2)
      if (rtn .ne. 0) then
        call message( 'chk_input',
     &  'Error converting processing start time to TAI.'
     &  // char(10) // ' [OPERATOR ACTION: Notify SDST]',
     &  0, 2 )
      endif
      rtn = pgs_td_utctotai(store10,stop2)
      if (rtn .ne. 0) then
        call message( 'chk_input',
     &  'Error converting processing stop time to TAI.'
     &  // char(10) // ' [OPERATOR ACTION: Notify SDST]',
     &  0, 2 )
      endif

c ....................................................................
      if (debug .gt. 0) then
        write(h_output,'(10x,'' Converted TAI time from granule '')')
        write(h_output,'(10x,'' Beginning Date/Time          Ending Date
     +/Time'',/,5x,d25.10,5x,d25.10/)') start1,stop1
        write(h_output,'(10x,'' Converted TAI time from PCF file '')')
        write(h_output,'(10x,'' Beginning Date/Time          Ending Date
     +/Time'',/,5x,d25.10,5x,d25.10/)') start2,stop2
      endif
c ....................................................................

c     Place TAI times from granule into variables to be passed out
      TAI_start = start1
      TAI_end = stop1

c ... Now compare the actual granule time with the interval that
c ... has been processed to make sure we are working on the
c ... the granule that we asked for.

      if (start1 .lt. start2  .or.  start1 .gt. stop2) then
c ....................................................................
         if (debug .gt. 0) then
          write(h_output,'(5x,'' Granule time does not fit into processi
     +ng time '',/3D25.10,/)') start1, start2, stop2
         endif
c ....................................................................

        call message( 'chk_input',
     &  'Error: Granule time does not match processing interval.'
     &  // char(10) // ' [OPERATOR ACTION: Load correct Geolocation file]',
     &  0, 2 )

      endif

c ... Added Riches test of other input granules

      version = 1
      if( FlagRA == 1) then
         itest = chk_input_L2(LRN_DSL1B_1km_RA,version)
      else
         itest = chk_input_L2(LRN_DSL1B_1km,version)
      endif
      if (itest .ne. 0) then
        call message( 'chk_input',
     &  'Granule L1B 1km time does not match processing interval.'
     &  // char(10) // ' [OPERATOR ACTION: Load correct L1B file]',
     &  0, 2 )
      endif
c ... Check 250 only if it was successfully opened
      if (.not. no_250) then
        itest = chk_input_L2(LRN_L1B_250,version)
        if (itest .ne. 0) then
          call message( 'chk_input',
     &    'Granule L1B 250m time does not match processing interval:' //
     &    ' 250m L1B data will not be used.' // char(10) //
     &    '[OPERATOR ACTION: In future, stage correct 250m L1B granule.' //
     &    ' If error persists, contact SDST]', 0, 0 )
          no_250 = .true.
        endif
      endif
      itest = chk_input_L2(LRN_Geo,version)
      if (itest .ne. 0) then
        call message( 'chk_input',
     &  'Granule Geo. time does not match processing interval.'
     &  // char(10) // ' [OPERATOR ACTION: Load correct Geo. file]',
     &  0, 2 )
      endif

      return
      end
