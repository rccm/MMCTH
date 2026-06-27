      integer function compare_times(dataset_time,tai_time,delta)

c-----------------------------------------------------------------------
c!F77
c
c!DESCRIPTION:
c    Using the output date/time from the get_ancillary subroutine,
c    Compare the time to the granule tai time.  Return success if
c    the time is within DELTA)
c
c!INPUT PARAMETERS:
c
c    dataset_time      Data set time as returned from the get_ancillary
c                        subroutine (array:year,month,day,hour) INTEGER
c    tai_time          TAI time of granule (REAL*8)
c    delta             Window of time that you will accept to use the
c                        ancillary data set in seconds (REAL*8)
c
c!OUTPUT PARAMETERS:
c    compare_times      0 => Success
c                      -1 => Failure
c
c!REVISION HISTORY:
c
c!TEAM-UNIQUE HEADER:
c    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c!DESIGN NOTES:
c    Original version by Kathy.Strabala@ssec.wisc.edu
c
c!END
c-----------------------------------------------------------------------

      implicit none

      include 'PGS_PC.f'

c ... Common statement for debug purposes
      common / bug / debug, h_output

c ... Arguments
      integer dataset_time(4)
      double precision tai_time,delta

c ... scalars
      integer iyr,imonth,iday,ihr,rtn,debug,h_output
      character*(PGSd_PC_VALUE_LENGTH_MAX) ccsds
      double precision End_Time


c ... external subroutines
      external message

c ... external functions
      integer pgs_td_utctotai,date_to_ccsds
      external pgs_td_utctotai,date_to_ccsds

c ... initialize
      ccsds = '  '
      End_Time = 0.0d0
      compare_times = -1

      iyr = dataset_time(1)
      imonth = dataset_time(2)
      iday = dataset_time(3)
      ihr = dataset_time(4)

c ... check inputs
      if (iyr .gt. 0 .and. imonth .gt. 0 .and. imonth .le. 12  .and.
     +  iday .gt. 0 .and. iday .le. 31 .and. ihr .gt. -1 .and.
     +  ihr .le. 23 .and. delta .gt. 0.0d0) then

c ...   convert the nise date/time to UTC
        rtn = date_to_ccsds(iyr,imonth,iday,ihr,0,0.0d0,ccsds)
        if (rtn .ne. 0) then
          call message( 'compare_times',
     +    'Error converting date/time from nise to CCSDS.'
     +    // char(10) // ' [OPERATOR ACTION: Notify SDST]',
     +    0, 1 )
        endif

c ...   now convert to tai
        rtn = pgs_td_utctotai(ccsds,End_Time)
        if (rtn .ne. 0) then
          call message( 'compare_times',
     +    'Error converting UTC Time to TAI.'
     +    // char(10) // ' [OPERATOR ACTION: Notify SDST]',
     +    0, 1 )
        endif

c ...   Make sure TAI's are not out of range
        if (End_Time .gt. 0.0d0 .and. tai_time .gt. 0.0d0) then

c ...     Now compare tai differences with delta
c rhucek 06/22/00: Replaced single test logic with 2-test logic to screen 
c separately at each ancillary data boundary. 
c         if (abs(End_Time - tai_time) .lt. (Delta + 1.0e-5)) then

          if ( (tai_time + Delta - End_Time .gt. -1) .and. 
     1         (End_Time + Delta - tai_time .gt.  1) ) then 
            compare_times = 0
          else
            compare_times = -1
          endif

        else

        call message( 'compare_times',
     +    'TAI times are less than zero.'
     +    // char(10) // ' [OPERATOR ACTION: Notify SDST]',
     +    0, 1 )
        endif

      else

        call message( 'compare_times',
     +  'Input date/time (yr,month,day,hour) is not valid.'
     +  // char(10) // ' [OPERATOR ACTION: Notify SDST]',
     +  0, 1 )
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(15x,'' compare_times.f '')')
        write(h_output,'(2x,'' date/time of anc data set: year, month,
     + day, hour'',/,4I10)') iyr,imonth,iday,ihr
        write(h_output,'(2x,'' tai times, delta'',/,3f20.5,/)')
     +       tai_time,End_Time,delta
      endif
c ................................................................



      END
