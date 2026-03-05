subroutine get_lut_xcal_modisa_ts( lut_file, status)
!
!f90
!description:
!   This subroutine reads core SDSs in the Aqua/MODIS RVS and polariztion 
!   LUT hdf file. The hdf file, containing data for Aqua over 
!   ocean, was provided by the NASA/GSFC OBPG (contact G. Meister & B. Franz).
!   file #1= 'xcal_modisa_axc05a_412.hdf'
!   file #2= 'xcal_modisa_axc05a_443.hdf'
!   file #3= 'xcal_modisa_axc05a_488.hdf'
!   includes:   rvs	
!   		coeff m12 (should not be used; OBPG advised to use pre-launch
!   		characterizations)
!               coeff m13 (should not be used due to same reason as for m12)
!   array dimension 
!  [time, bands, mirror side, detector, polinomial coeff(3rd-order)]
!  ---> [95, 2, 10, 4] : note M11, m12, m13 data only for 412nm (B8)
!  c.f. [200, 9, 2, 10, 4] --> Terra/MODIS M11, m12, m13 data structure (C051)
!  c.f. [117, 9, 2, 10, 4] --> Terra/MODIS M11, m12, m13 data structure (C006)
!
! The above three files were then merged into a single file:
!  "xcal_modisa_axc05a_Allbands_reformatted.hdf" to maintain the same structure
! as that for Terra/MODIS xcal file, except for different time spans 
! i.e., [95, 9, 2, 10, 4] --> This hdf file is listed as "412305" in the pcf file.
!
!   Written by Myeong-Jae Jeong (MJ) at GEST/UMBC & NASA GSFC
!              October 2010
!
!
   use GeneralAuxType
   use lut_arrays
   use hdf

   implicit none

   character(*),           intent (in)      :: lut_file
   integer,                intent (out)     :: status
   integer, dimension (2)      		    :: start, edge, stride

   integer(integer_fourbyte)		    :: rank, number_type, nattrs
   integer(integer_fourbyte)		    :: sds_id,sds_index,attr_index, hdfid
   integer(integer_fourbyte), dimension(5)  :: start5d, edge5d, stride5d
   character(len=64)			    :: sds_name, local_sds_name, attr_name, &
                                               scale_factor_name, offset_factor_name
   real(single)                             :: xhour      
   real(double)                             :: xcalsectmp 
   integer(integer_fourbyte)                :: it         
   integer                                  :: xxyear_tmp, xxday_tmp

   start5d  = (/ 0,0,0,0,0 /)
   edge5d   = (/ ntimer_aq,nbandlutr,nside,ndet,ncoeff /)
   stride5d = (/ 1,1,1,1,1 /)
    
!   write(*,*) 'Reading Aqua/MODIS polarization file from ',lut_file
!=======================================================================
! start the geolocation file
!=======================================================================
!   print *, 'lut_file',lut_file
   hdfid = sfstart(lut_file,DFACC_READ)
   status=sffinfo(hdfid, number_type, nattrs)
!   write(*,*) 'hdfid= ', hdfid
! wave 
   sds_name = 'wave'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
!   status=sfrdata(sds_id,start1d,stride1d,edge1d,angle_of_incident)
!  status=sfrdata(sds_id,0,1,9,xxwave) ! ori
   status=sfrdata(sds_id,0,1,nbandlutr,xxwave_aq) ! mj
   status = sfendacc(sds_id)
!  write(*,*) 'WVL= ', xxwave
! year
   sds_name = 'year'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
!  status=sfrdata(sds_id,0,1,200,xxyear) ! ori
   status=sfrdata(sds_id,0,1,ntimer_aq,xxyear_aq) ! mj
   status = sfendacc(sds_id)
!  write(*,*) 'year= ', xxyear
! day
   sds_name = 'day'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
!  status=sfrdata(sds_id,0,1,200,xxday) ! ori
   status=sfrdata(sds_id,0,1,ntimer_aq,xxday_aq) ! mj
   status = sfendacc(sds_id)
!  write(*,*) 'jday= ', xxday
! rvs
   sds_name = 'M11'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status=sfrdata(sds_id,start5d,stride5d,edge5d,xxrvs_aq)
!   print *, 'rvs status= ',status
   status = sfendacc(sds_id)

! m12  !--> data exists in the file, but will not be used !mj
   sds_name = 'm12'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status=sfrdata(sds_id,start5d,stride5d,edge5d,xxam12_aq)
!   print *, 'xxam12 status= ',status
   status = sfendacc(sds_id)

! m13  !--> data exists in the file, but will not be used !mj
   sds_name = 'm13'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   status=sfrdata(sds_id,start5d,stride5d,edge5d,xxam13_aq)
!   print *, 'xxam13 status= ',status
   status = sfendacc(sds_id)

! close this HDF file
   status = sfend(hdfid) 

   xhour=12.0  ! assume m12, m13, rvs were derived at 12UTC
   do it=1, ntimer_aq
!     if(it.gt.0.and.it.le.92) then  !(1st LUT applicable till Dec 2007) 
      if(it.gt.0.and.it.le.95) then
         if(xxyear(it).ge.1993.and.xxday(it).gt.0) then
            xxyear_tmp=xxyear(it)
            xxday_tmp=xxday(it)
            call xcal_yds2sec(xxyear_tmp,xxday_tmp,xhour,xcalsectmp)
            sectab_aq(it)=xcalsectmp
         else
!            print *, 'check-->', it, xxyear(it),xxday(it),xhour
         endif
!... ************************************************************
!... ************************************************************
!... Note: below is a temporary set up. should be revisited when 
!... processing data acquired after Oct 2010
!... vicarious calibration / polarization sensitivity parameters
!... are available for Y2002 - Y2010
!... ************************************************************
!... ************************************************************
!     else if(it.gt.92) then ! (1st LUT applicable till Dec 2007;Terra)
!        sectab(it)=sectab(92)+3600.0
      else if(it.gt.95) then
         sectab_aq(it)=sectab_aq(95)+3600.0
!        print *, it, sectab_aq(it)
      else
         sectab_aq(it)=sectab_aq(1)-3600.0
!        print *, it, sectab_aq(it)
      endif
   enddo

end subroutine get_lut_xcal_modisa_ts

