      SUBROUTINE THRESHOLDS_READ( VERSION )

C-------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Read all MOD35 thresholds from parameter file.  Also reads and
C     verifies version and satellite id information from thresholds 
C     file.
C
C !INPUT PARAMETERS:
C     None
C
C !OUTPUT PARAMETERS:
C     VERSION    Version string (RCS Id) for threshold parameter
C                file.
C
C     Threshold values defined in the following include files:
C     LandDay_desert_c_thr.inc
C     PolarDay_desert_c_thr.inc
C     LandDay_desert_thr.inc
C     PolarDay_desert_thr.inc
C     LandDay_thr.inc
C     PolarDay_land_thr.inc
C     LandDay_coast_thr.inc
C     PolarDay_coast_thr.inc
C     LandNite_thr.inc
C     PolarNite_land_thr.inc
C     PolarDay_snow_thr.inc
C     Day_snow_thr.inc
C     PolarNite_snow_thr.inc
C     Nite_snow_thr.inc
C     ocean_day_thr.inc
C     PolarDay_ocean_thr.inc
C     ocean_nite_thr.inc
C     PolarNite_ocean_thr.inc
C     shadows_thr.inc
C     snglntr_thr.inc
C     spatial_var_thr.inc
C     noncld_obs_chk.inc
C     snow_mask.inc
C     Antarctic_day_thr.inc
C     swc_ndvi.inc
C     land_restoral.inc
C
C !REVISION HISTORY:
C     Added version and satellite ID verification
C     02/07/02
C     Added variables (test thresholds) required for Collection 5
C     processing
C     06/04  R. Frey
C     Modified threshold (variable) lists for Collection 5b update.
C     10/04  R. Frey
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !DESIGN NOTES:
C     Original version by Liam.Gumley@ssec.wisc.edu
C     Version and satellite ID verification added by Rich Frey,
C     UW-Madison.
C
C !END
C--------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

c ... Arguments

      character*(*) version

c ... Local variables

      character*160 errmsg
      integer number

c ... MOD35 PCF number include file (defines LRN_THR_PAR)
      include 'mod35.inc'
                  
c ... Version and satellite ID internal verification
      include 'thresholds_ver.inc'
      include 'platform_name.inc'

c ... Threshold include files

      include 'LandDay_desert_c_thr.inc'
      include 'PolarDay_desert_c_thr.inc'
      include 'LandDay_desert_thr.inc'
      include 'PolarDay_desert_thr.inc'
      include 'LandDay_thr.inc'
      include 'PolarDay_land_thr.inc'
      include 'LandDay_coast_thr.inc'
      include 'PolarDay_coast_thr.inc'
      include 'LandNite_thr.inc'
      include 'PolarNite_land_thr.inc'
      include 'PolarDay_snow_thr.inc'
      include 'Day_snow_thr.inc'
      include 'PolarNite_snow_thr.inc'
      include 'Nite_snow_thr.inc'
      include 'ocean_day_thr.inc'
      include 'PolarDay_ocean_thr.inc'
      include 'ocean_nite_thr.inc'
      include 'PolarNite_ocean_thr.inc'
      include 'shadows_thr.inc'
      include 'snglntr_thr.inc'
      include 'spatial_var_thr.inc'
      include 'noncld_obs_chk.inc'
      include 'snow_mask.inc'
      include 'Antarctic_day_thr.inc'
      include 'swc_ndvi.inc'
      include 'land_restoral.inc'

c ... Threshold parameter file version (for output to file)
      call param_string( lrn_thr_par, 'rcs_id', version )
      
c ... Get threshold ID info and perform internal verification.
      call param_string(lrn_thr_par, ver_name, ver_string)

      if(platform_name(1:5) .eq. 'Terra' .or. 
     +  platform_name(1:5) .eq. 'terra' .or.
     +  platform_name(1:5) .eq. 'TERRA') then
        if((ver_string(8:12) .ne. thr_satid_terra) .or.  
     +     (ver_string(1:6) .ne. thr_ver_id_Terra)) then
          write(errmsg,'( ''Incorrect thresholds file in use'')')
          call message('thresholds_read',errmsg //
     +     ' [OPERATOR ACTION: Contact SDST]',0,2)
        end if
      else if(platform_name(1:4) .eq. 'Aqua' .or. 
     +  platform_name(1:4) .eq. 'aqua' .or.
     +  platform_name(1:4) .eq. 'AQUA') then
        if((ver_string(8:12) .ne. thr_satid_aqua) .or. 
     +     (ver_string(1:6) .ne. thr_ver_id_Aqua)) then
          write(errmsg,'( ''Incorrect thresholds file in use'')')
          call message('thresholds_read',errmsg //
     +     ' [OPERATOR ACTION: Contact SDST]',0,2)
        end if
      end if

c ... Daytime coastal desert thresholds.

      number = 4
      call param_real( lrn_thr_par, 'lds11_12hi_c', number, lds11_12hi_c )
      call param_real( lrn_thr_par, 'lds11_4hi_c',  number, lds11_4hi_c  )
      call param_real( lrn_thr_par, 'lds11_4lo_c',  number, lds11_4lo_c  )
      call param_real( lrn_thr_par, 'ldsco2_c',     number, ldsco2_c     )
      call param_real( lrn_thr_par, 'ldsh20_c',     number, ldsh20_c     )
      call param_real( lrn_thr_par, 'ldsref2_c',    number, ldsref2_c    )
      call param_real( lrn_thr_par, 'ldsref3_c',    number, ldsref3_c    )
      number = 2
      call param_real( lrn_thr_par, 'ldstci_c',    number, ldstci_c  )

c ... Daytime polar coastal desert thresholds.

      number = 4
      call param_real( lrn_thr_par, 'pds11_12hi_c', number, pds11_12hi_c )
      call param_real( lrn_thr_par, 'pds11_4hi_c',  number, pds11_4hi_c  )
      call param_real( lrn_thr_par, 'pds11_4lo_c',  number, pds11_4lo_c  )
      call param_real( lrn_thr_par, 'pdsh20_c',     number, pdsh20_c     )
      call param_real( lrn_thr_par, 'pdsref2_c',    number, pdsref2_c    )
      call param_real( lrn_thr_par, 'pdsref3_c',    number, pdsref3_c    )
      number = 2
      call param_real( lrn_thr_par, 'pdstci_c',    number, pdstci_c  )

c ... Daytime desert thresholds.

      number = 4
      call param_real( lrn_thr_par, 'lds11_12hi', number, lds11_12hi )
      call param_real( lrn_thr_par, 'lds11_4hi',  number, lds11_4hi  )
      call param_real( lrn_thr_par, 'lds11_4lo',  number, lds11_4lo  )
      call param_real( lrn_thr_par, 'ldsco2',     number, ldsco2     )
      call param_real( lrn_thr_par, 'ldsh20',     number, ldsh20     )
      call param_real( lrn_thr_par, 'ldsref2',    number, ldsref2    )
      call param_real( lrn_thr_par, 'ldsref3',    number, ldsref3    )
      number = 2
      call param_real( lrn_thr_par, 'ldstci',    number, ldstci    )

c ... Daytime polar desert thresholds.

      number = 4
      call param_real( lrn_thr_par, 'pds11_12hi', number, pds11_12hi )
      call param_real( lrn_thr_par, 'pds11_4hi',  number, pds11_4hi  )
      call param_real( lrn_thr_par, 'pds11_4lo',  number, pds11_4lo  )
      call param_real( lrn_thr_par, 'pdsh20',     number, pdsh20     )
      call param_real( lrn_thr_par, 'pdsref2',    number, pdsref2    )
      call param_real( lrn_thr_par, 'pdsref3',    number, pdsref3    )
      number = 2
      call param_real( lrn_thr_par, 'pdstci',    number, pdstci    )

c ... Daytime land thresholds.

      number = 1
      call param_real( lrn_thr_par, 'dl11_12hi', number, dl11_12hi )
      number = 4
      call param_real( lrn_thr_par, 'dl11_4lo',  number, dl11_4lo  )
      call param_real( lrn_thr_par, 'dlco2',     number, dlco2     )
      call param_real( lrn_thr_par, 'dlh20',     number, dlh20     )
      call param_real( lrn_thr_par, 'dlref1',    number, dlref1    )
      call param_real( lrn_thr_par, 'dlref3',    number, dlref3    )
      call param_real( lrn_thr_par, 'dlvrat',    number, dlvrat    )
      number = 2
      call param_real( lrn_thr_par, 'dltci',    number, dltci    )

c ... Daytime polar land thresholds.

      number = 1
      call param_real( lrn_thr_par, 'pdl11_12hi', number, pdl11_12hi )
      number = 4
      call param_real( lrn_thr_par, 'pdl11_4lo',  number, pdl11_4lo  )
      call param_real( lrn_thr_par, 'pdlh20',     number, pdlh20     )
      call param_real( lrn_thr_par, 'pdlref1',    number, pdlref1    )
      call param_real( lrn_thr_par, 'pdlref3',    number, pdlref3    )
      call param_real( lrn_thr_par, 'pdlvrat',    number, pdlvrat    )
      number = 2
      call param_real( lrn_thr_par, 'pdltci',    number, pdltci    )

c ... Daytime coastal land thresholds.

      number = 1
      call param_real( lrn_thr_par, 'dl11_12hi_t2', number, dl11_12hi_t2 )
      number = 4
      call param_real( lrn_thr_par, 'dl11_4lo_t2',  number, dl11_4lo_t2  )
      call param_real( lrn_thr_par, 'dlco2_t2',     number, dlco2_t2     )
      call param_real( lrn_thr_par, 'dlh20_t2',     number, dlh20_t2     )
      call param_real( lrn_thr_par, 'dlref1_t2',    number, dlref1_t2    )
      call param_real( lrn_thr_par, 'dlref3_t2',    number, dlref3_t2    )
      number = 2
      call param_real( lrn_thr_par, 'dltci_t2',    number, dltci_t2 )

c ... Daytime polar coastal land thresholds.

      number = 1
      call param_real( lrn_thr_par, 'pdl11_12hi_t2', number, pdl11_12hi_t2 )
      number = 4
      call param_real( lrn_thr_par, 'pdl11_4lo_t2',  number, pdl11_4lo_t2  )
      call param_real( lrn_thr_par, 'pdlh20_t2',     number, pdlh20_t2     )
      call param_real( lrn_thr_par, 'pdlref1_t2',    number, pdlref1_t2    )
      call param_real( lrn_thr_par, 'pdlref3_t2',    number, pdlref3_t2    )
      number = 2
      call param_real( lrn_thr_par, 'pdltci_t2',    number, pdltci_t2    )

c ... Nighttime land thresholds.

      number = 4
      call param_real( lrn_thr_par, 'nl4_12hi', number, nl4_12hi )
      call param_real( lrn_thr_par, 'nl4_12lo', number, nl4_12lo )
      call param_real( lrn_thr_par, 'nlco2',    number, nlco2    )
      call param_real( lrn_thr_par, 'nlh20',    number, nlh20    )
      call param_real( lrn_thr_par, 'nl7_11s',    number, nl7_11s  )
      call param_real( lrn_thr_par, 'nl_11_4l',    number, nl_11_4l  )
      call param_real( lrn_thr_par, 'nl_11_4h',    number, nl_11_4h  )
      call param_real( lrn_thr_par, 'nl_11_4m',    number, nl_11_4m  )
      number = 2
      call param_real( lrn_thr_par, 'bt_diff_bounds', number, bt_diff_bounds  )
      number = 1
      call param_real( lrn_thr_par, 'nl11_12hi', number, nl11_12hi)

c ... Nighttime polar land thresholds.

      number = 4
      call param_real( lrn_thr_par, 'pnlh20',    number, pnlh20    )
      number = 1
      call param_real( lrn_thr_par, 'pnl11_12hi', number, pnl11_12hi)

c ... Daytime polar snow thresholds.

      number = 4
      call param_real( lrn_thr_par, 'dpsh20',  number, dpsh20  )
      call param_real( lrn_thr_par, 'dpsref1', number, dpsref1 )
      call param_real( lrn_thr_par, 'dpsref3', number, dpsref3 )
      call param_real( lrn_thr_par, 'dps4_11l', number, dps4_11l )
      call param_real( lrn_thr_par, 'dps4_11h', number, dps4_11h )
      call param_real( lrn_thr_par, 'dps4_11m1', number, dps4_11m1 )
      call param_real( lrn_thr_par, 'dps4_11m2', number, dps4_11m2 )
      call param_real( lrn_thr_par, 'dps4_11m3', number, dps4_11m3 )
      call param_real( lrn_thr_par, 'bt_11_bnds3', number, bt_11_bnds3 )
      number = 2
      call param_real( lrn_thr_par, 'dpstci',    number, dpstci    )
      number = 1
      call param_real( lrn_thr_par, 'dps11_12hi', number, dps11_12hi    )
      call param_real( lrn_thr_par, 'dps11_12adj', number, dps11_12adj    )

c ... Antarctic day thresholds.

      number = 4
      call param_real( lrn_thr_par, 'ant4_11l', number, ant4_11l )
      call param_real( lrn_thr_par, 'ant4_11h', number, ant4_11h )
      call param_real( lrn_thr_par, 'ant4_11m1', number, ant4_11m1 )
      call param_real( lrn_thr_par, 'ant4_11m2', number, ant4_11m2 )
      call param_real( lrn_thr_par, 'ant4_11m3', number, ant4_11m3 )
      call param_real( lrn_thr_par, 'bt_11_bnds4', number, bt_11_bnds4 )
      call param_real( lrn_thr_par, 'anth20',  number, anth20  )

c ... Daytime snow thresholds.

      number = 4
      call param_real( lrn_thr_par, 'ds4_11', number, ds4_11 )
      call param_real( lrn_thr_par, 'ds4_11hel', number, ds4_11hel )
      call param_real( lrn_thr_par, 'dsco2',  number, dsco2  )
      call param_real( lrn_thr_par, 'dsh20',  number, dsh20  )
      call param_real( lrn_thr_par, 'dsref3', number, dsref3 )
      number = 2
      call param_real( lrn_thr_par, 'dstci',    number, dstci    )
      number = 1
      call param_real( lrn_thr_par, 'ds11_12hi', number, ds11_12hi )
      call param_real( lrn_thr_par, 'ds11_12adj', number, ds11_12adj )

c ... Nighttime polar snow thresholds.

      number = 4
      call param_real( lrn_thr_par, 'pn_4_12l',  number, pn_4_12l  )
      call param_real( lrn_thr_par, 'pn_4_12h',  number, pn_4_12h  )
      call param_real( lrn_thr_par, 'pn_4_12m1',  number, pn_4_12m1 )
      call param_real( lrn_thr_par, 'pn_4_12m2',  number, pn_4_12m2 )
      call param_real( lrn_thr_par, 'pn_4_12m3',  number, pn_4_12m3 )
      call param_real( lrn_thr_par, 'pn_7_11l',  number, pn_7_11l  )
      call param_real( lrn_thr_par, 'pn_7_11h',  number, pn_7_11h  )
      call param_real( lrn_thr_par, 'pn_7_11m1',  number, pn_7_11m1 )
      call param_real( lrn_thr_par, 'pn_7_11m2',  number, pn_7_11m2 )
      call param_real( lrn_thr_par, 'pn_7_11m3',  number, pn_7_11m3 )
      call param_real( lrn_thr_par, 'pn_7_11lw',  number, pn_7_11lw  )
      call param_real( lrn_thr_par, 'pn_7_11hw',  number, pn_7_11hw  )
      call param_real( lrn_thr_par, 'pn_7_11m1w',  number, pn_7_11m1w )
      call param_real( lrn_thr_par, 'pn_7_11m2w',  number, pn_7_11m2w )
      call param_real( lrn_thr_par, 'pn_7_11m3w',  number, pn_7_11m3w )
      call param_real( lrn_thr_par, 'pnsh20',    number, pnsh20    )
      call param_real( lrn_thr_par, 'pn_11_4l',  number, pn_11_4l  )
      call param_real( lrn_thr_par, 'pn_11_4h',  number, pn_11_4h  )
      call param_real( lrn_thr_par, 'pn_11_4m1',  number, pn_11_4m1 )
      call param_real( lrn_thr_par, 'pn_11_4m2',  number, pn_11_4m2 )
      call param_real( lrn_thr_par, 'pn_11_4m3',  number, pn_11_4m3 )
      call param_real( lrn_thr_par, 'bt_11_bounds', number, bt_11_bounds)
      call param_real( lrn_thr_par, 'bt_11_bnds2', number, bt_11_bnds2)
      number = 1
      call param_real( lrn_thr_par, 'pn65_11',    number, pn65_11   )
      call param_real( lrn_thr_par, 'pn13_11',    number, pn13_11   )
      call param_real( lrn_thr_par, 'pn7_11',    number, pn7_11   )
      call param_real( lrn_thr_par, 'pns11_12hi', number, pns11_12hi )
      call param_real( lrn_thr_par, 'pn11_12adj', number, pn11_12adj )

c ... Nighttime snow thresholds.

      number = 4
      call param_real( lrn_thr_par, 'ns11_4lo', number, ns11_4lo )
      call param_real( lrn_thr_par, 'ns4_12hi', number, ns4_12hi )
      call param_real( lrn_thr_par, 'nsco2',    number, nsco2    )
      call param_real( lrn_thr_par, 'nsh20',    number, nsh20    )
      number = 1
      call param_real( lrn_thr_par, 'n65_11',    number, n65_11   )
      call param_real( lrn_thr_par, 'ns11_12hi', number, ns11_12hi )
      call param_real( lrn_thr_par, 'ns11_12adj', number, ns11_12adj )

c ... Daytime ocean thresholds.

      number = 1
      call param_real( lrn_thr_par, 'do11_12hi', number, do11_12hi )
      number = 4
      call param_real( lrn_thr_par, 'do11_4lo',  number, do11_4lo  )
      call param_real( lrn_thr_par, 'dobt11',    number, dobt11    )
      call param_real( lrn_thr_par, 'doco2',     number, doco2     )
      call param_real( lrn_thr_par, 'doh20',     number, doh20     )
      call param_real( lrn_thr_par, 'doref2',    number, doref2    )
      call param_real( lrn_thr_par, 'doref3',    number, doref3    )
      call param_real( lrn_thr_par, 'dovrathi',  number, dovrathi  )
      call param_real( lrn_thr_par, 'dovratlo',  number, dovratlo  )
      number = 2
      call param_real( lrn_thr_par, 'dotci',    number, dotci    )

c ... Daytime polar ocean thresholds.

      number = 1
      call param_real( lrn_thr_par, 'pdo11_12hi', number, pdo11_12hi )
      number = 4
      call param_real( lrn_thr_par, 'pdo11_4lo',  number, pdo11_4lo  )
      call param_real( lrn_thr_par, 'pdobt11',    number, pdobt11    )
      call param_real( lrn_thr_par, 'pdoh20',     number, pdoh20     )
      call param_real( lrn_thr_par, 'pdoref2',    number, pdoref2    )
      call param_real( lrn_thr_par, 'pdoref3',    number, pdoref3    )
      call param_real( lrn_thr_par, 'pdovrathi',  number, pdovrathi  )
      call param_real( lrn_thr_par, 'pdovratlo',  number, pdovratlo  )
      number = 2
      call param_real( lrn_thr_par, 'pdotci',    number, pdotci    )

c ... Nighttime ocean thresholds.

      number = 1
      call param_real( lrn_thr_par, 'no11_12hi', number, no11_12hi )
      number = 4
      call param_real( lrn_thr_par, 'no11_4lo',  number, no11_4lo  )
      call param_real( lrn_thr_par, 'nobt11',    number, nobt11    )
      call param_real( lrn_thr_par, 'noco2',     number, noco2     )
      call param_real( lrn_thr_par, 'noh20',     number, noh20     )
      call param_real( lrn_thr_par, 'no86_73',   number, no86_73   )
      call param_real( lrn_thr_par, 'no_11var',   number, no_11var  )

c ... Nighttime polar ocean thresholds.

      number = 1
      call param_real( lrn_thr_par, 'pno11_12hi', number, pno11_12hi )
      number = 4
      call param_real( lrn_thr_par, 'pno11_4lo',  number, pno11_4lo  )
      call param_real( lrn_thr_par, 'pnobt11',    number, pnobt11    )
      call param_real( lrn_thr_par, 'pnoh20',     number, pnoh20     )
      call param_real( lrn_thr_par, 'pno86_73',   number, pno86_73   )
      call param_real( lrn_thr_par, 'pno_11var',   number, pno_11var  )

c ... Shadow Thresholds

      number = 2
      call param_real( lrn_thr_par, 'shadnir', number, shadnir )
      number = 1
      call param_real( lrn_thr_par, 'shavrat', number, shavrat )
      call param_real( lrn_thr_par, 'shad124', number, shad124 )

c ... Sun Glint Thresholds

      number = 2
      call param_real( lrn_thr_par, 'snglntv',   number, snglntv   )
      call param_real( lrn_thr_par, 'snglntvch', number, snglntvch )
      call param_real( lrn_thr_par, 'snglntvcl', number, snglntvcl )
      number = 1
      call param_real( lrn_thr_par, 'sg_tbdfl', number, sg_tbdfl )
      call param_real( lrn_thr_par, 'sg_tbdfh', number, sg_tbdfh )
      call param_real( lrn_thr_par, 'snglrat', number, snglrat )
      number = 4
      call param_real( lrn_thr_par, 'snglnt0', number, snglnt0 )
      call param_real( lrn_thr_par, 'snglnt10', number, snglnt10 )
      call param_real( lrn_thr_par, 'snglnt20', number, snglnt20 )
      call param_real( lrn_thr_par, 'snglnt_bounds', number, snglnt_bounds)

c ... Land Restoral Thresholds

      number = 1
      call param_real( lrn_thr_par, 'ldsr5_4_thr', number, ldsr5_4_thr )
      call param_real( lrn_thr_par, 'ldr5_4_thr', number, ldr5_4_thr )
      call param_real( lrn_thr_par, 'ld20m22', number, ld20m22 )
      call param_real( lrn_thr_par, 'ld22m31', number, ld22m31 )
      number = 3
      call param_real( lrn_thr_par, 'ldsbt11', number, ldsbt11 )
      call param_real( lrn_thr_par, 'ldsbt11bd', number, ldsbt11bd )
      call param_real( lrn_thr_par, 'lnbt11', number, lnbt11 )

c ... Day time ocean spatial variability threshold

      number = 1
      call param_real( lrn_thr_par, 'dovar11', number, dovar11 )

c ... Non-cloud Obstruction Thresholds

      number = 1
      call param_real( lrn_thr_par, 'nc_bt37',  number, nc_bt37 )
      call param_real( lrn_thr_par, 'nc37_11',  number, nc37_11 )
      call param_real( lrn_thr_par, 'nc21',     number, nc21 )
      call param_real( lrn_thr_par, 'nc11_12',  number, nc11_12 )
      call param_real( lrn_thr_par, 'ncrat',  number, ncrat )
      call param_real( lrn_thr_par, 'ncvrat',  number, ncvrat )
      call param_real( lrn_thr_par, 'ncsig',  number, ncsig )
     
c ... Snow Mask Thresholds

      number = 1
      call param_real( lrn_thr_par, 'sm_bt11',  number, sm_bt11 )
      call param_real( lrn_thr_par, 'sm_ndsi',  number, sm_ndsi )
      call param_real( lrn_thr_par, 'sm_ref2',  number, sm_ref2 )
      call param_real( lrn_thr_par, 'sm_ref3',  number, sm_ref3 )
      call param_real( lrn_thr_par, 'sm_co2',  number, sm_co2 )
      call param_real( lrn_thr_par, 'sm85_11',  number, sm85_11 )
      call param_real( lrn_thr_par, 'sm37_11',  number, sm37_11 )
      call param_real( lrn_thr_par, 'sm37_11hel',  number, sm37_11hel )
      call param_real( lrn_thr_par, 'sm_ndsi',  number, sm_ndsi )
      call param_real( lrn_thr_par, 'sm_mnir',  number, sm_mnir )

c ... Coast and shallow ocean ndvi thresholds.

      number = 2
      call param_real( lrn_thr_par, 'swc_ndvi', number, swc_ndvi )

      END
