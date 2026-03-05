/*
!C *********************************************************************
!Description:

  Integer function get_cm_thresholds.c
  Returns spectral test thresholds necessary for executing MODIS cloud
    mask.
  Called from collect_inputs.c

!Input arguments:
  none

  Inputs (definitions, constants, etc.) in thresholds.h

!Output arguments:
  none

  int return_code   successful completion is zero, otherwise non-zero
                    (returned through function call)

  Output thresholds stored in thresholds.h

!Revision History:
  R. Frey           05/2007


!Team-unique Header:

!References and Credits:

!END ******************************************************************/

/* Includes */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "thresholds.h"


int get_cm_thresholds(char *algData_dir_name) // GPC


{

/*  Declarations */

    extern int read_thresholds_file(char *);
    extern int get_thr_values(char thr_name[], int nvals, float *);

    int ret_code;
    int ret_code1;
    int return_code;
    int nvals;

/*  Initializations */

    return_code = 0;

/*  Get cloud mask thresholds */

/*  Read thresholds file and check the version. */

    /*ret_code1 = read_thresholds_file(sfc_dir_name);*/
    ret_code1 = read_thresholds_file(algData_dir_name); // GPC

    if(ret_code1 == 0) {

/*    Get values of desired thresholds. */

/*    Daytime land */

      nvals = 1;
      ret_code = get_thr_values("dl11_12hi", nvals, dl11_12hi);
      ret_code = get_thr_values("dl_ref3_tpw", nvals, dl_ref3_tpw);
      nvals = 4;
      ret_code = get_thr_values("dl11_4lo", nvals, dl11_4lo);
      ret_code = get_thr_values("dlco2", nvals, dlco2);
      ret_code = get_thr_values("dlh20", nvals, dlh20);
      ret_code = get_thr_values("dlref1", nvals, dlref1);
      ret_code = get_thr_values("dlref3", nvals, dlref3);
      ret_code = get_thr_values("dlvrat", nvals, dlvrat);
      nvals = 2;
      ret_code = get_thr_values("dltci", nvals, dltci);

/*    Daytime polar land */

      nvals = 1;
      ret_code = get_thr_values("pdl11_12hi", nvals, pdl11_12hi);
      ret_code = get_thr_values("pdl_ref3_tpw", nvals, pdl_ref3_tpw);
      nvals = 4;
      ret_code = get_thr_values("pdl11_4lo", nvals, pdl11_4lo);
      ret_code = get_thr_values("pdlh20", nvals, pdlh20);
      ret_code = get_thr_values("pdlref1", nvals, pdlref1);
      ret_code = get_thr_values("pdlref3", nvals, pdlref3);
      ret_code = get_thr_values("pdlvrat", nvals, pdlvrat);
      nvals = 2;
      ret_code = get_thr_values("pdltci", nvals, pdltci);

/*    Daytime coastal land */

      nvals = 1;
      ret_code = get_thr_values("dl11_12hi_t2", nvals, dl11_12hi_t2);
      ret_code = get_thr_values("dl_ref3_tpw_t2", nvals, dl_ref3_tpw_t2);
      nvals = 4;
      ret_code = get_thr_values("dl11_4lo_t2", nvals, dl11_4lo_t2);
      ret_code = get_thr_values("dlco2_t2", nvals, dlco2_t2);
      ret_code = get_thr_values("dlh20_t2", nvals, dlh20_t2);
      ret_code = get_thr_values("dlref1_t2", nvals, dlref1_t2);
      ret_code = get_thr_values("dlref3_t2", nvals, dlref3_t2);
      nvals = 2;
      ret_code = get_thr_values("dltci_t2", nvals, dltci_t2);

/*    Daytime polar coastal land */

      nvals = 1;
      ret_code = get_thr_values("pdl11_12hi_t2", nvals, pdl11_12hi_t2);
      ret_code = get_thr_values("pdl_ref3_tpw_t2", nvals, pdl_ref3_tpw_t2);
      nvals = 4;
      ret_code = get_thr_values("pdl11_4lo_t2", nvals, pdl11_4lo_t2);
      ret_code = get_thr_values("pdlh20_t2", nvals, pdlh20_t2);
      ret_code = get_thr_values("pdlref1_t2", nvals, pdlref1_t2);
      ret_code = get_thr_values("pdlref3_t2", nvals, pdlref3_t2);
      nvals = 2;
      ret_code = get_thr_values("pdltci_t2", nvals, pdltci_t2);

/*    Daytime ocean */

      nvals = 1;
      ret_code = get_thr_values("do11_12hi", nvals, do11_12hi);
      nvals = 4;
      ret_code = get_thr_values("do11_4lo", nvals, do11_4lo);
      ret_code = get_thr_values("dobt11", nvals, dobt11);
      ret_code = get_thr_values("doco2", nvals, doco2);
      ret_code = get_thr_values("doh20", nvals, doh20);
      ret_code = get_thr_values("doref2", nvals, doref2);
      ret_code = get_thr_values("doref3", nvals, doref3);
      ret_code = get_thr_values("dovrathi", nvals, dovrathi);
      ret_code = get_thr_values("dovratlo", nvals, dovratlo);
      ret_code = get_thr_values("do11_86", nvals, do11_86);
      ret_code = get_thr_values("dosst", nvals, dosst);
      nvals = 2;
      ret_code = get_thr_values("dotci", nvals, dotci);

/*    Daytime polar ocean */

      nvals = 1;
      ret_code = get_thr_values("pdo11_12hi", nvals, pdo11_12hi);
      nvals = 4;
      ret_code = get_thr_values("pdo11_4lo", nvals, pdo11_4lo);
      ret_code = get_thr_values("pdobt11", nvals, pdobt11);
      ret_code = get_thr_values("pdoh20", nvals, pdoh20);
      ret_code = get_thr_values("pdoref2", nvals, pdoref2);
      ret_code = get_thr_values("pdoref3", nvals, pdoref3);
      ret_code = get_thr_values("pdovrathi", nvals, pdovrathi);
      ret_code = get_thr_values("pdovratlo", nvals, pdovratlo);
      ret_code = get_thr_values("pdo11_86", nvals, pdo11_86);
      ret_code = get_thr_values("pdosst", nvals, pdosst);
      nvals = 2;
      ret_code = get_thr_values("pdotci", nvals, pdotci);

/*    Daytime ocean spatial variability */

      nvals = 1;
      ret_code = get_thr_values("dovar11", nvals, dovar11);

/*    Daytime desert */

      nvals = 4;
      ret_code = get_thr_values("lds11_12hi", nvals, lds11_12hi);
      ret_code = get_thr_values("lds11_4hi", nvals, lds11_4hi);
      ret_code = get_thr_values("lds11_4lo", nvals, lds11_4lo);
      ret_code = get_thr_values("ldsco2", nvals, ldsco2);
      ret_code = get_thr_values("ldsh20", nvals, ldsh20);
      ret_code = get_thr_values("ldsref2", nvals, ldsref2);
      ret_code = get_thr_values("ldsref3", nvals, ldsref3);
      ret_code = get_thr_values("ldsgemi0", nvals, ldsgemi0);
      ret_code = get_thr_values("ldsgemi1", nvals, ldsgemi1);
      ret_code = get_thr_values("ldsgemi2", nvals, ldsgemi2);
      nvals = 2;
      ret_code = get_thr_values("ldstci", nvals, ldstci);
      nvals = 1;
      ret_code = get_thr_values("lds_ref3_tpw", nvals, lds_ref3_tpw);
      ret_code = get_thr_values("ldsbt1", nvals, ldsbt1);
      ret_code = get_thr_values("lds_ndvi", nvals, lds_ndvi);

/*    Daytime polar desert */

      nvals = 4;
      ret_code = get_thr_values("pds11_12hi", nvals, pds11_12hi);
      ret_code = get_thr_values("pds11_4hi", nvals, pds11_4hi);
      ret_code = get_thr_values("pds11_4lo", nvals, pds11_4lo);
      ret_code = get_thr_values("pdsh20", nvals, pdsh20);
      ret_code = get_thr_values("pdsref2", nvals, pdsref2);
      ret_code = get_thr_values("pdsref3", nvals, pdsref3);
      ret_code = get_thr_values("pdsgemi0", nvals, pdsgemi0);
      ret_code = get_thr_values("pdsgemi1", nvals, pdsgemi1);
      ret_code = get_thr_values("pdsgemi2", nvals, pdsgemi2);
      nvals = 2;
      ret_code = get_thr_values("pdstci", nvals, pdstci);
      nvals = 1;
      ret_code = get_thr_values("pds_ref3_tpw", nvals, pds_ref3_tpw);
      ret_code = get_thr_values("pdsbt1", nvals, pdsbt1);

/*    Daytime coastal desert */

      nvals = 4;
      ret_code = get_thr_values("lds11_12hi_c", nvals, lds11_12hi_c);
      ret_code = get_thr_values("lds11_4hi_c", nvals, lds11_4hi_c);
      ret_code = get_thr_values("lds11_4lo_c", nvals, lds11_4lo_c);
      ret_code = get_thr_values("ldsco2_c", nvals, ldsco2_c);
      ret_code = get_thr_values("ldsh20_c", nvals, ldsh20_c);
      ret_code = get_thr_values("ldsref2_c", nvals, ldsref2_c);
      ret_code = get_thr_values("ldsref3_c", nvals, ldsref3_c);
      nvals = 2;
      ret_code = get_thr_values("ldstci_c", nvals, ldstci_c);
      nvals = 1;
      ret_code = get_thr_values("lds_ref3_tpw_c", nvals, lds_ref3_tpw_c);
      ret_code = get_thr_values("ldsbt1_c", nvals, ldsbt1_c);
      ret_code = get_thr_values("lds_ndvi_c", nvals, lds_ndvi_c);

/*    Daytime polar coastal desert */

      nvals = 4;
      ret_code = get_thr_values("pds11_12hi_c", nvals, pds11_12hi_c);
      ret_code = get_thr_values("pds11_4hi_c", nvals, pds11_4hi_c);
      ret_code = get_thr_values("pds11_4lo_c", nvals, pds11_4lo_c);
      ret_code = get_thr_values("pdsh20_c", nvals, pdsh20_c);
      ret_code = get_thr_values("pdsref2_c", nvals, pdsref2_c);
      ret_code = get_thr_values("pdsref3_c", nvals, pdsref3_c);
      nvals = 2;
      ret_code = get_thr_values("pdstci_c", nvals, pdstci_c);
      nvals = 1;
      ret_code = get_thr_values("pds_ref3_tpw_c", nvals, pds_ref3_tpw_c);
      ret_code = get_thr_values("pdsbt1_c", nvals, pdsbt1_c);

/*    Nighttime land */

      nvals = 4;
      ret_code = get_thr_values("nl4_12hi", nvals, nl4_12hi);
      ret_code = get_thr_values("nlco2", nvals, nlco2);
      ret_code = get_thr_values("nlh20", nvals, nlh20);
      ret_code = get_thr_values("nl7_11s", nvals, nl7_11s);
      ret_code = get_thr_values("nl_11_4l", nvals, nl_11_4l);
      ret_code = get_thr_values("nl_11_4h", nvals, nl_11_4h);
      ret_code = get_thr_values("nl_11_4m", nvals, nl_11_4m);
      nvals = 2;
      ret_code = get_thr_values("bt_diff_bounds", nvals, bt_diff_bounds);
      ret_code = get_thr_values("nl_df1", nvals, nl_df1);
      ret_code = get_thr_values("nl_df2", nvals, nl_df2);
      nvals = 1;
      ret_code = get_thr_values("nl11_12hi", nvals, nl11_12hi);
      ret_code = get_thr_values("nl_sfct1", nvals, nl_sfct1);
      ret_code = get_thr_values("nl_sfct2", nvals, nl_sfct2);
      ret_code = get_thr_values("nlbt1", nvals, nlbt1);
      ret_code = get_thr_values("nl_lat", nvals, nl_lat);
      ret_code = get_thr_values("nl_ndvi", nvals, nl_ndvi);
      ret_code = get_thr_values("nl_ndvi_Aust", nvals, nl_ndvi_Aust);
      ret_code = get_thr_values("nl_intadj", nvals, nl_intadj);
      ret_code = get_thr_values("nl_btd1", nvals, nl_btd1);

/*    Nighttime polar land */

      nvals = 4;
      ret_code = get_thr_values("pnlh20", nvals, pnlh20);
      nvals = 1;
      ret_code = get_thr_values("pnl11_12hi", nvals, pnl11_12hi);
      ret_code = get_thr_values("pnlbt1", nvals, pnlbt1);
      ret_code = get_thr_values("pnlbt2", nvals, pnlbt2);

/*    Nighttime ocean */

      nvals = 1;
      ret_code = get_thr_values("no11_12hi", nvals, no11_12hi);
      ret_code = get_thr_values("no_intadj", nvals, no_intadj);
      nvals = 4;
      ret_code = get_thr_values("no11_4lo", nvals, no11_4lo);
      ret_code = get_thr_values("nobt11", nvals, nobt11);
      ret_code = get_thr_values("noco2", nvals, noco2);
      ret_code = get_thr_values("noh20", nvals, noh20);
      ret_code = get_thr_values("no86_73", nvals, no86_73);
      ret_code = get_thr_values("no_11var", nvals, no_11var);
      ret_code = get_thr_values("no11_86", nvals, no11_86);
      ret_code = get_thr_values("nosst", nvals, nosst);

/*    Nighttime polar ocean */

      nvals = 1;
      ret_code = get_thr_values("pno11_12hi", nvals, pno11_12hi);
      ret_code = get_thr_values("pnobt1", nvals, pnobt1);
      ret_code = get_thr_values("pno_intadj", nvals, pno_intadj);
      nvals = 4;
      ret_code = get_thr_values("pno11_4lo", nvals, pno11_4lo);
      ret_code = get_thr_values("pnobt11", nvals, pnobt11);
      ret_code = get_thr_values("pnoh20", nvals, pnoh20);
      ret_code = get_thr_values("pno86_73", nvals, pno86_73);
      ret_code = get_thr_values("pno_11var", nvals, pno_11var);
      ret_code = get_thr_values("pno11_86", nvals, pno11_86);
      ret_code = get_thr_values("pnosst", nvals, pnosst);

/*    Daytime snow */

      nvals = 1;
      ret_code = get_thr_values("ds11_12hi", nvals, ds11_12hi);
      ret_code = get_thr_values("ds11_12adj", nvals, ds11_12adj);
      ret_code = get_thr_values("ds_ref3_tpw", nvals, ds_ref3_tpw);
      nvals = 4;
      ret_code = get_thr_values("ds4_11", nvals, ds4_11);
      ret_code = get_thr_values("ds4_11hel", nvals, ds4_11hel);
      ret_code = get_thr_values("dsco2", nvals, dsco2);
      ret_code = get_thr_values("dsh20", nvals, dsh20);
      ret_code = get_thr_values("dsref3", nvals, dsref3);
      nvals = 2;
      ret_code = get_thr_values("dstci", nvals, dstci);

/*    Daytime polar snow */

      nvals = 1;
      ret_code = get_thr_values("dps11_12hi", nvals, dps11_12hi);
      ret_code = get_thr_values("dps11_12adj", nvals, dps11_12adj);
      ret_code = get_thr_values("dps_ref3_tpw", nvals, dps_ref3_tpw);
      ret_code = get_thr_values("dpsbt1", nvals, dpsbt1);
      ret_code = get_thr_values("dpsbt2", nvals, dpsbt2);
      nvals = 4;
      ret_code = get_thr_values("dpsh20", nvals, dpsh20);
      ret_code = get_thr_values("dpsref3", nvals, dpsref3);
      ret_code = get_thr_values("dps4_11l", nvals, dps4_11l);
      ret_code = get_thr_values("dps4_11h", nvals, dps4_11h);
      ret_code = get_thr_values("dps4_11m1", nvals, dps4_11m1);
      ret_code = get_thr_values("dps4_11m2", nvals, dps4_11m2);
      ret_code = get_thr_values("dps4_11m3", nvals, dps4_11m3);
      ret_code = get_thr_values("bt_11_bnds1", nvals, bt_11_bnds1);
      ret_code = get_thr_values("bt_11_bnds3", nvals, bt_11_bnds3);
      nvals = 2;
      ret_code = get_thr_values("dpstci", nvals, dpstci);

/*    Antarctic day */

      nvals = 4;
      ret_code = get_thr_values("anth20", nvals, anth20);
      ret_code = get_thr_values("ant4_11l", nvals, ant4_11l);
      ret_code = get_thr_values("ant4_11h", nvals, ant4_11h);
      ret_code = get_thr_values("ant4_11m1", nvals, ant4_11m1);
      ret_code = get_thr_values("ant4_11m2", nvals, ant4_11m2);
      ret_code = get_thr_values("ant4_11m3", nvals, ant4_11m3);
      ret_code = get_thr_values("bt_11_bnds5", nvals, bt_11_bnds5);
      ret_code = get_thr_values("bt_11_bnds4", nvals, bt_11_bnds4);
      nvals = 1;
      ret_code = get_thr_values("antbt1", nvals, antbt1);
      ret_code = get_thr_values("antbt2", nvals, antbt2);

/*    Nighttime polar snow */

      nvals = 1;
      ret_code = get_thr_values("pns11_12hi", nvals, pns11_12hi);
      ret_code = get_thr_values("pn11_12adj", nvals, pn11_12adj);
      ret_code = get_thr_values("pn65_11", nvals, pn65_11);
      ret_code = get_thr_values("pn13_11", nvals, pn13_11);
      ret_code = get_thr_values("pn7_11", nvals, pn7_11);
      ret_code = get_thr_values("pnbt1", nvals, pnbt1);
      ret_code = get_thr_values("pnbt2", nvals, pnbt2);
      ret_code = get_thr_values("pnbt3", nvals, pnbt3);
      nvals = 4;
      ret_code = get_thr_values("pn_4_12l", nvals, pn_4_12l);
      ret_code = get_thr_values("pn_4_12h", nvals, pn_4_12h);
      ret_code = get_thr_values("pn_4_12m1", nvals, pn_4_12m1);
      ret_code = get_thr_values("pn_4_12m2", nvals, pn_4_12m2);
      ret_code = get_thr_values("pn_4_12m3", nvals, pn_4_12m3);
      ret_code = get_thr_values("pn_7_11l", nvals, pn_7_11l);
      ret_code = get_thr_values("pn_7_11h", nvals, pn_7_11h);
      ret_code = get_thr_values("pn_7_11m1", nvals, pn_7_11m1);
      ret_code = get_thr_values("pn_7_11m2", nvals, pn_7_11m2);
      ret_code = get_thr_values("pn_7_11m3", nvals, pn_7_11m3);
      ret_code = get_thr_values("pn_7_11lw", nvals, pn_7_11lw);
      ret_code = get_thr_values("pn_7_11hw", nvals, pn_7_11hw);
      ret_code = get_thr_values("pn_7_11m1w", nvals, pn_7_11m1w);
      ret_code = get_thr_values("pn_7_11m2w", nvals, pn_7_11m2w);
      ret_code = get_thr_values("pn_7_11m3w", nvals, pn_7_11m3w);
      ret_code = get_thr_values("pnsh20", nvals, pnsh20);
      ret_code = get_thr_values("pn_11_4l", nvals, pn_11_4l);
      ret_code = get_thr_values("pn_11_4h", nvals, pn_11_4h);
      ret_code = get_thr_values("pn_11_4m1", nvals, pn_11_4m1);
      ret_code = get_thr_values("pn_11_4m2", nvals, pn_11_4m2);
      ret_code = get_thr_values("pn_11_4m3", nvals, pn_11_4m3);
      ret_code = get_thr_values("bt_11_bounds", nvals, bt_11_bounds);
      ret_code = get_thr_values("bt_11_bnds2", nvals, bt_11_bnds2);

/*    Nighttime snow */

      nvals = 1;
      ret_code = get_thr_values("ns11_12hi", nvals, ns11_12hi);
      ret_code = get_thr_values("ns11_12adj", nvals, ns11_12adj);
      ret_code = get_thr_values("n65_11", nvals, n65_11);
      ret_code = get_thr_values("nsbt1", nvals, nsbt1);
      ret_code = get_thr_values("nsbt2", nvals, nsbt2);
      ret_code = get_thr_values("nsbt3", nvals, nsbt3);
      nvals = 4;
      ret_code = get_thr_values("ns11_4lo", nvals, ns11_4lo);
      ret_code = get_thr_values("ns4_12hi", nvals, ns4_12hi);
      ret_code = get_thr_values("nsco2", nvals, nsco2);
      ret_code = get_thr_values("nsh20", nvals, nsh20);

/*    Sun glint */

      nvals = 1;
      ret_code = get_thr_values("sg_tbdfl", nvals, sg_tbdfl);
      ret_code = get_thr_values("sg_tbdfh", nvals, sg_tbdfh);
      ret_code = get_thr_values("snglrat", nvals, snglrat);
      nvals = 2;
      ret_code = get_thr_values("snglntv", nvals, snglntv);
      ret_code = get_thr_values("snglntvch", nvals, snglntvch);
      ret_code = get_thr_values("snglntvcl", nvals, snglntvcl);
      nvals = 4;
      ret_code = get_thr_values("snglnt0", nvals, snglnt0);
      ret_code = get_thr_values("snglnt10", nvals, snglnt10);
      ret_code = get_thr_values("snglnt20", nvals, snglnt20);
      ret_code = get_thr_values("snglnt_bounds", nvals, snglnt_bounds);

/*    Land restoral */

      nvals = 1;
      ret_code = get_thr_values("ldsr5_4_thr", nvals, ldsr5_4_thr);
      ret_code = get_thr_values("ldr5_4_thr", nvals, ldr5_4_thr);
      ret_code = get_thr_values("ld20m22", nvals, ld20m22);
      ret_code = get_thr_values("ld22m31", nvals, ld22m31);
      nvals = 3;
      ret_code = get_thr_values("ldsbt11", nvals, ldsbt11);
      ret_code = get_thr_values("ldsbt11bd", nvals, ldsbt11bd);
      ret_code = get_thr_values("lnbt11", nvals, lnbt11);

/*    Non-cloud obstruction */

      nvals = 1;
      ret_code = get_thr_values("nc_bt37", nvals, nc_bt37);
      ret_code = get_thr_values("nc37_11", nvals, nc37_11);
      ret_code = get_thr_values("nc21", nvals, nc21);
      ret_code = get_thr_values("nc11_12", nvals, nc11_12);
      ret_code = get_thr_values("ncrat", nvals, ncrat);
      ret_code = get_thr_values("ncvrat", nvals, ncvrat);
      ret_code = get_thr_values("ncsig", nvals, ncsig);

      ret_code = get_thr_values("ncd_ldrat2", nvals, ncd_ldrat2);
      nvals = 2;
      ret_code = get_thr_values("ncd_lddf1", nvals, ncd_lddf1);
      ret_code = get_thr_values("ncd_ldm26", nvals, ncd_ldm26);
      ret_code = get_thr_values("ncd_ldrat1", nvals, ncd_ldrat1);
      nvals = 3;
      ret_code = get_thr_values("ncd_lddf2", nvals, ncd_lddf2);
      nvals= 1;
      ret_code = get_thr_values("ncd_wdsig2", nvals, ncd_wdsig2);
      ret_code = get_thr_values("ncd_wdm03", nvals, ncd_wdm03);
      nvals = 2;
      ret_code = get_thr_values("ncd_wdsst", nvals, ncd_wdsst);
      ret_code = get_thr_values("ncd_wdrat4", nvals, ncd_wdrat4);
      nvals = 4;
      ret_code = get_thr_values("ncd_wddf2", nvals, ncd_wddf2);
      ret_code = get_thr_values("ncd_wdrat3", nvals, ncd_wdrat3);
      ret_code = get_thr_values("ncd_wddf1", nvals, ncd_wddf1);
      nvals = 1;
      ret_code = get_thr_values("ncd_lnm31", nvals, ncd_lnm31);
      ret_code = get_thr_values("ncd_lnsig31", nvals, ncd_lnsig31);
      ret_code = get_thr_values("ncd_lndf1", nvals, ncd_lndf1);
      ret_code = get_thr_values("ncd_lndf2", nvals, ncd_lndf2);
      nvals = 2;
      ret_code = get_thr_values("ncd_lndf3", nvals, ncd_lndf3);
      nvals = 1;
      ret_code = get_thr_values("ncd_wnsig31", nvals, ncd_wnsig31);
      nvals = 2;
      ret_code = get_thr_values("ncd_wnsst", nvals, ncd_wnsst);
      nvals = 3;
      ret_code = get_thr_values("ncd_wndf1", nvals, ncd_wndf1);
      ret_code = get_thr_values("ncd_wndf3", nvals, ncd_wndf3);

/*    Snow mask */

      nvals = 1;
      ret_code = get_thr_values("sm_bt11", nvals, sm_bt11);
      ret_code = get_thr_values("sm_ndsi", nvals, sm_ndsi);
      ret_code = get_thr_values("sm_ref2", nvals, sm_ref2);
      ret_code = get_thr_values("sm_ref3", nvals, sm_ref3);
      ret_code = get_thr_values("sm_co2", nvals, sm_co2);
      ret_code = get_thr_values("sm85_11", nvals, sm85_11);
      ret_code = get_thr_values("sm37_11", nvals, sm37_11);
      ret_code = get_thr_values("sm37_11hel", nvals, sm37_11hel);
      ret_code = get_thr_values("sm_mnir", nvals, sm_mnir);

/*    Coast and shallow water NDVI */

      nvals = 2;
      ret_code = get_thr_values("swc_ndvi", nvals, swc_ndvi);

/*    Thin cirrus check on TPW */

      nvals = 1;
      ret_code = get_thr_values("tci_ref3_tpw", nvals, tci_ref3_tpw);

/*    Thresholds for calculating 0.66, 0.413 um thresholds */

      nvals = 1;
      ret_code = get_thr_values("des_ndvi", nvals, des_ndvi);
      ret_code = get_thr_values("fill_ndvi1", nvals, fill_ndvi1);
      ret_code = get_thr_values("fill_ndvi2", nvals, fill_ndvi2);
      ret_code = get_thr_values("ndvi_bnd1", nvals, ndvi_bnd1);
      ret_code = get_thr_values("ndvi_bnd2", nvals, ndvi_bnd2);
      ret_code = get_thr_values("thr_adj_fac_desert", nvals, thr_adj_fac_desert);
      ret_code = get_thr_values("thr_adj_fac_land", nvals, thr_adj_fac_land);

/*    Daytime land 0.66, 0.413 um NDVI and scattering angle dependent thresholds */

      nvals = 12;
      ret_code = get_thr_values("b1ndvi0", nvals, b1ndvi0);
      ret_code = get_thr_values("b1ndvi1", nvals, b1ndvi1);
      ret_code = get_thr_values("b1ndvi2", nvals, b1ndvi2);
      ret_code = get_thr_values("b1ndvi3", nvals, b1ndvi3);
      ret_code = get_thr_values("b1ndvi4", nvals, b1ndvi4);
      ret_code = get_thr_values("b1ndvi5", nvals, b1ndvi5);
      ret_code = get_thr_values("b1ndvi6", nvals, b1ndvi6);
      ret_code = get_thr_values("b1ndvi7", nvals, b1ndvi7);
      ret_code = get_thr_values("b1ndvi8", nvals, b1ndvi8);
      ret_code = get_thr_values("b1ndvi9", nvals, b1ndvi9);

      ret_code = get_thr_values("b8ndvi0", nvals, b8ndvi0);
      ret_code = get_thr_values("b8ndvi1", nvals, b8ndvi1);
      ret_code = get_thr_values("b8ndvi2", nvals, b8ndvi2);

    }

    else {

      return_code = ret_code1;

    }

    return return_code;

}
