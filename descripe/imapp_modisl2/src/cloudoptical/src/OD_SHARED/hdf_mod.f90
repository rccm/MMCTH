!****************************************************************************
!* NCSA HDF                                                            
!* Software Development Group                                          
!* National Center for Supercomputing Applications                     
!* University of Illinois at Urbana-Champaign                          
!* 605 E. Springfield, Champaign IL 61820                              
!*                                                                     
!* For conditions of distribution and use, see the accompanying        
!* hdf/COPYING file.                                                   
!*                                                                     
!****************************************************************************
!
! *-----------------------------------------------------------------------------
! * File:       hdf.inc
! * Purpose:    Fortran header file for HDF routines
! * Contents: 
! *     Tag definitions
! *     Error return codes
! *     Logical constants
! * Remarks: This file can be included with Fortran user programs.  As a
! *          general rule, don't use DFNT constants that don't include a
! *          number in their name.  E.g., don't use DFNT_FLOAT, use
! *          DFNT_FLOAT32 or DFNT_FLOAT64.  The DFNT constants that don't
! *          include numbers are for backward compatibility only.  Also,
! *          there are no current plans to support 128-bit number types.
! *          For more information about constants in this file, see the
! *          equivalent constant declarations in the C include file 'hdf.h'
! *---------------------------------------------------------------------XXXXXX*
! !F90
!
! !Description:
!   Fortran 90 version of the hdf.inc include file supplied with the 
!     HDF version 4.0.1 release. I've changed the syntax to look more like 
!     Fortran 90, but more importantly I've added interface definitions (like
!     function prototypes) for many of the functions in the SD API. 
!  
! !Input Parameters:
!   none
!
! !Output Parameters:
!    All the constants defined in hdf.inc
! 
! !Revision History:
!
! Revision 1.6  1999/10/12  12:35:14  EGMoody
! Added additional API definitions to this include file.
! Modified definition of SFscatt
!
! Revision 1.5  1997/12/30  15:29:09  pincus
! Moved definitions of netCDF parameters MAX_NC_DIMS, etc. into hdf.f90.
! Updated values to correspond with hlimits.h file in HDF include directory.
! Changed names.
!
! Revision 1.4  1997/11/20  21:52:55  pincus
! Added interface for SFscatt, SFsnatt.
!
! Revision 1.3  1997/11/03  23:26:12  pincus
! Cosmetic changes only: reformatted to comply with ANSI standard of
! 39 continuation lines max.
!
! Revision 1.2  1997/10/23  19:07:17  pincus
! Changed attr_index to intent( in) in SFgainfo interface; added interface
! for SFrcatt.
!
! Revision 1.1  1997/07/14  23:31:18  pincus
! Initial revision
!
!
! !Team-Unique Header:
!   Cloud Retrieval Group, NASA Goddard Space Flight Center
!
! !References and Credits:
!   Written by
!    Robert Pincus
!    Climate and Radiation Branch, Code 913
!    NASA/GSFC
!    Greenbelt MD 20771
!    Robert.Pincus@gsfc.nasa.gov
!
! !Design Notes:
!
! !END
module hdf_mod
  implicit none

!       Error Return Codes 
  integer, parameter :: DFE_NOERROR            =   0, &
                        DFE_NONE               =   0, &
                        DFE_FNF                =  -1, &
                        DFE_DENIED             =  -2, &
                        DFE_ALROPEN            =  -3, &
                        DFE_TOOMANY            =  -4, &
                        DFE_BADNAME            =  -5, &
                        DFE_BADACC             =  -6, &
                        DFE_BADOPEN            =  -7, &
                        DFE_NOTOPEN            =  -8, &
                        DFE_CANTCLOSE          =  -9, &
                        DFE_DFNULL             = -10, &
                        DFE_ILLTYPE            = -11, &
                        DFE_UNSUPPORTED        = -12, &
                        DFE_BADDDLIST          = -13, &
                        DFE_NOTDFFILE          = -14, &
                        DFE_SEEDTWICE          = -15, &
                        DFE_NOSPACE            = -16, &
                        DFE_NOSUCHTAG          = -17, &
                        DFE_READERROR          = -18, &
                        DFE_WRITEERROR         = -19, &
                        DFE_SEEKERROR          = -20, &
                        DFE_NOFREEDD           = -21, &
                        DFE_BADTAG             = -22, &
                        DFE_BADREF             = -23, &
                        DFE_RDONLY             = -24, &
                        DFE_BADCALL            = -25, &
                        DFE_BADPTR             = -26, &
                        DFE_BALEN              = -27, &
                        DFE_BADSEEK            = -28, &
                        DFE_NOMATCH            = -29, &
                        DFE_NOTINSET           = -30, &
                        DFE_BADDIM             = -31, &
                        DFE_BADOFFSET          = -32, &
                        DFE_BADSCHEME          = -33, &
                        DFE_NODIM              = -34, &
                        DFE_NOTENOUGH          = -35, &
                        DFE_NOVALS             = -36, &
                        DFE_CORRUPT            = -37, &
                        DFE_BADFP              = -38
                        
  integer, parameter :: DFE_NOREF              = -39, &
                        DFE_BADDATATYPE        = -40, &
                        DFE_BADMCTYPE          = -41, &
                        DFE_BADNUMTYPE         = -42, &
                        DFE_BADORDER           = -43, &
                        DFE_ARGS               = -44, &
                        DFE_INTERNAL           = -45, &
                        DFE_DUPDD              = -46, &
                        DFE_CANTMOD            = -47, &
                        DFE_RANGE              = -48, &
                        DFE_BADTABLE           = -49, &
                        DFE_BADSDG             = -50, &
                        DFE_BADNDG             = -51, &
                        DFE_BADFIELDS          = -52, &
                        DFE_NORESET            = -53, &
                        DFE_NOVS               = -54, &
                        DFE_VGSIZE             = -55, &
                        DFE_DIFFFILES          = -56, &
                        DFE_VTAB               = -57, &
                        DFE_BADAID             = -58, &
                        DFE_OPENAID            = -59, &
                        DFE_BADCONV            = -60, &
                        DFE_GENAPP             = -61, &
                        DFE_CANTFLUSH          = -62, &
                        DFE_BADTYPE            = -63, &
                        DFE_SYMSIZE            = -64, &
                        DFE_BADATTACH          = -65, &
                        DFE_CANTDETACH         = -66

! internal file access codes

  integer, parameter ::  DFACC_READ             = 1, &
                         DFACC_WRITE            = 2, &
                         DFACC_CREATE           = 4, &
                         DFACC_ALL              = 7, &
                         DFACC_RDONLY           = 1, &
                         DFACC_RDWR             = 3, &
                         DFACC_CLOBBER          = 4

!       Access types for SDsetaccesstype

  integer, parameter ::  DFACC_DEFAULT         = 0, &
                         DFACC_SERIAL          = 1, &
                         DFACC_PARALLEL        = 9

!       Constants for DFSDsetorder

  integer, parameter ::  DFO_FORTRAN            = 1, &
                         DFO_C                  = 2

!       Definitions of storage convention

  integer, parameter :: DFNTF_IEEE             = 1, &
                        DFNTF_VAX              = 2, &
                        DFNTF_CRAY             = 3, &
                        DFNTF_PC               = 4, &
                        DFNTF_CONVEX           = 5, &
                        DFNTF_VP               = 6

!       Masks for types

  integer, parameter ::  DFNT_HDF               = 0, &
                         DFNT_NATIVE            = 4096, &
                         DFNT_CUSTOM            = 8192, &
                         DFNT_LITEND            = 16384

!       Number type info codes 

  integer, parameter :: DFNT_NONE     = 0, &
                        DFNT_QUERY    = 0, &
                        DFNT_VERSION  = 1

  integer, parameter :: DFNT_FLOAT32  = 5, &
                        DFNT_FLOAT    = 5, &
                        DFNT_FLOAT64  = 6, &
                        DFNT_DOUBLE   = 6, &
                        DFNT_FLOAT128 = 7, &
              
                        DFNT_INT8     = 20, &
                        DFNT_UINT8    = 21, &
                        DFNT_INT16    = 22, &
                        DFNT_UINT16   = 23, &
                        DFNT_INT32    = 24, &
                        DFNT_UINT32   = 25, &
                        DFNT_INT64    = 26, &
                        DFNT_UINT64   = 27, &
                        DFNT_INT128   = 28, &
                        DFNT_UINT128  = 29, &
              
                        DFNT_UCHAR8   = 3, &
                        DFNT_UCHAR    = 3, &
                        DFNT_CHAR8    = 4, &
                        DFNT_CHAR     = 4, &
                        DFNT_CHAR16   = 42, &
                        DFNT_UCHAR16  = 43

  integer, parameter :: DFNT_NFLOAT32 = 4101, &
                        DFNT_NFLOAT   = 4101, &
                        DFNT_NFLOAT64 = 4102, &
                        DFNT_NDOUBLE  = 4102, &
                        DFNT_NFLOAT128= 4103, &
              
                        DFNT_NINT8    = 4116, &
                        DFNT_NUINT8   = 4117, &
                        DFNT_NINT16   = 4118, &
                        DFNT_NUINT16  = 4119, &
                        DFNT_NINT32   = 4120, &
                        DFNT_NUINT32  = 4121, &
                        DFNT_NINT64   = 4122, &
                        DFNT_NUINT64  = 4123, &
                        DFNT_NINT128  = 4124, &
                        DFNT_NUINT128 = 4125, &
              
                        DFNT_NUCHAR8  = 4099, &
                        DFNT_NUCHAR   = 4099, &
                        DFNT_NCHAR8   = 4100, &
                        DFNT_NCHAR    = 4100, &
                        DFNT_NCHAR16  = 4138, &
                        DFNT_NUCHAR16 = 4139

  integer, parameter :: DFNT_LFLOAT32 = 16389, &
                        DFNT_LFLOAT   = 16389, &
                        DFNT_LFLOAT64 = 16390, &
                        DFNT_LDOUBLE  = 16390, &
                        DFNT_LFLOAT128= 16391, &
              
                        DFNT_LINT8    = 16404, &
                        DFNT_LUINT8   = 16405, &
                        DFNT_LINT16   = 16406, &
                        DFNT_LUINT16  = 16407, &
                        DFNT_LINT32   = 16408, &
                        DFNT_LUINT32  = 16409, &
                        DFNT_LINT64   = 16410, &
                        DFNT_LUINT64  = 16411, &
                        DFNT_LINT128  = 16412, &
                        DFNT_LUINT128 = 16413, &
              
                        DFNT_LUCHAR8  = 16387, &
                        DFNT_LUCHAR   = 16387, &
                        DFNT_LCHAR8   = 16388, &
                        DFNT_LCHAR    = 16388, &
                        DFNT_LCHAR16  = 16426, &
                        DFNT_LUCHAR16 = 16427

!       tags and refs

  integer, parameter :: DFREF_WILDCARD   = 0, &
                        DFTAG_WILDCARD   = 0, &
                        DFTAG_NULL       = 1, &
                        DFTAG_LINKED     = 20, &
                        DFTAG_VERSION    = 30, &
                        DFTAG_COMPRESSED = 40


!       utility set

  integer, parameter :: DFTAG_FID     = 100, &
                        DFTAG_FD      = 101, &
                        DFTAG_TID     = 102, &
                        DFTAG_TD      = 103, &
                        DFTAG_DIL     = 104, &
                        DFTAG_DIA     = 105, &
                        DFTAG_NT      = 106, &
                        DFTAG_MT      = 107

!       raster-8 set 

  integer, parameter :: DFTAG_ID8     = 200, &
                        DFTAG_IP8     = 201, &
                        DFTAG_RI8     = 202, &
                        DFTAG_CI8     = 203, &
                        DFTAG_II8     = 204

!       Raster Image set

  integer, parameter :: DFTAG_ID      = 300, &
                        DFTAG_LUT     = 301, &
                        DFTAG_RI      = 302, &
                        DFTAG_CI      = 303

  integer, parameter :: DFTAG_RIG     = 306, &
                        DFTAG_LD      = 307, &
                        DFTAG_MD      = 308, &
                        DFTAG_MA      = 309, &
                        DFTAG_CCN     = 310, &
                        DFTAG_CFM     = 311, &
                        DFTAG_AR      = 312

  integer, parameter :: DFTAG_DRAW  = 400,  &
                        DFTAG_RUN   = 401,  &
                        DFTAG_XYP   = 500,  &
                        DFTAG_MTO   = 501

!       Tektronix 

  integer, parameter :: DFTAG_T14   = 602,  &
                        DFTAG_T105  = 603

!       Scientific Data set 

  integer, parameter :: DFTAG_SDG   = 700,  &
                        DFTAG_SDD   = 701,  &
                        DFTAG_SD    = 702,  &
                        DFTAG_SDS   = 703,  &
                        DFTAG_SDL   = 704,  &
                        DFTAG_SDU   = 705,  &
                        DFTAG_SDF   = 706,  &
                        DFTAG_SDM   = 707,  &
                        DFTAG_SDC   = 708,  &
                        DFTAG_SDT   = 709,  &
                        DFTAG_SDLNK = 710,  &
                        DFTAG_NDG   = 720,  &
                        DFTAG_CAL   = 731,  &
                        DFTAG_FV    = 732,  &
                        DFTAG_BREQ  = 799,  &
                        DFTAG_EREQ  = 780

!       VSets 

  integer, parameter :: DFTAG_VG    = 1965, &
                        DFTAG_VH    = 1962, &
                        DFTAG_VS    = 1963

!       compression schemes 

  integer, parameter :: DFTAG_RLE           =11, &
                        DFTAG_IMC           =12, &
                        DFTAG_IMCOMP        =12, &
                        DFTAG_JPEG          =13, &
                        DFTAG_GREYJPEG      =14

!       SPECIAL CODES 

  integer, parameter :: SPECIAL_LINKED      = 1, &
                        SPECIAL_EXT         = 2

!       PARAMETERS 

  integer, parameter :: DF_MAXFNLEN         = 256, &
                        SD_UNLIMITED        = 0,   &
                        SD_DIMVAL_BW_COMP   = 1,   &
                        SD_DIMVAL_BW_INCOMP = 0


!       Standard return codes       
 
  integer, parameter :: SUCCEED = 0,        &
                        FAIL    = -1


!       Compression Types 

  integer, parameter :: COMP_NONE   = 0,    &
                        COMP_RLE    = 11,   &
                        COMP_IMCOMP = 12,   &
                        COMP_JPEG   = 2

!       Interlace Types 

  integer, parameter :: MFGR_INTERLACE_PIXEL        = 0, &
                        INTERLACE_LINE              = 1, &
                        MFGR_INTERLACE_COMPONENT    = 2

!       Vdata fields packing types
  integer, parameter :: HDF_VSPACK   = 0, &
                        HDF_VSUNPACK = 1
                        
  ! Constants for the netCDF interface from the hlimits.h file. According to the HDF help
  !   desk, these are also the constants to use in the SD API. 
  integer, parameter :: MAX_NC_DIMS  = 5000, &
                        MAX_NC_ATTRS = 3000, &
                        MAX_NC_VARS  = 5000, &
                        MAX_NC_NAME  = 256, &
                        MAX_VAR_DIMS = 32

  ! Funtions in the SDAPI which use "type punning" - i.e. they can take
  !   arguments of arbitrary type and rank. Fortran 90 has no way
  !   to deal with this, so all we can do is declare the functions as external. 
  integer, external ::  SFrdata, SFrcdata, SFrnatt, SFwdata, SFsnatt, SFsdtstr 

  ! Interface definitions for selected procedures in the Scientific Dataset
  !  API (SDAPI).  By adding interface definitions, calling procedures 
  !  which USE this module can be checked to see that the type, kind, and rank
  !  of the arguments being passed to these procedures is consistent with 
  !  what the procedures expect. 
  interface   
    integer function SFstart(filename, access_mode)
      character(len = *),               intent( in) :: filename
      integer (kind = selected_int_kind(9)), intent( in) :: access_mode
    end function SFstart

    integer function SFend(sd_id)
      integer (kind = selected_int_kind(9)), intent( in) :: sd_id
    end function SFend

    integer function SFselect(sd_id, sds_index)
      integer (kind = selected_int_kind(9)), intent( in) :: sd_id, sds_index
    end function SFselect
    
    integer function SFcreate(sd_id, sds_name, number_type, rank, dimsizes)
      integer (kind = selected_int_kind(9)), intent( in) :: sd_id
      character(len = *),               intent( in) :: sds_name
      integer (kind = selected_int_kind(9)), intent( in) :: number_type, rank
      integer (kind = selected_int_kind(9)), &
                          dimension(*), intent( in) :: dimsizes
    end function SFcreate

    integer function SFendacc(sds_id)
      integer (kind = selected_int_kind(9)), intent( in) :: sds_id
    end function SFendacc

    integer function SFfinfo(sd_id, ndatasets, nglobal_attr)
      integer (kind = selected_int_kind(9)), intent( in) :: sd_id
      integer (kind = selected_int_kind(9)), intent(out) :: ndatasets, nglobal_attr
    end function SFfinfo

    integer function SFn2index(sd_id, sds_name)
      integer (kind = selected_int_kind(9)), intent( in) :: sd_id
      character (len = *),              intent( in) :: sds_name
    end function SFn2index

    integer function SFginfo(sds_id, sds_name, rank, dimsizes, number_type, nattrs)
      integer (kind = selected_int_kind(9)), intent( in) :: sds_id
      character(len = *),               intent(out) :: sds_name
      integer (kind = selected_int_kind(9)), intent(out) :: rank, number_type, nattrs
      integer (kind = selected_int_kind(9)), &
                          dimension(*), intent(out) :: dimsizes
    end function SFginfo

    integer function SFgdinfo(dim_id, name, count, data_type, nattrs)
      integer (kind = selected_int_kind(9)), intent( in) :: dim_id
      character (len = *),              intent(out) :: name
      integer (kind = selected_int_kind(9)), intent(out) :: count, data_type, nattrs
    end function SFgdinfo

    integer function SFgainfo(s_id, attr_index, attr_name, number_type, count)
      integer (kind = selected_int_kind(9)), intent( in) :: s_id, attr_index
      character(len = *),               intent(out) :: attr_name
      integer (kind = selected_int_kind(9)), intent(out) :: number_type, count
    end function SFgainfo

    integer function SFfattr(s_id, attr_name)
      integer (kind = selected_int_kind(9)), intent( in) :: s_id
      character(len = *),               intent( in) :: attr_name
    end function SFfattr

    integer function SFrcatt(s_id, attr_index, value)
      integer (kind = selected_int_kind(9)), intent( in) :: s_id, attr_index
      character(len = *),               intent(out) :: value
    end function SFrcatt
    
    integer function SFscatt(sd_id, attr_name, data_type, count, value)
      integer (kind = selected_int_kind(9)), intent( in) :: sd_id
      character (len = *),              intent( in) :: attr_name
      integer (kind = selected_int_kind(9)), intent( in) :: data_type, count
      !modified this definition of value.  EGMoody 11/16/1999
      character(len = *),               intent( in) :: value
      !Orig Def: character(len = *), dimension(*), intent( in) :: value
    end function SFscatt
    
    integer function SFdimid(sds_id, dim_number)
      integer (kind = selected_int_kind(9)), intent( in) :: sds_id, dim_number
    end function SFdimid

    integer function SFsdmname(dim_id, dim_name)
      integer (kind = selected_int_kind(9)), intent( in) :: dim_id
      character(len = *),               intent( in) :: dim_name
    end function SFsdmname

  end interface
end module hdf_mod
