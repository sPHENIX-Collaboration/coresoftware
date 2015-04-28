#ifndef __SUBEVT_CONSTANTS_H
#define __SUBEVT_CONSTANTS_H

/* the enum types for dump style */
#define EVT_DECIMAL     1
#define EVT_HEXADECIMAL 2
#define EVT_OCTAL       3

/* Misc. values  */
#define MAX_OUTLENGTH 80000

// the header length value
#define SEVTHEADERLENGTH 4



// ---------------------------------------------------------------------
//    IDCRAW requests the subevent to be copied without any decoding
#define IDCRAW   0 

// ---------------------------------------------------------------------
//    IDDGEN uses the standard decoding method imbedded in the subevent
//    header in the new data format
#define IDDGEN    1 

// ---------------------------------------------------------------------
//    IDHCPY requests only the subevent header (or the Event header) to be
//    copied:
#define IDHCPY    2 

// ---------------------------------------------------------------------
//    the next methods < 10 use what we consider standard methods by
//    now, i.e., no scheme proprietary to one particular hardware brand

#define ID1STR    3 
#define IDCSTR    4 
#define ID2EVT    5 
#define ID4EVT    6 
#define ID2SUP    7 

// ---------------------------------------------------------------------
// the next methods are for the hammond/g-2 board.

#define IDHAMMONDSET    31
#define IDHAMMOND       32

#define IDSAM           40

#define IDMIZNHC        41

#define IDDCFEM         51
#define IDTECFEM        52

#define IDSIS3300       55
#define IDCAENV792      56
#define IDCAENV785N     57

#define IDFIFOBOARD     58
#define IDRCPETDATA     59
#define IDBSPETDATA     60
#define IDUPPETDATA     61
#define IDUPPETDATA_V104     62
#define IDSIS3300R       65


#define IDSRSV01          70
#define IDUPPETPARAMS     71

#define IDDRS4V1        81
#define IDCAENV1742     85

#define IDPCONTAINER     89

#define IDFNALMWPC     90
#define IDFNALMWPCV2     91



// the "level 0", meaning the raw untreated FEM data 

#define IDBBC_DCM0    IDDCM0OFFSET + 1
#define IDMVD_DCM0    IDDCM0OFFSET + 2
#define IDDCH_DCM0    IDDCM0OFFSET + 3
#define IDPC_DCM0     IDDCM0OFFSET + 4
#define IDTEC_DCM0    IDDCM0OFFSET + 5
#define IDRICH_DCM0   IDDCM0OFFSET + 6
#define IDTOF_DCM0    IDDCM0OFFSET + 7
#define IDPBSC_DCM0   IDDCM0OFFSET + 8
#define IDPBGL_DCM0   IDDCM0OFFSET + 9
#define IDMUTA_DCM0   IDDCM0OFFSET + 10
#define IDMUTC_DCM0   IDDCM0OFFSET + 11
#define IDMUID_DCM0   IDDCM0OFFSET + 12

// the "level 1", FEM data zero-suppressed by the FPGA

#define IDBBC_DCM1    IDDCM1OFFSET + 1
#define IDMVD_DCM1    IDDCM1OFFSET + 2
#define IDDCH_DCM1    IDDCM1OFFSET + 3
#define IDPC_DCM1     IDDCM1OFFSET + 4
#define IDTEC_DCM1    IDDCM1OFFSET + 5
#define IDRICH_DCM1   IDDCM1OFFSET + 6
#define IDTOF_DCM1    IDDCM1OFFSET + 7
#define IDPBSC_DCM1   IDDCM1OFFSET + 8
#define IDPBGL_DCM1   IDDCM1OFFSET + 9
#define IDMUTA_DCM1   IDDCM1OFFSET + 10
#define IDMUTC_DCM1   IDDCM1OFFSET + 11
#define IDMUID_DCM1   IDDCM1OFFSET + 12

// the "level 2", data further compressed by the DSP

#define IDBBC_DCM2    IDDCM2OFFSET + 1
#define IDMVD_DCM2    IDDCM2OFFSET + 2
#define IDDCH_DCM2    IDDCM2OFFSET + 3
#define IDPC_DCM2     IDDCM2OFFSET + 4
#define IDTEC_DCM2    IDDCM2OFFSET + 5
#define IDRICH_DCM2   IDDCM2OFFSET + 6
#define IDTOF_DCM2    IDDCM2OFFSET + 7
#define IDPBSC_DCM2   IDDCM2OFFSET + 8
#define IDPBGL_DCM2   IDDCM2OFFSET + 9
#define IDMUTA_DCM2   IDDCM2OFFSET + 10
#define IDMUTC_DCM2   IDDCM2OFFSET + 11
#define IDMUID_DCM2   IDDCM2OFFSET + 12



#endif
