#ifndef __PACKETCONSTANTS_H__
#define __PACKETCONSTANTS_H__


/* Misc. values  */
#define MAX_OUTLENGTH 40000

// define some offset which takes us out well > 30,000 for our ID;s
#define IDOFFSET    30000

// normal pass through mode
#define IDDCM0OFFSET 400  
// normal fpga zero suppression mode
#define IDDCM1OFFSET 500  
// extra
#define IDDCM2OFFSET 600  
// level1 packets
#define IDLL1OFFSET  700  
#define IDL2OFFSET   750

// alternate (long) format pass through mode
#define IDDCM3OFFSET 800  
// short format (emcal, etc) 
#define IDDCMSOFFSET 900  

// ---------------------------------------------------------------------
//    IDCRAW requests the subevent to be copied without any decoding
#define IDCRAW   IDOFFSET + 0 

// ---------------------------------------------------------------------
//    IDDGEN uses the standard decoding method imbedded in the subevent
//    header in the new data format
#define IDDGEN    IDOFFSET + 1 

// ---------------------------------------------------------------------
//    IDHCPY requests only the subevent header (or the Event header) to be
//    copied:
#define IDHCPY    IDOFFSET + 2 

// ---------------------------------------------------------------------
//    the next methods < 10 use what we consider standard methods by
//    now, i.e., no scheme proprietary to one particular hardware brand

#define ID1STR    IDOFFSET + 3 
#define IDCSTR    IDOFFSET + 4 
#define ID2EVT    IDOFFSET + 5 
#define ID4EVT    IDOFFSET + 6 
#define ID2SUP    IDOFFSET + 7 
#define ID4SCALER IDOFFSET + 8 

// ---------------------------------------------------------------------
// the next methods are for the hammond/g-2 board.

#define IDHAMMONDSET    IDOFFSET + 31
#define IDHAMMOND       IDOFFSET + 32

#define IDSAM           IDOFFSET + 40

#define IDMIZNHC        IDOFFSET + 41

#define IDDCFEM         IDOFFSET + 51


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
#define IDZDC_DCM0    IDDCM0OFFSET + 13
#define IDPXL_DCM0    IDDCM0OFFSET + 24

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
#define IDZDC_DCM1    IDDCM1OFFSET + 13

// the "level 2", data further compressed by the DSP

#define IDBBC_DCM2    IDDCM2OFFSET + 1
#define IDMVD_DCM2    IDDCM2OFFSET + 2
#define IDDCH_DCM2    IDDCM2OFFSET + 3
#define IDPC_DCM2     IDDCM2OFFSET + 4
#define IDTEC_DCM2    IDDCM2OFFSET + 5
#define IDRICH_DCM2   IDDCM2OFFSET + 6
//#define IDTOF_Q1Q2T3T4    IDDCM2OFFSET + 7
#define IDTOF_DCM2    IDDCM2OFFSET + 7
#define IDEMC_OLDSTYLE  IDDCM2OFFSET + 58
#define IDPBGL_DCM2   IDDCM2OFFSET + 9
#define IDMUTA_DCM2   IDDCM2OFFSET + 10
#define IDMUTC_DCM2   IDDCM2OFFSET + 11
//#define IDMUID_DCM2   IDDCM2OFFSET + 12
#define IDZDC_DCM2    IDDCM2OFFSET + 13

// the "level 3", alternate (long) format in pass through mode

#define IDBBC_DCM3    IDDCM3OFFSET + 1
#define IDMVD_DCM3    IDDCM3OFFSET + 2
#define IDDCH_DCM3    IDDCM3OFFSET + 3
// moved to idpc_fpga #define IDPC_DCM3     IDDCM3OFFSET + 4
#define IDTEC_DCM3    IDDCM3OFFSET + 5
#define IDRICH_DCM3   IDDCM3OFFSET + 6
#define IDTOF_DCM3    IDDCM3OFFSET + 7
// emc FEM to DCM long format (192 channels, user words, ...)
#define IDPBSC_DCM3   IDDCM3OFFSET + 8
#define IDPBGL_DCM3   IDDCM3OFFSET + 9


// the emc short formats

#define IDPBSC_DCMS   908
#define IDPBGL_DCMS   909

// the emc zero-suppressed short formats (3 words per channel+address)

#define IDPBSC_DCMZS   608
#define IDPBGL_DCMZS   609


// the "pbsc 32 channel format" 

#define IDEMC_DCM32  808
#define IDPBGL_DCM32  809

// the emc non-suppressed format from the DCM (144 channels, no user words,...)
#define IDPBSC_DCM5  1008
#define IDPBGL_DCM5  1009

// the emc zero-suppressed format from the DCM (5 words per channel+address) 
#define IDPBSC_DCM05  1108
#define IDPBGL_DCM05  1109

// the fcal zero-suppressed formats (it will use the emcs 1008 1108 packets)
#define IDFCAL_FPGA       1016
#define IDFCAL_FPGA0SUP   1216
#define IDFCAL_FPGA3      1316
#define IDFCAL_FPGA0SUP3  1116


#define IDTOF_DCM16  307

// IDDCM3OFFSET = 800
#define IDMUTA_DCM3   IDDCM3OFFSET + 10
#define IDMUTC_DCM3   IDDCM3OFFSET + 11
#define IDMUID_DCM3   IDDCM3OFFSET + 12
#define IDZDC_DCM3    IDDCM3OFFSET + 13

#define IDFOCAL_FPGATEST 725

#define IDMUTRG_DCM0 791


// we start two new series -- 1000 : through fpga but not zero-supressed
//                         -- 1100 : through fpga AND zero-supressed




#define IDBBC_FPGA         1001
#define IDBBC_FPGA0SUP     1101

#define IDMVD_FPGA       1002
#define IDMVD_FPGA0SUP   1102

#define IDMVD_PED_FPGA0SUP	1502

#define IDPC_FPGA         804
#define IDPC_FPGA0SUP     1104

#define IDRICH_FPGA       1006
#define IDRICH_FPGA0SUP   1106

#define IDTOF_FPGA       1007
#define IDTOF_FPGA0SUP   1107

#define IDTOFW_FPGA       1057
#define IDTOFW_FPGA0SUP   1157

#define IDEMC_FPGA       1008
#define IDEMC_FPGA0SUP   1108

#define IDEMC_FPGASHORT       1208
#define IDEMC_FPGASHORT0SUP   1308

#define IDEMC_FPGA3WORDS     1408
#define IDEMC_FPGA3WORDS0SUP 1508

#define IDEMC_REFERENCE   1058
#define IDEMC_REFERENCE0SUP   1158

#define IDEMC_SHORTREFERENCE   1068
#define IDEMC_SHORTREFERENCE0SUP   1168

#define IDMUTC_FPGA          1011
#define IDMUTC_FPGA0SUP      1111
#define IDMUTC_FPGASHORT     1211
#define IDMUTC_FPGASHORTSUP  1311
#define IDMUTC_FPGANEW       1411
#define IDMUTC_FPGANEWSUP    1511

#define IDMUTC_15_FPGA       1051
#define IDMUTC_15_FPGA0SUP   1151

#define IDMUID_FPGA       1012
#define IDMUID_FPGA0SUP   1112

#define IDZDC_FPGA       1013
#define IDZDC_FPGA0SUP   1113

#define IDNTCT0_FPGA       1015
#define IDNTCT0_FPGA0SUP   1115

#define IDRPC_DCM0       1019
#define IDRPC_FPGA       1219
#define IDRPC_FPGA0SUP   1319


// HBD gets number 22
#define IDHBD_FPGA       1022
#define IDHBD_FPGA0SUP   1122
#define IDHBD_FPGASHORT       1222
#define IDHBD_FPGASHORT0SUP   1322
#define IDHBD_FPGA3SAMPLES    1422
#define IDHBD_FPGA3SAMPLES0SUP    1522

// RXNP  gets 23

#define  IDRXNP_FPGASHORT 1323
#define  IDRXNP_FPGASHORT0SUP 1423

// the "LL1", level 1 trigger info

#define IDBBC_LL1    IDLL1OFFSET + 1
#define IDMVD_LL1    IDLL1OFFSET + 2
#define IDRICH_LL1   IDLL1OFFSET + 6
#define IDTOF_LL1    IDLL1OFFSET + 7
#define IDPBSC_LL1   IDLL1OFFSET + 8
#define IDPBGL_LL1   IDLL1OFFSET + 9
#define IDMUIDH_LL1   IDLL1OFFSET + 12
#define IDMUIDV_LL1   IDLL1OFFSET + 13
#define IDGL1        IDLL1OFFSET + 14
#define IDGL1P       IDDCM3OFFSET + 14
#define IDGL1PSUM    914
#define IDGL1PSUMOBS 818
#define IDEMCRICH_LL1 IDLL1OFFSET + 15
#define IDNTCZDC_LL1  IDLL1OFFSET + 16
#define IDGL1_EVCLOCK   IDLL1OFFSET + 17
#define IDERT_E_LL1    IDLL1OFFSET + 18
#define IDERT_W_LL1    IDLL1OFFSET + 19
#define IDBIG_LL1      IDLL1OFFSET + 77

// L2 packets
//
#define IDL2DECISION IDL2OFFSET
#define IDL2PRIMITIVE IDL2OFFSET + 1


// the CDEV data formats, starting from 2000

#define IDCDEVIR            2001
#define IDCDEVDVM           2002
#define IDCDEVRING          2003
#define IDCDEVWCMHISTORY    2004
#define IDCDEVSIS           2005
#define IDCDEVPOLARIMETER   2006
#define IDCDEVPOLDATA       2007
#define IDCDEVPOLARIMETERTARGET 2008
#define IDCDEVBPM           2009
#define IDCDEVMADCH         2010
#define IDGL1PSCALER        2011
#define IDCDEVRINGPOL       2012
#define IDCDEVRINGFILL      2013
#define IDCDEVBUCKETS       2014
#define IDCDEVRINGNOPOL     2015
#define IDCDEVPOLARIMETERZ  2016
#define IDCDEVDESCR         2017
#define IDSTARSCALER        2098


// EMC data header and trailer length

#define EMC_SUPPRESSED_DATA_HEADER_LENGTH 8
#define EMC_DATA_TRAILER_LENGTH 10
#define EMC_SHORT_DATA_HEADER_LENGTH 9
#define EMC_LONG_DATA_HEADER_LENGTH 9
#define EMC_WORDS_PER_CH_SHORT 3 
#define EMC_WORDS_PER_CH_LONG 5 
#define EMC_DCMDATA_HEADER_LENGTH 8
#define EMC_DCMDATA_TRAILER_LENGTH 2
#endif /* __PACKETCONSTANTS_H__ */



