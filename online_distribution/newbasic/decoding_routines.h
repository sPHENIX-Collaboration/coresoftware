#ifndef __DECODING_ROUTINES_H
#define __DECODING_ROUTINES_H

#include "event_io.h"

#define WINDOWSEXPORT 

int WINDOWSEXPORT decode_id4evt( int *, int[]  ,int ,int ,int* );

int WINDOWSEXPORT decode_id2evt( int *, short[] ,int ,int ,int* );

int WINDOWSEXPORT decode_idhammondset( int *, int[] ,int ,int ,int* );

int WINDOWSEXPORT decode_idhammond( int *, int[] ,int ,int ,int* );

int WINDOWSEXPORT decode_idsam( int *, int[] ,int ,int ,int* );

int WINDOWSEXPORT decode_iddcfem( int *, int[] ,int ,int ,int* );

int WINDOWSEXPORT decode_idtecfem( int *, int[] ,int ,int ,int* );

int WINDOWSEXPORT decode_idmiznhc( int *, int[] ,int ,int ,int* );

// the "DCM" raw data format decoders

int WINDOWSEXPORT decode_bbc_dcm0( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_mvd_dcm0( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_dch_dcm0( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_pc_dcm0( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_tec_dcm0( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_rich_dcm0( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_tof_dcm0( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_pbsc_dcm0( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_pbgl_dcm0( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_muta_dcm0( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_mutc_dcm0( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_muid_dcm0( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_zdc_dcm0( int *, int[] ,int ,int ,int* );

int WINDOWSEXPORT decode_rich_ll1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_mvd_ll1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_bbc_ll1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_ntczdc_ll1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_big_ll1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_tof_ll1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_muid_ll1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_ert_ll1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_pbgl_ll1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_pbsc_ll1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_gl1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_gl1p( int *, int[] ,int ,int ,int* );

int WINDOWSEXPORT decode_bbc_dcm1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_mvd_dcm1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_dch_dcm1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_pc_dcm1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_tec_dcm1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_rich_dcm1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_tof_dcm1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_pbsc_dcm1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_pbgl_dcm1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_muta_dcm1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_mutc_dcm1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_muid_dcm1( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_zdc_dcm1( int *, int[] ,int ,int ,int* );


int WINDOWSEXPORT decode_bbc_dcm2( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_mvd_dcm2( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_dch_dcm2( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_pc_dcm2( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_tec_dcm2( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_rich_dcm2( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_tof_dcm2( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_pbsc_dcm2( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_pbgl_dcm2( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_muta_dcm2( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_mutc_dcm2( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_muid_dcm2( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_zdc_dcm2( int *, int[] ,int ,int ,int* );

int WINDOWSEXPORT decode_pc_dcm3( int *, int[] ,int ,int ,int* );
int WINDOWSEXPORT decode_mutc_dcm3( int *, int[] ,int ,int ,int* );

int WINDOWSEXPORT decode_pbsc_dcm32( int *, int[] ,int ,int ,int* );

int WINDOWSEXPORT decode_pbsc_dcms( int *, int[] ,int ,int ,int* );

int WINDOWSEXPORT decode_mutc_dcm3( int *, int[] ,int ,int ,int* );

// ----------

#endif
