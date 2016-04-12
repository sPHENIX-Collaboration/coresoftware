#ifndef __FEM_H__
#define __FEM_H__

namespace PROTOTYPE2_FEM {

 /*! Packet ID */
 const int PACKET_ID = 21101;

 /*! Number of ADC Samples per tower */
 const int NSAMPLES = 24;

 /*! Number of Inner HCAL towers */
 const int NCH_IHCAL_ROWS = 4;
 const int NCH_IHCAL_COLUMNS = 4;

 /*! Number of Outer HCAL towers */
 const int NCH_OHCAL_ROWS = 4;
 const int NCH_OHCAL_COLUMNS = 4;

 /*! Number of EMCAL towers */
 const int NCH_EMCAL_ROWS = 8;
 const int NCH_EMCAL_COLUMNS = 8;

 /*! Error assigned to saturated ADCs */
 const int SATURATED_ADC_ERROR = 100;

 /*! Error assigned to dead channels */
 const int DEAD_CHANNEL_ERROR = 300;

}

#endif
