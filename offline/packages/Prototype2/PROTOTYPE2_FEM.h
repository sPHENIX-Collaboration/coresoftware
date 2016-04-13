#ifndef __PROTOTYPE2_FEM_H__
#define __PROTOTYPE2_FEM_H__

#include <string>
#include <vector>

namespace PROTOTYPE2_FEM
{

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

  //! FEM mapping of channel -> calorimeter col and rows
  int
  GetHBDCh(std::string caloname, int i_column, int i_row);

  //! Abhisek's power-law + exp fit
  bool
  SampleFit_PowerLawExp(//
      const std::vector<double> & samples, //
      double & peak,//
      double & peak_sample,//
      double & pedstal, //
      const int verbosity = 0
      );

  // Abhisek's power-law + exp signal shape model
  double
  SignalShape_PowerLawExp(double *x, double *par);

}

#endif
