#ifndef __PROTOTYPE4_FEM_H__
#define __PROTOTYPE4_FEM_H__

#include <map>
#include <string>
#include <vector>

namespace PROTOTYPE4_FEM
{
/*! Packet ID */
const int PACKET_ID = 21351;

/*! Number of ADC Samples per tower */
const int NSAMPLES = 31;

//! Mask to obtain ADC from DWord data
const int ADC_DATA_MASK = (1 << 14) - 1;

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
int GetChannelNumber(const std::string &caloname, int i_column, int i_row);

//! Abhisek's power-law + exp fit
bool SampleFit_PowerLawExp(              //
    const std::vector<double> &samples,  //
    double &peak,                        //
    double &peak_sample,                 //
    double &pedstal,                     //
    const int verbosity = 0);

//! Power law double exp fit
bool SampleFit_PowerLawDoubleExp(          //
    const std::vector<double> &samples,    //
    double &peak,                          //! peak amplitude.
    double &peak_sample,                   //! peak sample position. Fixed to the input value if NOT NAN
    double &pedestal,                      //! pedestal
    std::map<int, double> &parameters_io,  //! IO for fullset of parameters. If a parameter exist and not an NAN, the fit parameter will be fixed to that value. The order of the parameters are ("Amplitude 1", "Sample Start", "Power", "Peak Time 1", "Pedestal", "Amplitude 2", "Peak Time 2")
    const int verbosity = 0);

//! Just return the max sample...
bool SampleFit_PeakSample(          //
    const std::vector<double> &samples,    //
    double &peak,                          //! peak amplitude.
    double &peak_sample,                   //! peak sample position. Fixed to the input value if NOT NAN
    double &pedestal,                      //! pedestal
    const int verbosity = 0);


// Abhisek's power-law + exp signal shape model
double
SignalShape_PowerLawExp(double *x, double *par);
double
SignalShape_PowerLawDoubleExp(double *x, double *par);

//! special treatment for EMCal tagging packet
//! See also https://wiki.bnl.gov/sPHENIX/index.php/2017_calorimeter_beam_test#What_is_new_in_the_data_structures_in_2017
//! Result stored in RUN_INFO node under variable EMCAL_Is_HighEta
const int PACKET_EMCAL_HIGHETA_FLAG = 905;
}

#endif
