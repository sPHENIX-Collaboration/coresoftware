#ifndef TPC_TPCDAQDEFS_H
#define TPC_TPCDAQDEFS_H

#include <stdint.h>
#include <map>
#include <string>
#include <utility>  // std::pair, std::make_pair
#include <vector>

class TCanvas;
class TPaveText;

namespace TPCDaqDefs
{
//! TPC v1 FEE test stand decoder
namespace FEEv1
{
static const unsigned int kPACKET_ID = 1024;
static const unsigned int kPACKET_LENGTH = 137;
static const unsigned int kN_CHANNELS = 256;
static const unsigned int kSAMPLE_LENGTH = 128;

static const unsigned int kMaxPadX = 50;
static const unsigned int kMaxPadY = 12;
std::pair<int, int> SAMPAChan2PadXY(uint32_t fee_channel);

class SampleFit_PowerLawDoubleExp_PDFMaker
{
public:
  SampleFit_PowerLawDoubleExp_PDFMaker();
  ~SampleFit_PowerLawDoubleExp_PDFMaker();
  void MakeSectionPage(const std::string & title);
private:
  TCanvas * m_canvas;
  TPaveText * m_pavedtext;
};

//! Power law double exp fit
bool SampleFit_PowerLawDoubleExp(          //
    const std::vector<double> &samples,    //
    double &peak,                          //! peak amplitude.
    double &peak_sample,                   //! peak sample position. Fixed to the input value if NOT NAN
    double &pedestal,                      //! pedestal
    std::map<int, double> &parameters_io,  //! IO for fullset of parameters. If a parameter exist and not an NAN, the fit parameter will be fixed to that value. The order of the parameters are ("Amplitude 1", "Sample Start", "Power", "Peak Time 1", "Pedestal", "Amplitude 2", "Peak Time 2")
    const int verbosity = 0);

// Abhisek's power-law + exp signal shape model
double
SignalShape_PowerLawExp(double *x, double *par);
double
SignalShape_PowerLawDoubleExp(double *x, double *par);

}  // namespace FEEv1

}  // namespace TPCDaqDefs

#endif  //TPC_TPCDAQDEFS_H
