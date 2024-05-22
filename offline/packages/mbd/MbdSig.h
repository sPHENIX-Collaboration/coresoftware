#ifndef __MBDSIG_H__
#define __MBDSIG_H__

#include "MbdRunningStats.h"

#include <TH1.h>

#include <fstream>
#include <vector>
#include <queue>


class TFile;
class TTree;
class TGraphErrors;
class TH2;
// class THnSparse;
class MbdCalib;

/**

MbdSig: Single Channel digital signal class, includes processing

*/

class MbdSig
{
 public:
  explicit MbdSig(const int chnum = 0, const int nsamp = 0);
  // explicit MbdSig(const MbdSig &obj);    // never used
  virtual ~MbdSig();

  // MbdSig& operator= (const MbdSig& obj) = delete; // never used

  void SetNSamples( const int s ) { _nsamples = s; }
  void SetY(const Float_t *y, const int invert = 1);
  void SetXY(const Float_t *x, const Float_t *y, const int invert = 1);

  void SetCalib(MbdCalib *mcal);

  TH1 *GetHist() { return hpulse; }
  TGraphErrors *GetGraph() { return gpulse; }
  Double_t GetAmpl() { return f_ampl; }
  Double_t GetTime() { return f_time; }
  Double_t GetIntegral() { return f_integral; }

  /**
   * Fill hists from data between minsamp and maxsamp bins
   * sample number starts from 0
   * Ped0 is the running pedestal and rms
   * */
  void FillPed0(const Int_t minsamp, const Int_t maxsamp);
  void FillPed0(const Double_t minsamp, const Double_t maxsamp);

  Double_t GetPed0() { return ped0; }
  Double_t GetPed0RMS() { return ped0rms; }

  /** Set the global pedestal. Once set, it is applied to the data for all events.  */
  void SetPed0(const Double_t mean, const Double_t rms = 0.);

  /** Use the event by event pedestal, */
  void SetEventPed0Range(const Int_t minsamp, const Int_t maxsamp)
  {
    minped0samp = minsamp;
    maxped0samp = maxsamp;
  }

  void SetEventPed0Range(const Double_t minx, const Double_t maxx)
  {
    minped0x = minx;
    maxped0x = maxx;
  }

  void SetEventPed0PreSamp(const Int_t presample, const Int_t nsamps = -1, const int max_samp = -1);

  void CalcEventPed0(const Int_t minsamp, const Int_t maxsamp);
  void CalcEventPed0(const Double_t minx, const Double_t maxx);
  void CalcEventPed0_PreSamp(const Int_t pre_samp, const Int_t nsamps = 1);

  TH1 *GetPedHist() { return hPed0; }

  /** Leading Edge Discriminator signal */
  Double_t LeadingEdge(const Double_t threshold);  // Leading Edge Discriminator Time

  /** digital CFD, threshold is in fraction of amplitude */
  Double_t dCFD(const Double_t fraction_threshold);

  /** MBD method to get time, max_samp is the sample to use */
  Double_t MBDTDC(const Int_t max_samp);

  /** Get pulse amplitude with spline fit */
  Double_t GetSplineAmpl();

  /** Simple integral to get total charge, etc */
  Double_t Integral(const Double_t xmin, const Double_t xmax);

  /** The maximum value from all samples */
  void LocMax(Double_t &x_at_max, Double_t &ymax, Double_t minxrange = 0., Double_t maxxrange = 0.);

  /** The minimum value from all samples (including negatives) */
  void LocMin(Double_t &x_at_min, Double_t &ymin, Double_t minxrange = 0., Double_t maxxrange = 0.);

  /** Use template fit to get ampl and time */
  Int_t FitTemplate(const Int_t sampmax = -1);
  // Double_t Ampl() { return f_ampl; }
  // Double_t Time() { return f_time; }

  /** Make template waveforms for later fits */
  void SetTemplateSize(const Int_t nptsx, const Int_t nptsy, const Double_t begt, const Double_t endt);
  Int_t SetTemplate(const std::vector<float> &shape, const std::vector<float> &sherr);

  // Double_t FitPulse();
  void SetTimeOffset(const Double_t o) { f_time_offset = o; }
  Double_t TemplateFcn(const Double_t *x, const Double_t *par);
  TF1 *GetTemplateFcn() { return template_fcn; }
  void SetMinMaxFitTime(const Double_t mintime, const Double_t maxtime);

  void WritePedHist();

  void PadUpdate();
  void Print();
  void Verbose(const int v) { _verbose = v; }

 private:
  void Init();

  int _ch;
  int _nsamples;
  int _status{0};

  int _evt_counter{0};
  MbdCalib *_mbdcal{nullptr};

  /** fit values*/
  // should make an array for the different methods
  Double_t f_ampl{0,}; /** best guess (from fit of spline or template, or max adc) */
  Double_t f_time{0.}; /** best guess (from fit of spline or template, or max adc) */

  Double_t f_time_offset{4.0}; /** time offset used in fit */

  Double_t f_integral{0.}; /** integral */

  TH1 *hRawPulse{nullptr};           //!
  TH1 *hSubPulse{nullptr};           //!
  TH1 *hpulse{nullptr};              //!
  TGraphErrors *gRawPulse{nullptr};  //!
  TGraphErrors *gSubPulse{nullptr};  //!
  TGraphErrors *gpulse{nullptr};     //!

  /** for CalcPed0 */
  //std::unique_ptr<MbdRunningStats> ped0stats{nullptr};    //!
  MbdRunningStats *ped0stats{nullptr};    //!
  TH1 *hPed0{nullptr};            //! all events
  TH1 *hPedEvt{nullptr};          //! evt-by-event pedestal
  TF1 *ped_fcn{nullptr};
  TF1 *ped_tail{nullptr};         //! tail of prev signal
  Double_t ped0{0.};                  //!
  Double_t ped0rms{0.};               //!
  int   use_ped0{0};                 //! whether to apply ped0
  Int_t minped0samp{-9999};       //! min sample for event-by-event ped, inclusive
  Int_t maxped0samp{-9999};       //! max sample for event-by-event ped, inclusive
  Double_t minped0x{0.};              //! min x for event-by-event ped, inclusive
  Double_t maxped0x{0.};              //! max x for event-by-event ped, inclusive
  Double_t ped_presamp{};         //! presamples for ped calculation
  Double_t ped_presamp_nsamps{};  //! num of presamples for ped calculation
  Double_t ped_presamp_maxsamp{-1}; //! a peak sample for ped calc (-1 = use max)

  /** for time calibration */
  // Double_t time_calib;

  /** For pulse template extraction */
  TH2 *h2Template{nullptr};
  TH2 *h2Residuals{nullptr};

  TH1 *hAmpl{nullptr};
  TH1 *hTime{nullptr};
  Int_t template_npointsx{0};
  Int_t template_npointsy{0};
  Double_t template_begintime{0.};
  Double_t template_endtime{0.};
  // Double_t template_min_good_amplitude{20.};    //! for template, in original units of waveform data
  // Double_t template_max_good_amplitude{4080.};  //! for template, in original units of waveform data
  // Double_t template_min_xrange{0.};             //! for template, in original units of waveform data
  // Double_t template_max_xrange{0.};             //! for template, in original units of waveform data
  std::vector<float> template_y;
  std::vector<float> template_yrms;
  TF1 *template_fcn{nullptr};
  Double_t fit_min_time{};  //! min time for fit, in original units of waveform data
  Double_t fit_max_time{};  //! max time for fit, in original units of waveform data

  int _verbose{0};
};

#endif  // __MBDSIG_H__
