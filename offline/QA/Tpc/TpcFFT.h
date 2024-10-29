// include guard to prevent TpcFFT.h from being preprocessed multiple times
#ifndef QA_TPC_TPCFFT_H
#define QA_TPC_TPCFFT_H
//

// includes
#include <fun4all/SubsysReco.h>
#include <TMath.h>

#include <trackbase/TrkrDefs.h>

#include <ffarawobjects/TpcRawHitContainer.h>
#include <ffarawobjects/TpcRawHit.h>

#include <string>
#include <vector>
//

// Call classes to be used in code
class PHCompositeNode;
class TFile;
class TH1;
class Fun4AllHistoManager;
//

class TpcFFT : public SubsysReco  // Inherit public parts of SubsysReco
{
  // list of public methods
 public:
  // Function that sets file name and allocated memory to adcSamples
  explicit TpcFFT(const std::string &name = "TpcFFT");

  ~TpcFFT() override = default;

  // Initialization of run: where you fetch data and place histograms that
  // care about run number
  int InitRun(PHCompositeNode *topNode) override;

  // Where you do the desired work associated with the run (analysis)
  int process_event(PHCompositeNode *topNode) override;

  // called at the end of the run when processing is over
  int End(PHCompositeNode *topNode) override;

  // List of private members
 private:

  TH1 *h_WF{nullptr};
  TH1 *h_FFT{nullptr};

  std::vector<TH1*> WF_clone{};
  std::vector<TH1*> FFT_clone{};
  std::vector<Int_t> evt_num{};
  std::vector<Float_t> pedestal{};
  std::vector<Float_t> pedestal_sigma{};

  Float_t pedestal_sum{0.0};     // average ADC value in channel
  Float_t pedestal_sigma_sum{0.0};     // average RMS in channel
  Float_t pedestal_samples{0.0};  // number of times channel appears

  bool dead{false};
  
  int m_nWaveformInFrame{0};  // Number of waveforms in a frame
  int m_Channel{0};           // Channel number
  int m_nSamples{0};          // Number of samples in waveform
  int m_FEE{0};

  int sector{0};
  std::string sec[24]{"00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23"};
  int ent_num{0};

  TFile *TpcFFTfile;

  std::vector<unsigned short> m_adcSamples;  // ADC values in waveform
};

#endif
