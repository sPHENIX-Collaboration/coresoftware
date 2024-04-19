#ifndef MBD_MBDEVENT_H
#define MBD_MBDEVENT_H

#include "MbdDefs.h"
#include "MbdSig.h"

#include <TFile.h>
#include <TTree.h>
#ifndef ONLINE
#include <fun4all/Fun4AllBase.h>
#endif
#include <vector>

class PHCompositeNode;
class Event;
class Packet;
class MbdPmtContainer;
class MbdOut;
class MbdCalib;
class MbdGeom;
class CDBUtils;
class TF1;
class TCanvas;

class MbdEvent
{
 public:
  MbdEvent(const int cal_pass = 0);
  virtual ~MbdEvent();

  int SetRawData(Event *event, MbdPmtContainer *mbdpmts);
  int Calculate(MbdPmtContainer *mbdpmts, MbdOut *mbdout);
  int InitRun();
  int End();
  void Clear();

  void SetSim(const int s) { _simflag = s; }

  float get_bbcz() { return m_bbcz; }
  float get_bbczerr() { return m_bbczerr; }
  float get_bbct0() { return m_bbct0; }
  float get_bbct0err() { return m_bbct0err; }

  int get_bbcn(const int iarm) { return m_bbcn[iarm]; }
  float get_bbcq(const int iarm) { return m_bbcq[iarm]; }
  float get_bbct(const int iarm) { return m_bbct[iarm]; }
  float get_bbcte(const int iarm) { return m_bbcte[iarm]; }

  int get_pmtq(const int ipmt) { return m_pmtq[ipmt]; }
  float get_pmttt(const int ipmt) { return m_pmttt[ipmt]; }
  float get_pmttq(const int ipmt) { return m_pmttq[ipmt]; }

  int get_EventNumber(void) const { return m_evt; }

  void set_debugintt(const int d) { _debugintt = d; }

  MbdSig *GetSig(const int ipmt) { return &_mbdsig[ipmt]; }

  int FillSampMaxCalib();

  int  calib_is_done() { return _calib_done; }
  int  Verbosity() { return _verbose; }
  void Verbosity(const int v) { _verbose = v; }

 private:
  static const int NCHPERPKT = 128;

  MbdGeom *_mbdgeom{nullptr};
  MbdCalib *_mbdcal{nullptr};

  int Read_Charge_Calib(const std::string &calfname);
  int Read_TQ_T0_Offsets(const std::string &calfname);
  int Read_TQ_CLK_Offsets(const std::string &calfname);
  int Read_TT_CLK_Offsets(const std::string &calfname);
  int DoQuickClockOffsetCalib();

  int _debugintt{0};
  void ReadSyncFile(const char *fname = "SYNC_INTTMBD.root");

  float gaincorr[MbdDefs::MBD_N_PMT]{};       // gain corrections
  float tq_t0_offsets[MbdDefs::MBD_N_PMT]{};  // t0 offsets in charge channels
  float tq_clk_offsets[MbdDefs::MBD_N_PMT]{};
  float tt_clk_offsets[MbdDefs::MBD_N_PMT]{};

  // float bz_offset{0.};

  int _verbose{0};
  int _runnum{0};
  int _simflag{0};
  int _nsamples{31};
  int _calib_done{0};
  int _is_online{0};
  Packet *p[2]{nullptr, nullptr};

  // alignment data
  Int_t m_evt{0};
  Short_t m_clk{0};
  Short_t m_femclk{0};

  // raw data
  Float_t m_adc[MbdDefs::MBD_N_FEECH][MbdDefs::MAX_SAMPLES]{};   // raw waveform, adc values
  Float_t m_samp[MbdDefs::MBD_N_FEECH][MbdDefs::MAX_SAMPLES]{};  // raw waveform, sample values
  Float_t m_ampl[MbdDefs::MBD_N_FEECH]{};                        // raw amplitude

  std::vector<MbdSig> _mbdsig;

  Float_t m_pmtq[MbdDefs::MBD_N_PMT]{};   // npe in each arm
  Float_t m_pmttt[MbdDefs::MBD_N_PMT]{};  // time in each arm
  Float_t m_pmttq[MbdDefs::MBD_N_PMT]{};  // time in each arm

  int do_templatefit{1};

  // output data
  Short_t m_bbcn[2]{};                                            // num hits for each arm (north and south)
  Float_t m_bbcq[2]{};                                            // total charge (currently npe) in each arm
  Float_t m_bbct[2]{};                                            // time in arm
  Float_t m_bbcte[2]{};                                           // earliest hit time in arm
  Float_t m_bbctl[2]{};                                           // latest hit time in arm
  Float_t m_bbcz{std::numeric_limits<Float_t>::quiet_NaN()};      // z-vertex
  Float_t m_bbczerr{std::numeric_limits<Float_t>::quiet_NaN()};   // z-vertex error
  Float_t m_bbct0{std::numeric_limits<Float_t>::quiet_NaN()};     // start time
  Float_t m_bbct0err{std::numeric_limits<Float_t>::quiet_NaN()};  // start time error
  Float_t _tres = std::numeric_limits<Float_t>::quiet_NaN();      // time resolution of one channel

  TH1 *hevt_bbct[2]{};  // time in each bbc, per event
  TF1 *gausfit[2]{nullptr, nullptr};

  float TRIG_SAMP[16]{};  // [board]

  // Calibration Data
  int _calpass{0};
  TString _caldir;
  //std::string _caldir;

  // sampmax
  int CalcSampMaxCalib();
  std::unique_ptr<TFile> _smax_tfile{nullptr};
  TH1 *h_tmax[256]{};     // [feech], max sample in event
  TH2 *h2_tmax[2]{};      // [0 == time ch, 1 == chg ch], max sample in evt vs ch
  TH2 *h2_wave[2]{};      // [0 == time ch, 1 == chg ch], all samples in evt vs ch
  TH2 *h2_trange_raw{};   // raw tdc at maxsamp vs ch
  TH2 *h2_trange{};       // subtracted tdc at maxsamp vs ch
  //TH1 *h_trange[2]{};     // subtracted tdc at maxsamp, [S/N]

  TCanvas *ac{nullptr};  // for plots used during debugging

  // debug stuff
  std::unique_ptr<TFile> _synctfile{nullptr};
  TTree *_syncttree{nullptr};
  int _syncevt{0};
  std::vector<Int_t> bbevt;
  std::vector<UShort_t> bbclk;
  std::vector<Float_t> mybbz;
  std::vector<Long64_t> bco;
  std::vector<Double_t> intz;
  std::vector<Double_t> bbz;
};

#endif /* MBD_MBDEVENT_H */
