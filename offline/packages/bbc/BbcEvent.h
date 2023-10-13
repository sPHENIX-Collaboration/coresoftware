#ifndef __BBCEVENT_H__
#define __BBCEVENT_H__

#include "BbcDefs.h"
#include "BbcSig.h"

class PHCompositeNode;
class Event;
class Packet;
class BbcPmtInfoContainerV1;
class BbcOut;
class BbcCalib;
class BbcGeom;
class CDBUtils;
class TF1;
class TCanvas;

class BbcEvent
{
public:
  BbcEvent ();
  ~BbcEvent ();

  int SetRawData(Event *event, BbcPmtInfoContainerV1 *bbcpmts);
  int Calculate(BbcPmtInfoContainerV1 *bbcpmts, BbcOut *bbcout);
  int InitRun();
  void Clear();

  float get_bbcz() { return m_bbcz; }
  float get_bbczerr() { return m_bbczerr; }
  float get_bbct0() { return m_bbct0; }
  float get_bbct0err() { return m_bbct0err; }

  int   get_bbcn(const int iarm) { return m_bbcn[iarm]; }
  float get_bbcq(const int iarm) { return m_bbcq[iarm]; }
  float get_bbct(const int iarm) { return m_bbct[iarm]; }
  float get_bbcte(const int iarm) { return m_bbcte[iarm]; }

  int   get_pmtq(const int ipmt) { return m_pmtq[ipmt]; }
  float get_pmttt(const int ipmt) { return m_pmttt[ipmt]; }
  float get_pmttq(const int ipmt) { return m_pmttq[ipmt]; }

  int get_EventNumber (void) const
  {
    return _eventnum;
  }

  BbcSig* GetSig(const int ipmt) { return &bbcsig[ipmt]; }


private:
  static const int NCHPERPKT = 128;

  BbcGeom *_bbcgeom {nullptr};
  BbcCalib *_bbccal {nullptr};

  int Read_Charge_Calib( const char *calfname );
  int Read_TQ_T0_Offsets( const char *calfname );
  int Read_TQ_CLK_Offsets( const char *calfname );
  int Read_TT_CLK_Offsets( const char *calfname );
  int DoQuickClockOffsetCalib();

  float gaincorr[BbcDefs::BBC_N_PMT] {};          // gain corrections
  float tq_t0_offsets[BbcDefs::BBC_N_PMT] {};     // t0 offsets in charge channels
  float tq_clk_offsets[BbcDefs::BBC_N_PMT] {};
  float tt_clk_offsets[BbcDefs::BBC_N_PMT] {};

  //float bz_offset{0.};

  int _verbose;
  int _runnum;
  int _eventnum;
  Packet *p[2] {nullptr,nullptr} ;

  //raw data
  Float_t m_adc[BbcDefs::BBC_N_FEECH][BbcDefs::MAX_SAMPLES];     // raw waveform, adc values
  Float_t m_samp[BbcDefs::BBC_N_FEECH][BbcDefs::MAX_SAMPLES];    // raw waveform, sample values
  Float_t m_ampl[BbcDefs::BBC_N_FEECH];              // raw amplitude

  std::vector<BbcSig> bbcsig;
  
  Float_t m_pmtq[BbcDefs::BBC_N_PMT]{};   // npe in each arm
  Float_t m_pmttt[BbcDefs::BBC_N_PMT]{};  // time in each arm
  Float_t m_pmttq[BbcDefs::BBC_N_PMT]{};  // time in each arm

  Short_t m_bbcn[2]{};      // num hits for each arm (north and south)
  Float_t m_bbcq[2]{};      // total charge (currently npe) in each arm
  Float_t m_bbct[2]{};      // time in arm
  Float_t m_bbcte[2]{};     // earliest hit time in arm
  Float_t m_bbcz {NAN};     // z-vertex
  Float_t m_bbczerr {NAN};  // z-vertex error
  Float_t m_bbct0 {NAN};    // start time
  Float_t m_bbct0err {NAN}; // start time error
  TH1 *hevt_bbct[2]{};  // time in each bbc, per event
  TF1 *gaussian {nullptr};
  Float_t _tres = NAN;  // time resolution of one channel

  TH2 *h2_tmax[2] = {};  // [0 == time ch, 1 == chg ch], max sample in evt vs ch

  float TRIG_SAMP[16];  // [board]

  TCanvas *ac;  // for plots used during debugging

};

#endif /* __BBCEVENT_H__ */
