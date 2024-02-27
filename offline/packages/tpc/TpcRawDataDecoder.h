// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCRAWDATADECODER_H
#define TPCRAWDATADECODER_H

//#include "TpcMap.h"

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>

#include <fun4all/SubsysReco.h>

#include <memory>
#include <string>

class CDBTTree;
class CDBInterface;
class PHCompositeNode;
class Fun4AllHistoManager;
// class TrkrHitSetContainer;
// class TrkrHitSet;
// class TrkrHit;
class TH2;
class TH3;
class TNtuple;

class TpcRawDataDecoder : public SubsysReco
{
 public:
  TpcRawDataDecoder(const std::string &name = "TpcRawDataDecoder");

  ~TpcRawDataDecoder() override;

  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
  */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
  */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  // int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  // int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  // int Reset(PHCompositeNode * /*topNode*/) override;

  // void Print(const std::string &what = "ALL") const override;
  // void setHistoFileName(const std::string &what = "./outputfile.root");

 protected:
  Fun4AllHistoManager *hm = nullptr;
  std::string _filename;

  static const int layercount = 16;
  static const int layeroffset = 7 + 16;
  int _ievent = 0;
  // TrkrHitSetContainer *m_hits = nullptr;
  // TrkrHitSet *m_hitset;
  // TrkrHit *m_hit = nullptr;

  // RawHitSetContainer *m_rawhits __attribute__ ((unused)) = nullptr;

  // TpcMap M;
  TNtuple *h_Alive = nullptr;
  CDBTTree *m_cdbttree = nullptr;
  CDBInterface *m_cdb = nullptr;

  struct ped_tpc_map
  {
    unsigned int CHN_ID;
    unsigned int FEE_ID;
    unsigned int MOD_ID;
    double PedMean;
    double PedStdi;
    unsigned int SEC_ID;
  };

  // std::map<unsigned int, struct ped_tpc_map> tmap;

  int starting_BCO;
  int rollover_value;
  int current_BCOBIN;

 private:
  int m_Debug = 0;
  // TH3*   _h_hit_XYT = nullptr;
  // TH3*   _h_hit_PT_ADCcut = nullptr;
  TH2 *_h_hit_XY = nullptr;
  TH2 *_h_hit_XY_ADCcut = nullptr;
};

#endif  // TPC_RAWDATADECODER_H
