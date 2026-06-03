#ifndef MBD_MBDCALIBRECO_H
#define MBD_MBDCALIBRECO_H

#include "MbdDefs.h"

#include <fun4all/SubsysReco.h>

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <THnSparse.h>

class PHCompositeNode;
class MbdCalib;
class MbdPmtContainer;
class MbdOut;
class MbdGeom;
class Gl1Packet;
class EventHeader;
class RunHeader;
class TH1;
class TH2;
class TFile;
class TGraphErrors;

class MbdCalibReco : public SubsysReco
{
 public:
  MbdCalibReco(const std::string& name = "MbdCalibReco");
  ~MbdCalibReco() override = default;

  int Init(PHCompositeNode* topNode) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int EndRun(const int runnumber) override;

  void SetSubPass(const int s)          { _subpass = s; }
  void SetCalDir(const std::string& d)  { _caldir = d; }
  void SetCDBTag(const std::string& t)  { _cdbtag = t; }

 private:
  int  getNodes(PHCompositeNode* topNode);
  void BookHistograms();
  void DeleteHistograms();

  uint64_t    _mbias_trigger_mask{0};

  int         _subpass{0};
  int         _runnumber{0};
  std::string _caldir{"results"};
  std::string _rundir;   // _caldir/<runnumber>/
  std::string _cdbtag{}; // non-empty → download from CDB instead of local files

  std::unique_ptr<MbdCalib> _mbdcal;
  MbdPmtContainer* _mbdpmts{nullptr};
  MbdOut*          _mbdout{nullptr};
  MbdGeom*         _mbdgeom{nullptr};
  EventHeader*     _evtheader{nullptr};
  RunHeader*       _runheader{nullptr};
  Gl1Packet*       _gl1packet{nullptr};

  std::array<TH1*, MbdDefs::MBD_N_PMT> h_tt{};
  std::array<TH1*, MbdDefs::MBD_N_PMT> h_tq{};
  std::array<TH1*, MbdDefs::MBD_N_PMT> h_qp{};
  std::array<THnSparseF*, MbdDefs::MBD_N_PMT> h2_slew{};
  TH2* h2_tt{nullptr};
  TH2* h2_tq{nullptr};

  std::unique_ptr<TFile> _outfile{nullptr};
};

#endif  // MBD_MBDCALIBRECO_H
