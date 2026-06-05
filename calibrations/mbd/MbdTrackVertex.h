#ifndef MBDTRACKVERTEX_H
#define MBDTRACKVERTEX_H

#include <fun4all/SubsysReco.h>

#include <TFile.h>
#include <TH1.h>
#include <THnSparse.h>
#include <TTree.h>

#include <cstdint>
#include <string>

class PHCompositeNode;

class MbdTrackVertex : public SubsysReco
{
 public:

  MbdTrackVertex(const std::string &name = "MbdTrackVertex");

  ~MbdTrackVertex() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void setOutputName(const std::string& name) { outFileName = name; };
  void SetTreeFlag(bool flag) { _treeflag = flag; }
  void SetTriggerMask(uint64_t mask) { _gl1_trigmask = mask; }

 private:

  TFile* outFile {nullptr};
  TTree* outTree {nullptr};
  TH1F* h_mbdtrkz {nullptr};
  TH1F* h_bz {nullptr};
  TH1F* h_trkz {nullptr};
  THnSparseF* h2_mbdtrkz {nullptr};
  std::string outFileName = "mbdtrk_vertex.root";

  Float_t _mbdVertex {std::numeric_limits<float>::quiet_NaN()};
  Float_t _trackerVertex {std::numeric_limits<float>::quiet_NaN()};
  UInt_t  _nTracks {std::numeric_limits<unsigned int>::quiet_NaN()};
  UInt_t  _nMBDVertex {std::numeric_limits<unsigned int>::quiet_NaN()};
  UInt_t  _nTRKVertex {std::numeric_limits<unsigned int>::quiet_NaN()};

  bool _hasMBD {false};
  bool _hasTRK {false};

  bool _treeflag {true};
  uint64_t _gl1_trigmask {0};
  int _counter{0};
  int _evt{0};
};

#endif // MBDTRACKVERTEX_H
