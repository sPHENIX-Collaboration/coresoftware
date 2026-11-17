#ifndef MBDTRACKVERTEX_H
#define MBDTRACKVERTEX_H

#include <fun4all/SubsysReco.h>

#include <THnSparse.h>  // THnSparseF is a typedef, can't forward declare

#include <Rtypes.h>

#include <cstdint>
#include <limits>
#include <string>

class PHCompositeNode;
class TFile;
class TTree;
class TH1F;
class TH2F;

/**
 * Compare the MBD z-vertex with the tracker (SvtxVertexMap) z-vertex.
 *
 * Inputs (per event):
 *   MbdOut, MbdVertexMap  : from MbdReco run on DST_CALOFITTING
 *   SvtxVertexMap         : from DST_TRKR_TRACKS (polyseed + Acts fit + vertexing)
 *   Gl1Packet             : node 14001 (DST_CALOFITTING) or GL1RAWHIT (DST_TRKR_TRACKS)
 *
 * The tracker vertex used is the one with the most tracks in the selected
 * beam crossing (default 0, the triggered crossing).
 */
class MbdTrackVertex : public SubsysReco
{
 public:
  MbdTrackVertex(const std::string &name = "MbdTrackVertex");

  ~MbdTrackVertex() override = default;

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  void setOutputName(const std::string &name) { outFileName = name; }
  void SetTreeFlag(bool flag) { _treeflag = flag; }
  void SetTriggerMask(uint64_t mask) { _gl1_trigmask = mask; }
  void SetBeamCrossing(short crossing) { _trk_crossing = crossing; }
  void SetMinTracks(unsigned int ntrk) { _min_ntracks = ntrk; }

 private:
  TFile *outFile{nullptr};
  TTree *outTree{nullptr};
  TH1F *h_mbdtrkz{nullptr};
  TH1F *h_bz{nullptr};
  TH1F *h_trkz{nullptr};
  TH1F *h_bt0{nullptr};
  THnSparseF *h2_mbdtrkz{nullptr};
  TH2F *h2_dzvstrkz{nullptr};
  std::string outFileName = "mbdtrk_vertex.root";

  bool _hasMBD{false};
  bool _hasTRK{false};

  // Tree variables
  Int_t _run{0};
  Int_t _evt{0};
  ULong64_t _strig{0};  // GL1 scaled trigger vector
  Float_t _mbdVertex{std::numeric_limits<float>::quiet_NaN()};
  Float_t _mbdT0{std::numeric_limits<float>::quiet_NaN()};
  Float_t _trackerVertex{std::numeric_limits<float>::quiet_NaN()};
  Float_t _trackerX{std::numeric_limits<float>::quiet_NaN()};
  Float_t _trackerY{std::numeric_limits<float>::quiet_NaN()};
  Float_t _trackerChi2ndf{std::numeric_limits<float>::quiet_NaN()};
  UInt_t _nMBDVertex{0};
  UInt_t _nTRKVertex{0};     // number of tracker vertices in selected crossing
  UInt_t _nTRKVertexAll{0};  // number of tracker vertices in all crossings
  UInt_t _nTracks{0};
  Short_t _nMBD[2]{0};
  Float_t _qMBD[2]{0};

  bool _treeflag{true};
  uint64_t _gl1_trigmask{0};
  short _trk_crossing{0};
  unsigned int _min_ntracks{1};
  int _counter{0};
};

#endif  // MBDTRACKVERTEX_H
