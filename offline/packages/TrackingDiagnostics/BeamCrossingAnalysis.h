#ifndef BEAMCROSSINGANALYSIS_H
#define BEAMCROSSINGANALYSIS_H

#include <fun4all/SubsysReco.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

#include <globalvertex/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrackMap.h>

class TFile;
class TH1D;
class TNtuple;

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class TrackVertexCrossingAssoc;
class TNTuple;
class TH1D;

class BeamCrossingAnalysis : public SubsysReco
{
 public:
  BeamCrossingAnalysis(const std::string& name = "BeamCrossingAnalysis");
  virtual ~BeamCrossingAnalysis() {}

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* /*topNode*/) override;

  void set_output_file(const std::string& outputfile) { filepath = outputfile; }

 private:

  int getNodes(PHCompositeNode* topNode);

  SvtxTrackMap* m_svtxTrackMap = nullptr;
  SvtxVertexMap* m_vertexMap = nullptr;
  TrackVertexCrossingAssoc* m_track_vertex_crossing_map{nullptr};

  TNtuple *ntp_vertex{nullptr};
  TNtuple *ntp_track{nullptr};
  TH1D* hcross{nullptr};

  std::string filepath = "";
  bool _require_mvtx = true;
  double _qual_cut = 10.0;
  TFile* fout = nullptr;

};

#endif  // BEAMCROSSINGANALYSIS_H
