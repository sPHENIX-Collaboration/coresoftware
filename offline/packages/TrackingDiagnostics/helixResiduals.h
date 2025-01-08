#ifndef HELIXRESIDUALS_H
#define HELIXRESIDUALS_H

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackermillepedealignment/HelicalFitter.h>

#include <fun4all/SubsysReco.h>

class TFile;
class TNtuple;

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;

class helixResiduals : public SubsysReco
{
 public:
  helixResiduals(const std::string& name = "helixResiduals");
  ~helixResiduals() {}

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* /*topNode*/) override;
  int End(PHCompositeNode* /*topNode*/) override;
  void set_output_file(const std::string& outputfile) { filepath = outputfile; }

 private:
  int getNodes(PHCompositeNode* topNode);

  void fill_residuals(TrackSeed* seed, int seed_id, bool isTpc);

  HelicalFitter* _fitter = nullptr;
  ActsGeometry* tGeometry = nullptr;

  TNtuple* ntp_residuals = nullptr;

  SvtxTrackMap* _tracks = nullptr;
  TrkrClusterContainer* _clusters = nullptr;
  TrackSeedContainer* _tpc_seeds = nullptr;
  TrackSeedContainer* _si_seeds = nullptr;

  TFile* fout = nullptr;

  std::string filepath = "";
};

#endif  // HELIXRESIDUALS_H
