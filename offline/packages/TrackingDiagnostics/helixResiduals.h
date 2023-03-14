#ifndef HELIXRESIDUALS_H
#define HELIXRESIDUALS_H

#include <fun4all/SubsysReco.h>

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase_historic/TrackSeed_v1.h>
#include <trackermillepedealignment/HelicalFitter.h>

class TFile;
class TNtuple;

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;

class helixResiduals : public SubsysReco
{
public:
  helixResiduals(const std::string &name = "helixResiduals", const std::string &outputFile = "residuals.root");
  ~helixResiduals(){}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode * /*topNode*/) override;
  int End(PHCompositeNode * /*topNode*/) override;

private:  

  int getNodes(PHCompositeNode *topNode);

  void fill_residuals(TrackSeed* seed, int seed_id, bool isTpc);

  HelicalFitter* _fitter;
  ActsGeometry* tGeometry;

  TNtuple* ntp_residuals;

  SvtxTrackMap* _tracks = nullptr;
  TrkrClusterContainer* _clusters = nullptr;
  TrackSeedContainer* _tpc_seeds = nullptr;
  TrackSeedContainer* _si_seeds = nullptr;

  TFile *fout;

  std::string _outputfile;
};

#endif  // HELIXRESIDUALS_H 
