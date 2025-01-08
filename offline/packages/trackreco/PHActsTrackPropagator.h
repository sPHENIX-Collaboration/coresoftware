// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHACTSTRACKPROPAGATOR_H
#define PHACTSTRACKPROPAGATOR_H

#include "ActsPropagator.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>

#include <trackbase/ActsGeometry.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Utilities/Result.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>

#include <string>

class PHCompositeNode;
class SvtxTrackMap;
class ActsGeometry;
class SvtxVertexMap;

class PHActsTrackPropagator : public SubsysReco
{
 public:
  using BoundTrackParam =
      const Acts::BoundTrackParameters;
  using BoundTrackParamResult = ActsPropagator::BTPPairResult;
  using SurfacePtr = std::shared_ptr<const Acts::Surface>;
  using Trajectory = ActsExamples::Trajectories;

  PHActsTrackPropagator(const std::string &name = "PHActsTrackPropagator");

  ~PHActsTrackPropagator() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;
  void setPropagationLayer(unsigned int layer) { m_sphenixLayer = layer; }

 private:
  int getNodes(PHCompositeNode *topNode);
  BoundTrackParamResult propagateTrack(
      const Acts::BoundTrackParameters &params);
  void convertsPHENIXLayerToActsLayer(unsigned int &actsvolume,
                                      unsigned int &actslayer);
  void addTrackState(BoundTrackParamResult &params,
                     SvtxTrack *svtxTrack);

  /// Objects containing the Acts track fit results
  ActsGeometry *m_tGeometry = nullptr;
  SvtxTrackMap *m_trackMap = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;

  unsigned int m_sphenixLayer = std::numeric_limits<unsigned int>::max();
};

#endif  // PHACTSTRACKPROPAGATOR_H
