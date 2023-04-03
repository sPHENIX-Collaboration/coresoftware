// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4INTT_PHG4INTTHITRECO_H
#define G4INTT_PHG4INTTHITRECO_H

#include "TruthInttClusterBuilder.h"

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_vector.h>  // for gsl_vector

#include <cmath>
#include <string>

class PHCompositeNode;

class TrkrTruthTrackContainer;
class TrkrClusterContainer;


class PHG4InttHitReco : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4InttHitReco(const std::string &name = "PHG4InttHitReco");

  ~PHG4InttHitReco() override;
  //! module initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  //! set default parameter values
  void SetDefaultParameters() override;

  void Detector(const std::string &d) { m_Detector = d; }

 protected:
  std::string m_Detector = "INTT";
  std::string m_HitNodeName;
  std::string m_CellNodeName;
  std::string m_GeoNodeName;

  TrkrTruthTrackContainer* m_truthtracks { nullptr };
  TrkrClusterContainer*    m_truthclusters { nullptr };

  double m_Tmin;
  double m_Tmax;
  double m_crossingPeriod;

  TruthInttClusterBuilder* m_truth_clusterer { nullptr };

  gsl_vector *m_LocalOutVec = nullptr;
  gsl_vector *m_PathVec = nullptr;
  gsl_vector *m_SegmentVec = nullptr;
};

#endif
