// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4INTT_PHG4INTTHITRECO_H
#define G4INTT_PHG4INTTHITRECO_H

#include <phparameter/PHParameterInterface.h>

#include <fun4all/SubsysReco.h>

#include <gsl/gsl_vector.h>  // for gsl_vector

#include <cmath>
#include <string>

class PHCompositeNode;
class PHG4InttTruthClusterizer;

// Process to cluster tracks:
class PHG4InttDigitizer; // digitize (no pruning that I can see)
class InttClusterizer;   // reconstruct

class PHG4InttHitReco : public SubsysReco, public PHParameterInterface
{
 public:
  PHG4InttHitReco(const std::string &name = "PHG4InttHitReco");

  ~PHG4InttHitReco() override;

  void SetTruthClusterizer (
        PHG4InttDigitizer* _digitiser
      , InttClusterizer*   _clusterizer
      , int                _verbosity =-1); // if left at default, it will use the value from Verbosity()

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

  double m_Tmin;
  double m_Tmax;
  double m_crossingPeriod;

  gsl_vector *m_LocalOutVec = nullptr;
  gsl_vector *m_PathVec = nullptr;
  gsl_vector *m_SegmentVec = nullptr;

  PHG4InttTruthClusterizer *m_truthclusterizer { nullptr };
};

#endif
