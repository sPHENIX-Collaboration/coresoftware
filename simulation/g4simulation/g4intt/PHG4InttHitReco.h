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

class TrkrTruthTrackContainer;
class TrkrClusterContainer;
class PHG4InttTruthClusterizer;


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

  PHG4InttTruthClusterizer* get_truth_clusterizer() { return m_truthclusterizer; };
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

  PHG4InttTruthClusterizer* m_truthclusterizer;

  gsl_vector *m_LocalOutVec = nullptr;
  gsl_vector *m_PathVec = nullptr;
  gsl_vector *m_SegmentVec = nullptr;

  // needed for clustering truth tracks
  private:
  TrkrTruthTrackContainer* m_truthtracks     { nullptr }; // output truth tracks
  TrkrClusterContainer*    m_truthclusters   { nullptr }; // output clusters indexed to TrkrDefs::cluskeys in m_truthtracks
  PHG4TruthInfoContainer*  m_truthinfo       { nullptr };
  int                      m_trkid           { -1      };
  bool                     m_is_emb          { false   };
  TrkrTruthTrack*          m_current_track   { nullptr };
  const int                m_cluster_version { 4 };
  TrkrHitSetContainer*     m_truth_hits; // generate and delete a container for each truth track
  std::map<TrkrDefs::hitsetkey,unsigned int> m_hitsetkey_cnt {}; // counter for making ckeys form hitsetkeys

  void truthcheck_g4hit       ( PHG4Hit*, PHCompositeNode* topNode );
  void addtruthhitset         ( TrkrDefs::hitsetkey, TrkrDefs::hitkey, float neffelectrons );
  void clusterize_truthtrack  ( PHCompositeNode* topNode );
  void end_event_truthcluster ( PHCompositeNode* topNode );

  double m_truth_pixelthreshold { 0.01 };
  public:
  void set_truth_pixelthreshold (double val) { m_truth_pixelthreshold = val; };

};

#endif
