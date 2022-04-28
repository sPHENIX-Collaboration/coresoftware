/*!
 *  \file		PHTruthTrackSeeding.h
 *  \brief		Vertexing using truth info
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHTRUTHTRACKSEEDING_H
#define TRACKRECO_PHTRUTHTRACKSEEDING_H

#include "PHTrackSeeding.h"
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/TrkrDefs.h>
#include <string>  // for string
#include <vector>

// forward declarations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4HitContainer;
class TrkrHitTruthAssoc;
class TrkrClusterContainer;
class TrkrClusterCrossingAssoc;
class SvtxClusterEval;
class SvtxTrack;

//class SvtxHitMap;
//class PHG4CellContainer;

/// \class PHTruthTrackSeeding
///
/// \brief Vertexing using truth info
///

class PHTruthTrackSeeding : public PHTrackSeeding
{
 public:
  PHTruthTrackSeeding(const std::string& name = "PHTruthTrackSeeding");

  unsigned int get_min_clusters_per_track() const
  {
    return _min_clusters_per_track;
  }

  void set_min_clusters_per_track(unsigned int minClustersPerTrack)
  {
    _min_clusters_per_track = minClustersPerTrack;
  }

  void set_min_layer(unsigned int minLayer)
  {
    _min_layer = minLayer;
  }

  void set_max_layer(unsigned int maxLayer)
  {
    _max_layer = maxLayer;
  }

  //! minimal truth momentum cut
  double get_min_momentum() const
  {
    return _min_momentum;
  }

  //! minimal truth momentum cut
  void set_min_momentum(double m)
  {
    _min_momentum = m;
  }
void helicalTrackFit(const bool helicalTrackFit)
{m_helicalTrackFit = helicalTrackFit; }
 protected:
  int Setup(PHCompositeNode* topNode) override;

  int Process(PHCompositeNode* topNode) override;

  int End() override;

 private:
  /// fetch node pointers
  int GetNodes(PHCompositeNode* topNode);

  /// get crossing id from intt clusters associated to track
  /* this is a copy of the code in PHTruthSiliconAssociation */
  std::set<short int> getInttCrossings(SvtxTrack*) const;

  void circleFitSeed(std::vector<TrkrDefs::cluskey> clusters,
		     double& x, double& y, double&z,
		       double& px, double& py, double& pz, int charge) const;
  std::vector<Acts::Vector3> circleFitByTaubin(std::vector<TrkrDefs::cluskey>& clusters,
						double& R, double& X0, double& Y0) const ;
  void findRoot(const double& R, const double& X0, const double& Y0,
		double& x, double& y) const ;
  void lineFit(std::vector<Acts::Vector3>& clusterPositions,
	       double& A, double& B) const ;
      
  PHG4TruthInfoContainer* m_g4truth_container = nullptr;
  TrkrClusterContainer *m_clusterMap = nullptr;
  TrkrClusterCrossingAssoc *m_cluster_crossing_map = nullptr;
  PHG4HitContainer* phg4hits_tpc = nullptr;
  PHG4HitContainer* phg4hits_intt = nullptr;
  PHG4HitContainer* phg4hits_mvtx = nullptr;
  PHG4HitContainer* phg4hits_micromegas = nullptr;

  TrkrHitTruthAssoc* hittruthassoc = nullptr;
  SvtxClusterEval* _clustereval;

  unsigned int _min_clusters_per_track = 3;
  unsigned int _min_layer = 0;
  unsigned int _max_layer = 60;

  /// Option to perform a helical track fit to get track parameters for
  /// truth seeded track
  bool m_helicalTrackFit = false;

  //! minimal truth momentum cut (GeV)
  double _min_momentum = 50e-3;

  ActsTrackingGeometry *tgeometry = nullptr;
  ActsSurfaceMaps *surfmaps = nullptr;

};

#endif
