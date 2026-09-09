#ifndef TPC_TPCCLUSTERMOVER_H
#define TPC_TPCCLUSTERMOVER_H

/*!
 * \file TpcClusterMover.h
 * \Moves TPC clusters to the readout TPC surface after distortion corrections
 * \author Tony Frawley, May 2022
 */
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/ActsTransformations.h>
#include <vector>

class PHG4TpcGeomContainer;
class PHCompositeNode;
class TpcClusterMover
{
 public:

  //! constructor
  TpcClusterMover() = default;

  void set_verbosity(int verb) { _verbosity = verb; }

  std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> processTrack(const std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> &global_in) const;

  //! Updates the assumed default geometry below to that contained in the
  //! cell geo
  void initialize_geometry(PHG4TpcGeomContainer*, ActsGeometry*, PHCompositeNode*);

 private:
  int get_circle_circle_intersection(double target_radius, double R, double X0, double Y0, double xclus, double yclus, double &x, double &y) const;

  bool get_moved_position(TrkrDefs::cluskey cluskey, TrkrCluster *cluster, std::vector<float> &fitpars, Acts::Vector3 &global, Acts::Vector3 &global_new, TrkrDefs::subsurfkey &new_subsurfkey) const;

  //! verbosity
  int _verbosity = 0;

  //! pointer to acts geometry container
  ActsGeometry *_tGeometry = nullptr;

  //! pointer to main node
  PHCompositeNode *_topNode = nullptr;
};

#endif
