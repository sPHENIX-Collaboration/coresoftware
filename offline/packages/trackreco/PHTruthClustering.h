/*!
 *  \file		PHTruthClustering.h
 *  \brief		Clustering using truth info
 *  \author		Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TRACKRECO_PHTRUTHCLUSTERING_H
#define TRACKRECO_PHTRUTHCLUSTERING_H

#include <fun4all/SubsysReco.h>

// rootcint barfs with this header so we need to hide it
#include <gsl/gsl_rng.h>

class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHG4CylinderGeomContainer;
class PHG4CylinderCellGeomContainer;
//class PHG4VtxPoint;
class TrkrCluster;


#include <string>             // for string
#include <vector>
#include <map>
#include <set>
#include <memory>

// forward declarations
class PHCompositeNode;
class PHG4TruthInfoContainer;

class PHTruthClustering  : public SubsysReco
{
public:
  PHTruthClustering(const std::string &name = "PHTruthClustering");
  virtual ~PHTruthClustering();

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);


private:
/// fetch node pointers
int GetNodes(PHCompositeNode *topNode);

std::map<unsigned int, TrkrCluster* > all_truth_clusters(PHG4Particle* particle);
std::set<PHG4Hit*> all_truth_hits(PHG4Particle* particle);

  void LayerClusterG4Hits(std::set<PHG4Hit*> truth_hits, std::vector<PHG4Hit*> &contributing_hits, std::vector<double> &contributing_hits_energy, std::vector<std::vector<double>> &contributing_hits_entry, std::vector<std::vector<double>> &contributing_hits_exit, float layer, float &x, float &y, float &z,  float &t, float &e);
  
  void G4ClusterSize(unsigned int layer, std::vector<std::vector<double>> contributing_hits_entry,std::vector<std::vector<double>> contributing_hits_exit, float &g4phisize, float &g4zsize);

  float line_circle_intersection(float x[], float y[], float z[], float radius);

  int iclus = 0;

PHG4TruthInfoContainer *_g4truth_container{nullptr};

  PHG4HitContainer* _g4hits_svtx{nullptr};
  PHG4HitContainer* _g4hits_mms{nullptr};
  PHG4HitContainer* _g4hits_tracker{nullptr};
  PHG4HitContainer* _g4hits_maps{nullptr};

  PHG4CylinderCellGeomContainer* _tpc_geom_container{nullptr};
  PHG4CylinderGeomContainer *_intt_geom_container{nullptr};
  PHG4CylinderGeomContainer* _mvtx_geom_container{nullptr};
  PHG4CylinderGeomContainer* _mms_geom_container{nullptr};

 const unsigned int _nlayers_maps = 3;
  const unsigned int _nlayers_intt = 4;
  const unsigned int _nlayers_tpc = 48;
  const unsigned int _nlayers_mms = 2;

  double clus_err_rphi[57] = {0};
  double clus_err_z[57] = {0};

  double mvtx_clus_err_rphi = 5e-04;
  double mvtx_clus_err_z = 5e-04;
  double intt_clus_err_rphi = 25e-04;
  double intt_clus_err_z = 1.0;
  double tpc_inner_clus_err_rphi = 200e-04;
  double tpc_inner_clus_err_z = 750e-04;
  double tpc_outer_clus_err_rphi = 150e-04;
  double tpc_outer_clus_err_z = 500e-04;
  double mms_layer55_clus_err_rphi = 100e-04;
  double mms_layer55_clus_err_z = 25.0;
  double mms_layer56_clus_err_rphi = 12.5;
  double mms_layer56_clus_err_z = 200e-04;
};

#endif
