#ifndef TRACKBASE_CLUSTERERRORPARA_H
#define TRACKBASE_CLUSTERERRORPARA_H

#include "TrkrDefs.h"

#include <tuple>
#include <utility>
#include <vector>
#include <string>
#include <TF1.h>
//class TF1;
class TrackSeed;
class TrkrCluster;

class ClusterErrorPara
{

  public:
  ClusterErrorPara();
 // delete copy ctor and assignment operator (cppcheck)
  explicit ClusterErrorPara(const ClusterErrorPara&) = delete;
  ClusterErrorPara& operator=(const ClusterErrorPara&) = delete;

  ~ClusterErrorPara(){
    delete f0;
    delete f1;
    delete f2;
    delete fz;
  };
  
  using error_t = std::pair<double, double>;

  error_t get_cluster_error(TrackSeed *seed, TrkrCluster* cluster, double cluster_r, TrkrDefs::cluskey key);
  error_t get_simple_cluster_error(TrkrCluster* cluster, double cluster_r, TrkrDefs::cluskey key);
  error_t get_fix_tpc_cluster_error(TrkrCluster* cluster, TrkrDefs::cluskey key);
  error_t get_si_cluster_error(const TrkrCluster* cluster, TrkrDefs::cluskey key);

 private:

  TF1 *f0 = nullptr;
  TF1 *f1 = nullptr;
  TF1 *f2 = nullptr;
  TF1 *fz = nullptr;
  double pitcherr_phi_mvtx;
  double pitcherr_z_mvtx;

  double pitcherr_phi_intt;
  double pitcherr_z_intt;

  double pitcherr_phi_mm1;
  double pitcherr_z_mm1;

  double pitcherr_phi_mm2;
  double pitcherr_z_mm2;

};

#endif
