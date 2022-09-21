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
    delete fz0;
    delete fz;
    delete fmm_55_2;
    delete fmm_56_2;
    delete fmm_3;
    delete fadcphi0;
    delete fadcphi1;
    delete fadcphi2;
    delete fadcphi2fine;
    delete fadcz0;
    delete fadcz1;
    delete fadcz2;
  };
  
  using error_t = std::pair<double, double>;

  error_t get_cluster_error(TrackSeed *seed, TrkrCluster* cluster, double cluster_r, TrkrDefs::cluskey key);
  error_t get_cluster_error(TrkrCluster* cluster,  TrkrDefs::cluskey key, double alpha, double beta);

  error_t get_simple_cluster_error(TrkrCluster* cluster, double cluster_r, TrkrDefs::cluskey key);
  error_t get_fix_tpc_cluster_error(TrkrCluster* cluster, TrkrDefs::cluskey key);
  error_t get_si_cluster_error(const TrkrCluster* cluster, TrkrDefs::cluskey key);

 private:

  TF1 *f0 = nullptr;
  TF1 *f1 = nullptr;
  TF1 *f2 = nullptr;
  TF1 *fz0 = nullptr;
  TF1 *fz = nullptr;
  TF1 *fmm_55_2 = nullptr;
  TF1 *fmm_56_2 = nullptr;
  TF1 *fmm_3 = nullptr;
  TF1 *fadcz0 = nullptr;
  TF1 *fadcz1 = nullptr;
  TF1 *fadcz2 = nullptr;
  TF1 *fadcphi0 = nullptr;
  TF1 *fadcphi1 = nullptr;
  TF1 *fadcphi2 = nullptr;
  TF1 *fadcphi2fine = nullptr;
  double pitcherr_phi_mvtx;
  double pitcherr_z_mvtx;

  double pitcherr_phi_intt;
  double pitcherr_z_intt;

  double pitcherr_phi_mm1;
  double pitcherr_z_mm1;

  double pitcherr_phi_mm2;
  double pitcherr_z_mm2;
  double scale_mvtx = 1.2;
  double scale_mvtx_z = 0.9;
  double scale_intt_3 = 0.97;
  double scale_intt_4 = 0.964;
  double scale_intt_5 = 0.894;
  double scale_intt_6 = 0.893;
  /*  double scale_tpc_0 = 1.1;
  double scale_tpc_1 = 1.2;
  double scale_tpc_2 = 1.3;
  double scale_tpc_0_z = 1.1;
  double scale_tpc_1_z = 1.07;
  double scale_tpc_2_z = 1.04;
  */
  double scale_mm_0 = 1.1;
  double scale_mm_1 = 1.5; 
};

#endif
