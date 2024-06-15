#ifndef TRACKBASE_CLUSTERERRORPARA_H
#define TRACKBASE_CLUSTERERRORPARA_H

#include "TrkrDefs.h"

#include <TF1.h>
#include <trackbase/TrkrClusterv5.h>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
// class TF1;
class TrkrCluster;

class ClusterErrorPara
{
 public:
  ClusterErrorPara();
  // delete copy ctor and assignment operator (cppcheck)
  explicit ClusterErrorPara(const ClusterErrorPara &) = delete;
  ClusterErrorPara &operator=(const ClusterErrorPara &) = delete;

  ~ClusterErrorPara()
  {
    delete f0;
    delete f1;
    delete f2;
    delete f0fine;
    delete f1fine;
    delete f2fine;
    delete f2fine2;
    delete fz0;
    delete fz1;
    delete fz2;
    delete fz0fine;
    delete fz1fine;
    delete fz2fine;
    delete fmm_55_2;
    delete fmm_56_2;
    delete fmm_3;
    delete fadcphi0;
    delete fadcphi0fine;
    delete fadcphi1;
    delete fadcphi1fine;
    delete fadcphi2;
    delete fadcphi2fine1;
    delete fadcphi2fine2;
    delete fadcz0;
    delete fadcz1;
    delete fadcz2;
    delete fadcz0fine;
    delete fadcz1fine;
    delete fadcz2fine;
  };

  using error_t = std::pair<double, double>;

  error_t get_clusterv5_modified_error(TrkrCluster *cluster, double cluster_r, TrkrDefs::cluskey key);
  error_t get_cluster_error(TrkrCluster *cluster, double cluster_r, TrkrDefs::cluskey key, float qOverR, float slope);
  error_t get_cluster_error(TrkrCluster *cluster, TrkrDefs::cluskey key, double alpha, double beta);

  error_t get_simple_cluster_error(TrkrCluster *cluster, double cluster_r, TrkrDefs::cluskey key);
  error_t get_fix_tpc_cluster_error(TrkrCluster *cluster, TrkrDefs::cluskey key);
  error_t get_si_cluster_error(const TrkrCluster *cluster, TrkrDefs::cluskey key);
  double mm_phi_error(int layer, double alpha, TrkrCluster *cluster);
  double mm_z_error(int layer, double beta, TrkrCluster *cluster);
  double mvtx_phi_error(TrkrCluster *cluster);
  double mvtx_phi_error(const TrkrCluster *cluster);
  double mvtx_z_error(TrkrCluster *cluster);
  double mvtx_z_error(const TrkrCluster *cluster);
  double intt_phi_error(int layer, TrkrCluster *cluster);
  double intt_z_error(TrkrCluster *cluster);
  double intt_phi_error(int layer, const TrkrCluster *cluster);
  double intt_z_error(const TrkrCluster *cluster);
  double tpc_phi_error(int layer, double alpha, TrkrCluster *cluster);
  double tpc_z_error(int layer, double beta, TrkrCluster *cluster);

 private:
  TF1 *f0 = nullptr;
  TF1 *f1 = nullptr;
  TF1 *f2 = nullptr;
  TF1 *f0fine = nullptr;
  TF1 *f1fine = nullptr;
  TF1 *f2fine = nullptr;
  TF1 *f2fine2 = nullptr;
  TF1 *fz0 = nullptr;
  TF1 *fz1 = nullptr;
  TF1 *fz2 = nullptr;
  TF1 *fz0fine = nullptr;
  TF1 *fz1fine = nullptr;
  TF1 *fz2fine = nullptr;
  TF1 *fmm_55_2 = nullptr;
  TF1 *fmm_56_2 = nullptr;
  TF1 *fmm_3 = nullptr;
  TF1 *fadcz0 = nullptr;
  TF1 *fadcz1 = nullptr;
  TF1 *fadcz2 = nullptr;
  TF1 *fadcz0fine = nullptr;
  TF1 *fadcz1fine = nullptr;
  TF1 *fadcz2fine = nullptr;
  TF1 *fadcphi0 = nullptr;
  TF1 *fadcphi0fine = nullptr;
  TF1 *fadcphi1 = nullptr;
  TF1 *fadcphi1fine = nullptr;
  TF1 *fadcphi2 = nullptr;
  TF1 *fadcphi2fine1 = nullptr;
  TF1 *fadcphi2fine2 = nullptr;
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
  double pull_fine_phi[60]{};
  double pull_fine_z[60]{};
};

#endif
