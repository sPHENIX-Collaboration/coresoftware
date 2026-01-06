#ifndef G4EVAL_TRKRCLUSTERISMATCHER_H
#define G4EVAL_TRKRCLUSTERISMATCHER_H
// Principle use:
//   Does comparison of two clusters
//   For now, does comparison by relative separation in T and RxPhi space
//   Option to compare by a relative number of bins|pixels, or by relative
//   overall size of the clusters (using the largest clustersize).
//

#include <trackbase/TrkrDefs.h>

#include <array>
#include <string>

class ActsGeometry;
class TrkrCluster;
class PHCompositeNode;
class TrkrClusterContainer;

class TrkrClusterIsMatcher
{
 public:
  TrkrClusterIsMatcher() = default;
  ~TrkrClusterIsMatcher() = default;

  bool operator()(TrkrDefs::cluskey key_T, TrkrDefs::cluskey key_R);

  int init(PHCompositeNode* topNode,
           const std::string& name_phg4_clusters = "TRKR_TRUTHCLUSTERCONTAINER",
           const std::string& name_reco_clusters = "TRKR_CLUSTER");

  void set_tol_phi_MVTX(float _val);
  void set_tol_phi_INTT(float _val);
  void set_tol_phi_TPC(float _val);

  void set_tol_z_MVTX(float _val);
  void set_tol_t_TPC(float _val);

  // options for if tolerance is multipled by width of one pixel|bin, or by total number of pixel|bin's
  // in the larger of two clusters compared
  bool single_pixel_phi_MVTX{false};  // default to pitch*max(N_pixels_M,N_pixels_T)*tol_MVTX
  bool single_pixel_phi_INTT{false};  // ... same as for MVTX
  bool single_bin_phi_TPC{true};      // default to pitch*tol_TPC

  bool single_pixel_z_MVTX{false};  // default to pitch*max(N_pixels_M,N_pixels_T)*tol_MVTX
  bool single_pixel_z_INTT{false};  // ... same as for MVTX
  bool single_bin_t_TPC{true};      // default to pitch*tol_TPC

  // For efficiency and convenience, hang on to minimal information about most recent comparison:
  // --------------------------------------------------------------------------------------------
  TrkrCluster* clus_T{nullptr};  // _T for truth cluster
  TrkrCluster* clus_R{nullptr};  // _R for reco  cluster
  int layer{0};
  int det_0123{0};
  bool is_match{false};  // also the return value of operator()
  float dphi{0.};

  // multipliers for pixels and bins to cm and times
  std::array<float, 56> pitch_phi{0.};  // the phistep squared
  float pitch_z_MVTX{0.};               // pixel width for MVTX
  float step_t_TPC{0.};                 // time bin width for PTC

  std::array<float, 56> tol_pitch_phi{0.};  // pitched times tol_phi_{MVTX,INTT,TPC}
  float tol_pitch_z_MVTX{0.};               // tol * pixel width for MVTX
  float tol_step_t_TPC{0.};                 // tol * time bin width for PTC

  // containers to get the clusters from
  TrkrClusterContainer* m_TruthClusters{nullptr};
  TrkrClusterContainer* m_RecoClusters{nullptr};
  ActsGeometry* m_ActsGeometry{nullptr};

  /* public: */
  /* TrkrClusterContainer* get_TruthClusters() { return m_TruthClusters; }; */
  /* TrkrClusterContainer* get_RecoClusters () { return m_RecoClusters; }; */
  /* ActsGeometry* get_ActsGeometry () { return m_ActsGeometry; }; */
  /* private: */

  // tolerances for comparison in MVTX, INTT, TPC, and TPOT
  float tol_phi_MVTX{0.5};
  float tol_phi_INTT{0.5};
  float tol_phi_TPC{1.0};

  float tol_z_MVTX{0.5};
  float tol_t_TPC{1.0};
};

#endif
