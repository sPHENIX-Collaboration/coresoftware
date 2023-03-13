#ifndef ACTSEVALUATOR_H
#define ACTSEVALUATOR_H

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>

#include <Acts/Utilities/Helpers.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>
#include <ActsExamples/EventData/TrkrClusterSourceLink.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <ActsExamples/TrackFitting/TrkrClusterFittingAlgorithm.hpp>
#pragma GCC diagnostic pop

class TTree;
class TFile;
class PHG4Particle;
class SvtxTrack;
class SvtxVertexMap;
class SvtxEvalStack;
class SvtxTrackMap;
class PHG4TruthInfoContainer;
class TrkrClusterContainer;
class SvtxEvaluator;

#include <map>
#include <string>
#include <vector>

using SourceLink = ActsExamples::TrkrClusterSourceLink;
using FitResult = Acts::KalmanFitterResult;
using Trajectory = ActsExamples::Trajectories;
using Measurement = Acts::Measurement<Acts::BoundIndices,2>;
using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::theta;


/**
 * This class is an analyzing class for the Acts track fitting, and produces
 * a tree with many branches useful for debugging what Acts is doing. 
 * The truth G4Particle, reconstructed Acts track fit result from 
 * PHActsTrkFitter, and the comparison between truth and reconstructed states
 * throughout the track fit are written out to the tree.
 * Note that this module works and outputs in Acts units of mm and GeV, 
 */
class ActsEvaluator : public SubsysReco
{
 public:
  ActsEvaluator(const std::string &name = "ActsEvaluator.root",
                     SvtxEvaluator *svtxEvaluator = nullptr);
  ~ActsEvaluator() override;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void setEvalCKF(bool evalCKF) {m_evalCKF = evalCKF;}

 private:
  int getNodes(PHCompositeNode *topNode);
  
  /// Function to evaluate Trajectories fit results from the KF
  void evaluateTrackFits(PHCompositeNode *topNode);
  
  void initializeTree();

  void fillG4Particle(PHG4Particle *part);

  void fillProtoTrack(SvtxTrack* track, PHCompositeNode *topNode);

  void fillFittedTrackParams(const Trajectory traj,
			     const size_t &trackTip,
			     const Acts::Vector3 vertex);

  void visitTrackStates(const Acts::MultiTrajectory& traj,
			const size_t &trackTip, 
			PHCompositeNode *topNode);

  void clearTrackVariables();
  
  void calculateDCA(const Acts::BoundTrackParameters param, 
		    const Acts::Vector3 vertex);

  Surface getSurface(TrkrDefs::cluskey cluskey,TrkrDefs::subsurfkey surfkey);
  Surface getSiliconSurface(TrkrDefs::hitsetkey hitsetkey);
  Surface getTpcSurface(TrkrDefs::hitsetkey hitsetkey, TrkrDefs::subsurfkey surfkey);
  Surface getMMSurface(TrkrDefs::hitsetkey hitsetkey);


  Acts::Vector3 getGlobalTruthHit(PHCompositeNode *topNode, 
				   TrkrDefs::cluskey cluskey,
				   float &_gt);

  ActsSurfaceMaps *m_surfMaps{nullptr};

  SvtxEvaluator *m_svtxEvaluator{nullptr};
  PHG4TruthInfoContainer *m_truthInfo{nullptr};
  SvtxTrackMap *m_trackMap{nullptr}, *m_actsProtoTrackMap{nullptr};
  SvtxEvalStack *m_svtxEvalStack{nullptr};
  std::map<const unsigned int, std::map<const size_t, 
    const unsigned int>> *m_actsTrackKeyMap{nullptr};
  std::map<const unsigned int, Trajectory> *m_actsFitResults{nullptr};

  ActsTrackingGeometry *m_tGeometry{nullptr};
  SvtxVertexMap *m_vertexMap{nullptr};
  TrkrClusterContainer *m_clusterContainer{nullptr};

  /// boolean indicating whether or not to evaluate the CKF or
  /// the KF. Must correspond with what was run to do fitting
  /// i.e. PHActsTrkFitter or PHActsTrkProp
  bool m_evalCKF;

  TFile *m_trackFile{nullptr};
  TTree *m_trackTree{nullptr};

  /// Acts tree values
  int m_eventNr{0};
  int m_trajNr{0};
  int m_trackNr{0};

  unsigned long m_t_barcode{0};  /// Truth particle barcode
  int m_t_charge{0};             /// Truth particle charge
  float m_t_time{0};             /// Truth particle time
  float m_t_vx{-99.};            /// Truth particle vertex x
  float m_t_vy{-99.};            /// Truth particle vertex y
  float m_t_vz{-99.};            /// Truth particle vertex z
  float m_t_px{-99.};            /// Truth particle initial momentum px
  float m_t_py{-99.};            /// Truth particle initial momentum py
  float m_t_pz{-99.};            /// Truth particle initial momentum pz
  float m_t_theta{-99.};         /// Truth particle initial momentum theta
  float m_t_phi{-99.};           /// Truth particle initial momentum phi
  float m_t_pT{-99.};            /// Truth particle initial momentum pT
  float m_t_eta{-99.};           /// Truth particle initial momentum eta

  std::vector<float> m_t_x;  /// Global truth hit position x
  std::vector<float> m_t_y;  /// Global truth hit position y
  std::vector<float> m_t_z;  /// Global truth hit position z
  std::vector<float> m_t_r;  /// Global truth hit position r
  std::vector<float>
      m_t_dx;  /// Truth particle direction x at global hit position
  std::vector<float>
      m_t_dy;  /// Truth particle direction y at global hit position
  std::vector<float>
      m_t_dz;  /// Truth particle direction z at global hit position

  std::vector<float> m_t_eLOC0;   /// truth parameter eLOC_0
  std::vector<float> m_t_eLOC1;   /// truth parameter eLOC_1
  std::vector<float> m_t_ePHI;    /// truth parameter ePHI
  std::vector<float> m_t_eTHETA;  /// truth parameter eTHETA
  std::vector<float> m_t_eQOP;    /// truth parameter eQOP
  std::vector<float> m_t_eT;      /// truth parameter eT

  int m_nSharedHits{0};             /// number of shared hits in the track fit
  int m_nHoles{0};                  /// number of holes in the track fit
  int m_nOutliers{0};               /// number of outliers in the track fit
  int m_nStates{0};                 /// number of all states
  int m_nMeasurements{0};           /// number of states with measurements
  std::vector<int> m_volumeID;      /// volume identifier
  std::vector<int> m_layerID;       /// layer identifier
  std::vector<int> m_moduleID;      /// surface identifier
  std::vector<float> m_lx_hit;      /// uncalibrated measurement local x
  std::vector<float> m_ly_hit;      /// uncalibrated measurement local y
  std::vector<float> m_x_hit;       /// uncalibrated measurement global x
  std::vector<float> m_y_hit;       /// uncalibrated measurement global y
  std::vector<float> m_z_hit;       /// uncalibrated measurement global z
  std::vector<float> m_res_x_hit;   /// hit residual x
  std::vector<float> m_res_y_hit;   /// hit residual y
  std::vector<float> m_err_x_hit;   /// hit err x
  std::vector<float> m_err_y_hit;   /// hit err y
  std::vector<float> m_pull_x_hit;  /// hit pull x
  std::vector<float> m_pull_y_hit;  /// hit pull y
  std::vector<int> m_dim_hit;       /// dimension of measurement

  bool m_hasFittedParams{false};  /// if the track has fitted parameter
  float m_eLOC0_fit{-99.};        /// fitted parameter eLOC_0
  float m_eLOC1_fit{-99.};        /// fitted parameter eLOC_1
  float m_ePHI_fit{-99.};         /// fitted parameter ePHI
  float m_eTHETA_fit{-99.};       /// fitted parameter eTHETA
  float m_eQOP_fit{-99.};         /// fitted parameter eQOP
  float m_eT_fit{-99.};           /// fitted parameter eT
  float m_err_eLOC0_fit{-99.};    /// fitted parameter eLOC_-99.err
  float m_err_eLOC1_fit{-99.};    /// fitted parameter eLOC_1 err
  float m_err_ePHI_fit{-99.};     /// fitted parameter ePHI err
  float m_err_eTHETA_fit{-99.};   /// fitted parameter eTHETA err
  float m_err_eQOP_fit{-99.};     /// fitted parameter eQOP err
  float m_err_eT_fit{-99.};       /// fitted parameter eT err
  float m_px_fit{-99.};           /// fitted parameter global px
  float m_py_fit{-99.};           /// fitted parameter global py
  float m_pz_fit{-99.};           /// fitted parameter global pz
  float m_x_fit{-99.};            /// fitted parameter global PCA x
  float m_y_fit{-99.};            /// fitted parameter global PCA y
  float m_z_fit{-99.};            /// fitted parameter global PCA z
  float m_chi2_fit{-99.};         /// fitted parameter chi2
  float m_quality{-99.};          /// SvtxTrack quality parameter
  float m_ndf_fit{-99.};          /// fitted parameter ndf
  float m_dca3Dxy{-99.};          /// fitted parameter 3D DCA in xy plane
  float m_dca3Dz{-99.};           /// fitted parameter 3D DCA in z plane
  float m_dca3DxyCov{-99.};       /// fitted parameter 3D DCA covariance in xy
  float m_dca3DzCov{-99.};        /// fitted parameter 3D DCA covariance in z
  int m_charge_fit{-99};          /// fitted parameter charge

  int m_nPredicted{0};                   /// number of states with predicted parameter
  std::vector<bool> m_prt;               /// predicted status
  std::vector<float> m_eLOC0_prt;        /// predicted parameter eLOC0
  std::vector<float> m_eLOC1_prt;        /// predicted parameter eLOC1
  std::vector<float> m_ePHI_prt;         /// predicted parameter ePHI
  std::vector<float> m_eTHETA_prt;       /// predicted parameter eTHETA
  std::vector<float> m_eQOP_prt;         /// predicted parameter eQOP
  std::vector<float> m_eT_prt;           /// predicted parameter eT
  std::vector<float> m_res_eLOC0_prt;    /// predicted parameter eLOC0 residual
  std::vector<float> m_res_eLOC1_prt;    /// predicted parameter eLOC1 residual
  std::vector<float> m_res_ePHI_prt;     /// predicted parameter ePHI residual
  std::vector<float> m_res_eTHETA_prt;   /// predicted parameter eTHETA residual
  std::vector<float> m_res_eQOP_prt;     /// predicted parameter eQOP residual
  std::vector<float> m_res_eT_prt;       /// predicted parameter eT residual
  std::vector<float> m_err_eLOC0_prt;    /// predicted parameter eLOC0 error
  std::vector<float> m_err_eLOC1_prt;    /// predicted parameter eLOC1 error
  std::vector<float> m_err_ePHI_prt;     /// predicted parameter ePHI error
  std::vector<float> m_err_eTHETA_prt;   /// predicted parameter eTHETA error
  std::vector<float> m_err_eQOP_prt;     /// predicted parameter eQOP error
  std::vector<float> m_err_eT_prt;       /// predicted parameter eT error
  std::vector<float> m_pull_eLOC0_prt;   /// predicted parameter eLOC0 pull
  std::vector<float> m_pull_eLOC1_prt;   /// predicted parameter eLOC1 pull
  std::vector<float> m_pull_ePHI_prt;    /// predicted parameter ePHI pull
  std::vector<float> m_pull_eTHETA_prt;  /// predicted parameter eTHETA pull
  std::vector<float> m_pull_eQOP_prt;    /// predicted parameter eQOP pull
  std::vector<float> m_pull_eT_prt;      /// predicted parameter eT pull
  std::vector<float> m_x_prt;            /// predicted global x
  std::vector<float> m_y_prt;            /// predicted global y
  std::vector<float> m_z_prt;            /// predicted global z
  std::vector<float> m_px_prt;           /// predicted momentum px
  std::vector<float> m_py_prt;           /// predicted momentum py
  std::vector<float> m_pz_prt;           /// predicted momentum pz
  std::vector<float> m_eta_prt;          /// predicted momentum eta
  std::vector<float> m_pT_prt;           /// predicted momentum pT

  int m_nFiltered{0};                    /// number of states with filtered parameter
  std::vector<bool> m_flt;               /// filtered status
  std::vector<float> m_eLOC0_flt;        /// filtered parameter eLOC0
  std::vector<float> m_eLOC1_flt;        /// filtered parameter eLOC1
  std::vector<float> m_ePHI_flt;         /// filtered parameter ePHI
  std::vector<float> m_eTHETA_flt;       /// filtered parameter eTHETA
  std::vector<float> m_eQOP_flt;         /// filtered parameter eQOP
  std::vector<float> m_eT_flt;           /// filtered parameter eT
  std::vector<float> m_res_eLOC0_flt;    /// filtered parameter eLOC0 residual
  std::vector<float> m_res_eLOC1_flt;    /// filtered parameter eLOC1 residual
  std::vector<float> m_res_ePHI_flt;     /// filtered parameter ePHI residual
  std::vector<float> m_res_eTHETA_flt;   /// filtered parameter eTHETA residual
  std::vector<float> m_res_eQOP_flt;     /// filtered parameter eQOP residual
  std::vector<float> m_res_eT_flt;       /// filtered parameter eT residual
  std::vector<float> m_err_eLOC0_flt;    /// filtered parameter eLOC0 error
  std::vector<float> m_err_eLOC1_flt;    /// filtered parameter eLOC1 error
  std::vector<float> m_err_ePHI_flt;     /// filtered parameter ePHI error
  std::vector<float> m_err_eTHETA_flt;   /// filtered parameter eTHETA error
  std::vector<float> m_err_eQOP_flt;     /// filtered parameter eQOP error
  std::vector<float> m_err_eT_flt;       /// filtered parameter eT error
  std::vector<float> m_pull_eLOC0_flt;   /// filtered parameter eLOC0 pull
  std::vector<float> m_pull_eLOC1_flt;   /// filtered parameter eLOC1 pull
  std::vector<float> m_pull_ePHI_flt;    /// filtered parameter ePHI pull
  std::vector<float> m_pull_eTHETA_flt;  /// filtered parameter eTHETA pull
  std::vector<float> m_pull_eQOP_flt;    /// filtered parameter eQOP pull
  std::vector<float> m_pull_eT_flt;      /// filtered parameter eT pull
  std::vector<float> m_x_flt;            /// filtered global x
  std::vector<float> m_y_flt;            /// filtered global y
  std::vector<float> m_z_flt;            /// filtered global z
  std::vector<float> m_px_flt;           /// filtered momentum px
  std::vector<float> m_py_flt;           /// filtered momentum py
  std::vector<float> m_pz_flt;           /// filtered momentum pz
  std::vector<float> m_eta_flt;          /// filtered momentum eta
  std::vector<float> m_pT_flt;           /// filtered momentum pT
  std::vector<float> m_chi2;             /// chisq from filtering

  int m_nSmoothed{0};                    /// number of states with smoothed parameter
  std::vector<bool> m_smt;               /// smoothed status
  std::vector<float> m_eLOC0_smt;        /// smoothed parameter eLOC0
  std::vector<float> m_eLOC1_smt;        /// smoothed parameter eLOC1
  std::vector<float> m_ePHI_smt;         /// smoothed parameter ePHI
  std::vector<float> m_eTHETA_smt;       /// smoothed parameter eTHETA
  std::vector<float> m_eQOP_smt;         /// smoothed parameter eQOP
  std::vector<float> m_eT_smt;           /// smoothed parameter eT
  std::vector<float> m_res_eLOC0_smt;    /// smoothed parameter eLOC0 residual
  std::vector<float> m_res_eLOC1_smt;    /// smoothed parameter eLOC1 residual
  std::vector<float> m_res_ePHI_smt;     /// smoothed parameter ePHI residual
  std::vector<float> m_res_eTHETA_smt;   /// smoothed parameter eTHETA residual
  std::vector<float> m_res_eQOP_smt;     /// smoothed parameter eQOP residual
  std::vector<float> m_res_eT_smt;       /// smoothed parameter eT residual
  std::vector<float> m_err_eLOC0_smt;    /// smoothed parameter eLOC0 error
  std::vector<float> m_err_eLOC1_smt;    /// smoothed parameter eLOC1 error
  std::vector<float> m_err_ePHI_smt;     /// smoothed parameter ePHI error
  std::vector<float> m_err_eTHETA_smt;   /// smoothed parameter eTHETA error
  std::vector<float> m_err_eQOP_smt;     /// smoothed parameter eQOP error
  std::vector<float> m_err_eT_smt;       /// smoothed parameter eT error
  std::vector<float> m_pull_eLOC0_smt;   /// smoothed parameter eLOC0 pull
  std::vector<float> m_pull_eLOC1_smt;   /// smoothed parameter eLOC1 pull
  std::vector<float> m_pull_ePHI_smt;    /// smoothed parameter ePHI pull
  std::vector<float> m_pull_eTHETA_smt;  /// smoothed parameter eTHETA pull
  std::vector<float> m_pull_eQOP_smt;    /// smoothed parameter eQOP pull
  std::vector<float> m_pull_eT_smt;      /// smoothed parameter eT pull
  std::vector<float> m_x_smt;            /// smoothed global x
  std::vector<float> m_y_smt;            /// smoothed global y
  std::vector<float> m_z_smt;            /// smoothed global z
  std::vector<float> m_px_smt;           /// smoothed momentum px
  std::vector<float> m_py_smt;           /// smoothed momentum py
  std::vector<float> m_pz_smt;           /// smoothed momentum pz
  std::vector<float> m_eta_smt;          /// smoothed momentum eta
  std::vector<float> m_pT_smt;           /// smoothed momentum pT

  float m_protoTrackPx{-9999.};          /// Proto track px
  float m_protoTrackPy{-9999.};          /// Proto track py
  float m_protoTrackPz{-9999.};          /// Proto track pz
  float m_protoTrackX{-9999.};           /// Proto track PCA x
  float m_protoTrackY{-9999.};           /// Proto track PCA y
  float m_protoTrackZ{-9999.};           /// Proto track PCA z
  float m_protoD0Cov{-9999.};            /// Proto track loc0 covariance
  float m_protoZ0Cov{-9999.};            /// Proto track loc1 covariance
  float m_protoPhiCov{-9999.};           /// Proto track phi covariance
  float m_protoThetaCov{-9999.};         /// Proto track theta covariance
  float m_protoQopCov{-9999.};           /// Proto track q/p covariance
  
  
  std::vector<float> m_SL_lx;            /// Proto track source link local x pos
  std::vector<float> m_SL_ly;            /// Proto track source link local y pos
  std::vector<float> m_SLx;              /// Proto track source link global x pos
  std::vector<float> m_SLy;              /// Proto track source link global y pos
  std::vector<float> m_SLz;              /// Proto track source link global z pos
  std::vector<float> m_t_SL_lx;          /// Proto track truth hit local x
  std::vector<float> m_t_SL_ly;          /// Proto track truth hit local y
  std::vector<float> m_t_SL_gx;          /// Proto track truth hit global x
  std::vector<float> m_t_SL_gy;          /// Proto track truth hit global y
  std::vector<float> m_t_SL_gz;          /// Proto track truth hit global z

};

#endif
