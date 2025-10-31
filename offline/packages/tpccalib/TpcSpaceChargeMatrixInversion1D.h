#ifndef TPCCALIB_TPCSPACECHARGEMATRIXINVERSION1D_H
#define TPCCALIB_TPCSPACECHARGEMATRIXINVERSION1D_H
/**
 * \file TpcSpaceChargeMatrixInversion1D.h
 * \brief aggregate space charge reconstruction matrices from several jobs and perform the inversion 1D
 * \author Xudong Yu <xyu3@bnl.gov>
 */

#include "TpcSpaceChargeMatrixContainer.h"

#include <fun4all/Fun4AllBase.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <memory>
#include <map>

class TH1;

/**
 * \class TpcSpaceChargeMatrixInversion1D
 * \brief performs space charge distortion reconstruction using tracks
 */

class TpcSpaceChargeMatrixInversion1D : public Fun4AllBase
{
 public:
  /// constructor
  TpcSpaceChargeMatrixInversion1D(const std::string& = "TPCSPACECHARGEMATRIXINVERSION1D");

  ///@name modifiers
  //@{

  /// load central membrane distortion correction
  void load_cm_distortion_corrections(const std::string& /*filename*/);

  /// load average distortion correction
  void load_average_distortion_corrections(const std::string& /*filename*/);

  /// add space charge correction matrix to current. Returns true on success
  bool add(const TpcSpaceChargeMatrixContainer&, const TpcSpaceChargeMatrixContainer&, const TpcSpaceChargeMatrixContainer&, const TpcSpaceChargeMatrixContainer&);

  /// add space charge correction matrix, loaded from file, to current. Returns true on success
  bool add_from_file(const std::string& /*filename*/, const std::string& /*objectname*/ = "TpcSpaceChargeMatrixContainer");

  /// calculate distortions by inverting stored matrices, and save relevant histograms
  void calculate_distortion_corrections();

  /// extrapolate distortions
  void extrapolate_distortion_corrections();

  /// save distortions
  void save_distortion_corrections(const std::string& /*filename*/ = "DistortionCorrections.root");

  TH1* transform_dphi_rdphi_layer(const TH1* source, const TString& name);
  TH1* transform_dphi_rdphi_radius(const TH1* source, const TString& name);

  //@}

 private:
  /// matrix container
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container_layer_negz;
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container_layer_posz;
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container_radius_negz;
  std::unique_ptr<TpcSpaceChargeMatrixContainer> m_matrix_container_radius_posz;

  /// output distortion container
  std::unique_ptr<TpcDistortionCorrectionContainer> m_dcc_average_layer;
  std::unique_ptr<TpcDistortionCorrectionContainer> m_dcc_average_radius;

  /// central membrane distortion container
  std::unique_ptr<TpcDistortionCorrectionContainer> m_dcc_cm;

  TH1* h_rdphi_layer_negz = nullptr;
  TH1* h_rdphi_layer_posz = nullptr;
  TH1* h_rdphi_radius_negz = nullptr;
  TH1* h_rdphi_radius_posz = nullptr;

  std::map<int, float> TpcRadiusMap;
  void loadTpcRadius()
  {
  TpcRadiusMap[6] = 31.372; // for interpolate
  TpcRadiusMap[7] = 31.372;
  TpcRadiusMap[8] = 31.944;
  TpcRadiusMap[9] = 32.5159;
  TpcRadiusMap[10] = 33.0879;
  TpcRadiusMap[11] = 33.6599;
  TpcRadiusMap[12] = 34.2318;
  TpcRadiusMap[13] = 34.8038;
  TpcRadiusMap[14] = 35.3758;
  TpcRadiusMap[15] = 35.9477;
  TpcRadiusMap[16] = 36.5197;
  TpcRadiusMap[17] = 37.0917;
  TpcRadiusMap[18] = 37.6636;
  TpcRadiusMap[19] = 38.2356;
  TpcRadiusMap[20] = 38.8076;
  TpcRadiusMap[21] = 39.3795;
  TpcRadiusMap[22] = 39.9515;
  TpcRadiusMap[23] = 41.659;
  TpcRadiusMap[24] = 42.6797;
  TpcRadiusMap[25] = 43.7003;
  TpcRadiusMap[26] = 44.721;
  TpcRadiusMap[27] = 45.7417;
  TpcRadiusMap[28] = 46.7623;
  TpcRadiusMap[29] = 47.783;
  TpcRadiusMap[30] = 48.8037;
  TpcRadiusMap[31] = 49.8243;
  TpcRadiusMap[32] = 50.845;
  TpcRadiusMap[33] = 51.8657;
  TpcRadiusMap[34] = 52.8863;
  TpcRadiusMap[35] = 53.907;
  TpcRadiusMap[36] = 54.9277;
  TpcRadiusMap[37] = 55.9483;
  TpcRadiusMap[38] = 56.969;
  TpcRadiusMap[39] = 58.911;
  TpcRadiusMap[40] = 60.0081;
  TpcRadiusMap[41] = 61.1051;
  TpcRadiusMap[42] = 62.2022;
  TpcRadiusMap[43] = 63.2993;
  TpcRadiusMap[44] = 64.3963;
  TpcRadiusMap[45] = 65.4934;
  TpcRadiusMap[46] = 66.5905;
  TpcRadiusMap[47] = 67.6875;
  TpcRadiusMap[48] = 68.7846;
  TpcRadiusMap[49] = 69.8817;
  TpcRadiusMap[50] = 70.9787;
  TpcRadiusMap[51] = 72.0758;
  TpcRadiusMap[52] = 73.1729;
  TpcRadiusMap[53] = 74.2699;
  TpcRadiusMap[54] = 75.367;
  TpcRadiusMap[55] = 75.367; // for interpolate
  }

};

#endif
