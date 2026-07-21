// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCTRACKFIT_H
#define TPCTRACKRECO_TPCTRACKFIT_H

#include <array>
#include <string>
#include <vector>

class PHField;

struct TpcTrackVec3
{
  double x{0.0};
  double y{0.0};
  double z{0.0};
};

struct TpcTrackPoint
{
  int track_id{0};
  int shower_id{0};
  int layer{0};
  TpcTrackVec3 position;
  TpcTrackVec3 momentum;
  double t{0.0};
  double path{0.0};
};

struct TpcTrackHelix
{
  double cx{0.0};
  double cy{0.0};
  double radius{0.0};
  double z0{0.0};
  double pitch{0.0};
  double theta_first{0.0};
  double theta_last{0.0};
  double theta_min{0.0};
  double theta_max{0.0};
  double direction{1.0};
  double bfield_t{1.4};
};

struct TpcTrackLinePca
{
  TpcTrackVec3 pca1;
  TpcTrackVec3 pca2;
  double dca{0.0};
  double step1{0.0};
  double step2{0.0};
};

struct TpcTrackHelixPca
{
  TpcTrackVec3 pca1;
  TpcTrackVec3 pca2;
  double dca{0.0};
  double theta1{0.0};
  double theta2{0.0};
};

struct TpcTrackHelixSearchRange
{
  bool valid{false};
  int anchor_point_index{-1};
  double anchor_theta{0.0};
  double anchor_path_cm{0.0};
  double anchor_residual_cm{0.0};
  double theta_min{0.0};
  double theta_max{0.0};
  double upstream_cm{0.0};
  double downstream_cm{0.0};
};

enum class TpcTrackPointOrder
{
  Path,
  Input,
  Radius,
  ThetaZ,
  Auto
};

struct TpcTrackState
{
  TpcTrackVec3 position;
  TpcTrackVec3 momentum;
  int charge{0};
  double chi2{0.0};
  int ndof{0};
  bool valid{false};
};

struct TpcKalmanConfig
{
  double bfield_t{1.4};
  const PHField *magnetic_field{nullptr};
  bool analytic_uniform_propagation{false};
  double rkn_max_step_cm{5.0};
  double rkn_step_tolerance{1.0e-4};
  int rkn_max_step_trials{12};
  int rkn_max_total_steps{2000};
  bool rkn_fast_field_jacobian{true};
  bool rkn_fast_field_pca{true};
  int rkn_field_pca_refine_iterations{6};
  TpcTrackPointOrder point_order{TpcTrackPointOrder::Radius};
  double meas_sigma_rphi_cm{0.03};
  double meas_sigma_r_cm{0.03};
  double meas_sigma_z_cm{0.05};
  double min_measurement_sigma_cm{1.0e-6};
  double initial_sigma_pos_cm{0.1};
  double initial_sigma_phi{0.2};
  double initial_sigma_qop_t{0.2};
  double initial_sigma_tanl{0.2};
  double process_sigma_pos_cm{1.0e-4};
  double process_sigma_phi{1.0e-5};
  double process_sigma_qop_t{1.0e-6};
  double process_sigma_tanl{1.0e-6};
  bool collect_innovation_components{false};
  double material_x0_per_cm{0.0};
  double multiple_scattering_scale{1.0};
  double energy_loss_gev_per_cm{0.0};
  double energy_loss_sigma_fraction{0.0};
  double min_pt_gev{0.05};
};

struct TpcKalmanResult
{
  bool success{false};
  std::string message;
  int charge{0};
  double bfield_t{1.4};
  const PHField *magnetic_field{nullptr};
  bool analytic_uniform_propagation{false};
  TpcTrackHelix seed;
  std::vector<double> path_s;
  std::vector<std::array<double, 6>> states_filtered;
  std::vector<std::array<double, 36>> covs_filtered;
  std::vector<std::array<double, 6>> states_smoothed;
  std::vector<std::array<double, 36>> covs_smoothed;
  std::vector<double> measurement_chi2;
  std::vector<unsigned char> measurement_used;
  std::vector<unsigned char> measurement_in_seed;
  std::vector<double> innovation_residual_r;
  std::vector<double> innovation_residual_rphi;
  std::vector<double> innovation_residual_z;
  std::vector<double> prediction_sigma_r;
  std::vector<double> prediction_sigma_rphi;
  std::vector<double> prediction_sigma_z;
  std::vector<double> innovation_sigma_r;
  std::vector<double> innovation_sigma_rphi;
  std::vector<double> innovation_sigma_z;
  std::vector<double> innovation_rho_r_rphi;
  std::vector<double> innovation_rho_r_z;
  std::vector<double> innovation_rho_rphi_z;
  std::vector<double> innovation_whitened_0;
  std::vector<double> innovation_whitened_1;
  std::vector<double> innovation_whitened_2;
  std::size_t naccepted{0};
  std::size_t nrejected{0};
  double chi2{0.0};
  int ndof{0};
  double mass_gev{0.13957039};
  std::size_t rkn_propagations{0};
  std::size_t rkn_accepted_steps{0};
  std::size_t rkn_rejected_trials{0};
  std::size_t rkn_max_trial_accepts{0};
  std::size_t rkn_failures{0};
  double rkn_seconds{0.0};
};

#endif
