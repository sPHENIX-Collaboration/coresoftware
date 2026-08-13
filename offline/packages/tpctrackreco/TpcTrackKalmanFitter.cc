#include "TpcTrackKalmanFitter.h"

#include "TpcTrackHelixFitter.h"

#include <phfield/PHField.h>

#include <CLHEP/Units/SystemOfUnits.h>

#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>

namespace
{
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kCurvaturePerCm = 0.003;

  template <class T>
  constexpr T square(const T &value)
  {
    return value * value;
  }

  double normalize_phi(const double phi)
  {
    return std::atan2(std::sin(phi), std::cos(phi));
  }

  using StateVector = Eigen::Matrix<double, 6, 1>;
  using StateMatrix = Eigen::Matrix<double, 6, 6>;
  using Vector3 = Eigen::Matrix<double, 3, 1>;
  using MeasurementVector = Eigen::Matrix<double, 3, 1>;
  using MeasurementMatrix = Eigen::Matrix<double, 3, 6>;
  using MeasurementCov = Eigen::Matrix<double, 3, 3>;

  StateVector to_eigen(const std::array<double, 6> &state)
  {
    StateVector output;
    for (int i = 0; i < 6; ++i)
    {
      output(i) = state[static_cast<std::size_t>(i)];
    }
    return output;
  }

  std::array<double, 6> to_array(const StateVector &state)
  {
    std::array<double, 6> output{};
    for (int i = 0; i < 6; ++i)
    {
      output[static_cast<std::size_t>(i)] = state(i);
    }
    return output;
  }

  std::array<double, 36> to_array(const StateMatrix &matrix)
  {
    std::array<double, 36> output{};
    for (int row = 0; row < 6; ++row)
    {
      for (int col = 0; col < 6; ++col)
      {
        const auto offset = static_cast<std::size_t>(row) * 6U + static_cast<std::size_t>(col);
        output[offset] = matrix(row, col);
      }
    }
    return output;
  }

  StateVector residual(const StateVector &lhs, const StateVector &rhs)
  {
    StateVector diff = lhs - rhs;
    diff(TpcTrackKalmanFitter::Phi) = normalize_phi(diff(TpcTrackKalmanFitter::Phi));
    return diff;
  }

  double omega_from_state(const StateVector &state, const double bfield_t)
  {
    return kCurvaturePerCm * bfield_t * state(TpcTrackKalmanFitter::QOverPt);
  }

  StateVector apply_mean_energy_loss_path3d(const StateVector &state,
                                            const double signed_path3d_cm,
                                            const TpcKalmanConfig &config,
                                            const double mass_gev)
  {
    if (config.energy_loss_gev_per_cm <= 0.0 || signed_path3d_cm == 0.0)
    {
      return state;
    }

    StateVector output = state;
    const double qop_t = output(TpcTrackKalmanFitter::QOverPt);
    if (std::abs(qop_t) < 1.0e-12)
    {
      return output;
    }

    const double tanl = output(TpcTrackKalmanFitter::TanLambda);
    const double path3d_cm = std::abs(signed_path3d_cm);
    if (path3d_cm <= 0.0)
    {
      return output;
    }

    const double pt = 1.0 / std::abs(qop_t);
    const double momentum = pt * std::sqrt(1.0 + tanl * tanl);
    const double energy = std::sqrt(momentum * momentum + mass_gev * mass_gev);
    const double signed_loss = std::copysign(config.energy_loss_gev_per_cm * path3d_cm,
                                             signed_path3d_cm);
    const double new_energy = std::max(mass_gev + 1.0e-9, energy - signed_loss);
    const double new_momentum = std::sqrt(std::max(0.0, new_energy * new_energy - mass_gev * mass_gev));
    if (new_momentum <= 0.0)
    {
      return output;
    }

    const double min_pt = std::max(config.min_pt_gev, 1.0e-6);
    const double new_pt = std::max(min_pt, new_momentum / std::sqrt(1.0 + tanl * tanl));
    output(TpcTrackKalmanFitter::QOverPt) = std::copysign(1.0 / new_pt, qop_t);
    return output;
  }

  StateVector apply_mean_energy_loss(const StateVector &state,
                                     const double ds_cm,
                                     const TpcKalmanConfig &config,
                                     const double mass_gev)
  {
    const double tanl = state(TpcTrackKalmanFitter::TanLambda);
    return apply_mean_energy_loss_path3d(state, ds_cm * std::sqrt(1.0 + tanl * tanl),
                                         config, mass_gev);
  }

  Vector3 direction_from_state(const StateVector &state)
  {
    Vector3 direction(std::cos(state(TpcTrackKalmanFitter::Phi)),
                      std::sin(state(TpcTrackKalmanFitter::Phi)),
                      state(TpcTrackKalmanFitter::TanLambda));
    const double norm = direction.norm();
    if (norm <= 0.0 || !std::isfinite(norm))
    {
      return Vector3(1.0, 0.0, 0.0);
    }
    return direction / norm;
  }

  Vector3 magnetic_field_tesla(const Vector3 &position_cm,
                               const TpcKalmanConfig &config)
  {
    if (config.magnetic_field != nullptr)
    {
      const double point[4] = {position_cm.x() * CLHEP::cm,
                               position_cm.y() * CLHEP::cm,
                               position_cm.z() * CLHEP::cm,
                               0.0};
      double field[3] = {0.0, 0.0, 0.0};
      config.magnetic_field->GetFieldValue(point, field);
      Vector3 output(field[0] / CLHEP::tesla,
                     field[1] / CLHEP::tesla,
                     field[2] / CLHEP::tesla);
      if (std::isfinite(output.x()) && std::isfinite(output.y()) &&
          std::isfinite(output.z()))
      {
        return output;
      }
      return Vector3::Zero();
    }

    return Vector3(0.0, 0.0, config.bfield_t);
  }

  Vector3 rkn_acceleration(const Vector3 &direction,
                           const Vector3 &bfield_tesla,
                           const double q_over_p)
  {
    return kCurvaturePerCm * q_over_p * direction.cross(bfield_tesla);
  }

  struct RknStep
  {
    Vector3 position;
    Vector3 direction;
    double error{0.0};
  };

  struct RknDiagnostics
  {
    std::size_t propagations{0};
    std::size_t accepted_steps{0};
    std::size_t rejected_trials{0};
    std::size_t max_trial_accepts{0};
    std::size_t failures{0};
    double seconds{0.0};
  };

  RknStep try_rkn4_step(const Vector3 &position_cm,
                        const Vector3 &direction,
                        const double q_over_p,
                        const double h_cm,
                        const TpcKalmanConfig &config)
  {
    const double h2 = h_cm * h_cm;
    const double half_h = 0.5 * h_cm;

    const Vector3 b_first = magnetic_field_tesla(position_cm, config);
    const Vector3 k1 = rkn_acceleration(direction, b_first, q_over_p);

    const Vector3 pos1 = position_cm + half_h * direction + 0.125 * h2 * k1;
    const Vector3 b_middle = magnetic_field_tesla(pos1, config);
    const Vector3 k2 = rkn_acceleration(direction + half_h * k1, b_middle, q_over_p);
    const Vector3 k3 = rkn_acceleration(direction + half_h * k2, b_middle, q_over_p);

    const Vector3 pos2 = position_cm + h_cm * direction + 0.5 * h2 * k3;
    const Vector3 b_last = magnetic_field_tesla(pos2, config);
    const Vector3 k4 = rkn_acceleration(direction + h_cm * k3, b_last, q_over_p);

    RknStep output;
    output.position = position_cm + h_cm * direction + h2 / 6.0 * (k1 + k2 + k3);
    output.direction = direction + h_cm / 6.0 * (k1 + 2.0 * (k2 + k3) + k4);
    const double direction_norm = output.direction.norm();
    if (direction_norm > 0.0 && std::isfinite(direction_norm))
    {
      output.direction /= direction_norm;
    }
    else
    {
      output.direction = direction;
    }
    output.error = std::max(1.0e-20,
                            h2 * (k1 - k2 - k3 + k4).template lpNorm<1>());
    return output;
  }

  double rkn_step_scale(const double tolerance, const double error)
  {
    if (error <= 0.0 || !std::isfinite(error))
    {
      return 1.0;
    }
    return std::clamp(std::sqrt(std::sqrt(tolerance / error)), 0.25, 4.0);
  }

  bool advance_rkn4(Vector3 &position_cm,
                    Vector3 &direction,
                    const double q_over_p,
                    const double signed_path3d_cm,
                    const TpcKalmanConfig &config,
                    RknDiagnostics *diagnostics)
  {
    double remaining = signed_path3d_cm;
    if (remaining == 0.0)
    {
      return true;
    }

    constexpr double min_step_cm = 1.0e-5;
    const std::size_t max_total_steps =
        static_cast<std::size_t>(std::max(config.rkn_max_total_steps, 1));
    const double max_step_cm = std::max(std::abs(config.rkn_max_step_cm), min_step_cm);
    const double tolerance = std::max(config.rkn_step_tolerance, 1.0e-20);
    const int max_trials = std::max(config.rkn_max_step_trials, 1);
    const bool adaptive = config.rkn_step_tolerance > 0.0;

    double h = std::copysign(std::min(std::abs(remaining), max_step_cm), remaining);
    for (std::size_t nsteps = 0;
         std::abs(remaining) > min_step_cm && nsteps < max_total_steps;
         ++nsteps)
    {
      if (std::abs(h) > std::abs(remaining))
      {
        h = remaining;
      }

      RknStep trial;
      int ntrials = 0;
      while (true)
      {
        trial = try_rkn4_step(position_cm, direction, q_over_p, h, config);
        if (!trial.position.allFinite() || !trial.direction.allFinite() ||
            !std::isfinite(trial.error))
        {
          if (diagnostics != nullptr)
          {
            ++diagnostics->failures;
          }
          return false;
        }
        if (!adaptive || trial.error <= 4.0 * tolerance ||
            std::abs(h) <= min_step_cm || ntrials >= max_trials)
        {
          if (adaptive && trial.error > 4.0 * tolerance && ntrials >= max_trials &&
              diagnostics != nullptr)
          {
            ++diagnostics->max_trial_accepts;
          }
          break;
        }
        h *= rkn_step_scale(tolerance, trial.error);
        ++ntrials;
        if (diagnostics != nullptr)
        {
          ++diagnostics->rejected_trials;
        }
      }

      position_cm = trial.position;
      direction = trial.direction;
      remaining -= h;
      if (diagnostics != nullptr)
      {
        ++diagnostics->accepted_steps;
      }

      if (std::abs(remaining) <= min_step_cm)
      {
        break;
      }

      double next_h = std::abs(h);
      if (adaptive)
      {
        next_h *= rkn_step_scale(tolerance, trial.error);
      }
      next_h = std::min({next_h, max_step_cm, std::abs(remaining)});
      h = std::copysign(std::max(next_h, min_step_cm), remaining);
    }

    const bool complete = std::abs(remaining) <= min_step_cm;
    if (!complete && diagnostics != nullptr)
    {
      ++diagnostics->failures;
    }
    return complete;
  }

  StateVector propagate_rkn4(const StateVector &state,
                             const double ds_cm,
                             const TpcKalmanConfig &config,
                             const double mass_gev,
                             RknDiagnostics *diagnostics = nullptr)
  {
    const auto propagation_start = std::chrono::steady_clock::now();
    const auto record_seconds = [&]()
    {
      if (diagnostics != nullptr)
      {
        diagnostics->seconds += std::chrono::duration<double>(
                                    std::chrono::steady_clock::now() - propagation_start)
                                    .count();
      }
    };
    if (diagnostics != nullptr)
    {
      ++diagnostics->propagations;
    }
    StateVector output = state;
    if (ds_cm == 0.0)
    {
      output = apply_mean_energy_loss(output, ds_cm, config, mass_gev);
      record_seconds();
      return output;
    }

    Vector3 position_cm(state(TpcTrackKalmanFitter::X),
                        state(TpcTrackKalmanFitter::Y),
                        state(TpcTrackKalmanFitter::Z));
    Vector3 direction = direction_from_state(state);
    const double transverse = std::hypot(direction.x(), direction.y());
    if (transverse <= 1.0e-12 || !std::isfinite(transverse))
    {
      output = apply_mean_energy_loss(output, ds_cm, config, mass_gev);
      record_seconds();
      return output;
    }

    const double q_over_p = -state(TpcTrackKalmanFitter::QOverPt) * transverse;
    const double signed_path3d_cm = ds_cm / transverse;
    if (!advance_rkn4(position_cm, direction, q_over_p, signed_path3d_cm,
                      config, diagnostics))
    {
      record_seconds();
      return StateVector::Constant(std::numeric_limits<double>::quiet_NaN());
    }

    const double final_transverse = std::hypot(direction.x(), direction.y());
    output(TpcTrackKalmanFitter::X) = position_cm.x();
    output(TpcTrackKalmanFitter::Y) = position_cm.y();
    output(TpcTrackKalmanFitter::Z) = position_cm.z();
    if (final_transverse > 1.0e-12 && std::isfinite(final_transverse))
    {
      output(TpcTrackKalmanFitter::Phi) =
          normalize_phi(std::atan2(direction.y(), direction.x()));
      output(TpcTrackKalmanFitter::TanLambda) = direction.z() / final_transverse;
      output(TpcTrackKalmanFitter::QOverPt) = -q_over_p / final_transverse;
      const double max_abs_qop_t = 1.0 / std::max(config.min_pt_gev, 1.0e-6);
      if (std::abs(output(TpcTrackKalmanFitter::QOverPt)) > max_abs_qop_t)
      {
        output(TpcTrackKalmanFitter::QOverPt) =
            std::copysign(max_abs_qop_t, output(TpcTrackKalmanFitter::QOverPt));
      }
    }

    output = apply_mean_energy_loss_path3d(output, signed_path3d_cm, config, mass_gev);
    record_seconds();
    return output;
  }

  StateVector propagate_uniform_bz(const StateVector &state,
                                   const double ds_cm,
                                   const double bfield_t,
                                   const TpcKalmanConfig &config,
                                   const double mass_gev)
  {
    StateVector output = state;
    const double phi = state(TpcTrackKalmanFitter::Phi);
    const double omega = kCurvaturePerCm * bfield_t *
                         state(TpcTrackKalmanFitter::QOverPt);

    if (std::abs(omega) < 1.0e-12)
    {
      output(TpcTrackKalmanFitter::X) += ds_cm * std::cos(phi);
      output(TpcTrackKalmanFitter::Y) += ds_cm * std::sin(phi);
    }
    else
    {
      const double phi2 = phi + omega * ds_cm;
      output(TpcTrackKalmanFitter::X) +=
          (std::sin(phi2) - std::sin(phi)) / omega;
      output(TpcTrackKalmanFitter::Y) +=
          -(std::cos(phi2) - std::cos(phi)) / omega;
      output(TpcTrackKalmanFitter::Phi) = normalize_phi(phi2);
    }
    output(TpcTrackKalmanFitter::Z) +=
        state(TpcTrackKalmanFitter::TanLambda) * ds_cm;
    return apply_mean_energy_loss(output, ds_cm, config, mass_gev);
  }

  StateMatrix transport_jacobian(const StateVector &state,
                                 const double ds_cm,
                                 const TpcKalmanConfig &config,
                                 const double mass_gev,
                                 RknDiagnostics *diagnostics)
  {
    StateMatrix jac = StateMatrix::Identity();
    const bool use_analytic_uniform =
        config.magnetic_field == nullptr && config.analytic_uniform_propagation;
    const bool use_fast_field_jacobian =
        config.magnetic_field != nullptr && config.rkn_fast_field_jacobian;
    double local_bz_t = config.bfield_t;
    if (use_fast_field_jacobian)
    {
      const Vector3 position_cm(state(TpcTrackKalmanFitter::X),
                                state(TpcTrackKalmanFitter::Y),
                                state(TpcTrackKalmanFitter::Z));
      local_bz_t = magnetic_field_tesla(position_cm, config).z();
      if (!std::isfinite(local_bz_t))
      {
        local_bz_t = config.bfield_t;
      }
    }

    const auto propagate_for_jacobian = [&](const StateVector &input)
    {
      if (use_analytic_uniform || use_fast_field_jacobian)
      {
        return propagate_uniform_bz(input, ds_cm, local_bz_t, config, mass_gev);
      }
      return propagate_rkn4(input, ds_cm, config, mass_gev, diagnostics);
    };
    const double scales[6] = {1.0e-4, 1.0e-4, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-6};
    for (int col = 0; col < 6; ++col)
    {
      const double step = scales[col] * std::max(1.0, std::abs(state(col)));
      StateVector plus = state;
      StateVector minus = state;
      plus(col) += step;
      minus(col) -= step;
      if (col == TpcTrackKalmanFitter::Phi)
      {
        plus(col) = normalize_phi(plus(col));
        minus(col) = normalize_phi(minus(col));
      }

      const StateVector f_plus = propagate_for_jacobian(plus);
      const StateVector f_minus = propagate_for_jacobian(minus);
      if (!f_plus.allFinite() || !f_minus.allFinite())
      {
        return StateMatrix::Constant(std::numeric_limits<double>::quiet_NaN());
      }
      jac.col(col) = residual(f_plus, f_minus) / (2.0 * step);
    }
    return jac;
  }

  double multiple_scattering_theta0(const StateVector &state,
                                    const double ds_cm,
                                    const TpcKalmanConfig &config,
                                    const double mass_gev)
  {
    if (config.material_x0_per_cm <= 0.0 || ds_cm == 0.0)
    {
      return 0.0;
    }

    const double qop_t = state(TpcTrackKalmanFitter::QOverPt);
    if (std::abs(qop_t) < 1.0e-12)
    {
      return 0.0;
    }

    const double tanl = state(TpcTrackKalmanFitter::TanLambda);
    const double path3d_cm = std::abs(ds_cm) * std::sqrt(1.0 + tanl * tanl);
    const double x_over_x0 = config.material_x0_per_cm * path3d_cm;
    if (x_over_x0 <= 0.0)
    {
      return 0.0;
    }

    const double pt = 1.0 / std::abs(qop_t);
    const double momentum = pt * std::sqrt(1.0 + tanl * tanl);
    const double energy = std::sqrt(momentum * momentum + mass_gev * mass_gev);
    const double beta = (energy > 0.0) ? momentum / energy : 0.0;
    if (beta <= 0.0 || momentum <= 0.0)
    {
      return 0.0;
    }

    const double log_term = 1.0 + 0.038 * std::log(x_over_x0);
    return config.multiple_scattering_scale * 0.0136 / (beta * momentum) *
           std::sqrt(x_over_x0) * log_term;
  }

  StateMatrix process_noise(const StateVector &state,
                            const double ds_cm,
                            const TpcKalmanConfig &config,
                            const double mass_gev)
  {
    const double scale = std::max(1.0, std::abs(ds_cm));
    StateMatrix noise = StateMatrix::Zero();
    noise(TpcTrackKalmanFitter::X, TpcTrackKalmanFitter::X) = square(config.process_sigma_pos_cm * scale);
    noise(TpcTrackKalmanFitter::Y, TpcTrackKalmanFitter::Y) = square(config.process_sigma_pos_cm * scale);
    noise(TpcTrackKalmanFitter::Z, TpcTrackKalmanFitter::Z) = square(config.process_sigma_pos_cm * scale);
    noise(TpcTrackKalmanFitter::Phi, TpcTrackKalmanFitter::Phi) = square(config.process_sigma_phi * scale);
    noise(TpcTrackKalmanFitter::QOverPt, TpcTrackKalmanFitter::QOverPt) = square(config.process_sigma_qop_t * scale);
    noise(TpcTrackKalmanFitter::TanLambda, TpcTrackKalmanFitter::TanLambda) = square(config.process_sigma_tanl * scale);

    const double theta0 = multiple_scattering_theta0(state, ds_cm, config, mass_gev);
    if (theta0 > 0.0)
    {
      const double tanl = state(TpcTrackKalmanFitter::TanLambda);
      const double sec2_lambda = 1.0 + square(tanl);
      noise(TpcTrackKalmanFitter::Phi, TpcTrackKalmanFitter::Phi) +=
          square(theta0 * std::sqrt(sec2_lambda));
      noise(TpcTrackKalmanFitter::TanLambda, TpcTrackKalmanFitter::TanLambda) +=
          square(theta0 * sec2_lambda);
    }

    if (config.energy_loss_sigma_fraction > 0.0 && config.energy_loss_gev_per_cm > 0.0)
    {
      const double tanl = state(TpcTrackKalmanFitter::TanLambda);
      const double path3d_cm = std::abs(ds_cm) * std::sqrt(1.0 + tanl * tanl);
      const double loss_sigma = config.energy_loss_sigma_fraction *
                                config.energy_loss_gev_per_cm * path3d_cm;
      const double qop_t = state(TpcTrackKalmanFitter::QOverPt);
      if (std::abs(qop_t) > 1.0e-12)
      {
        const double pt = 1.0 / std::abs(qop_t);
        const double momentum = pt * std::sqrt(1.0 + tanl * tanl);
        if (momentum > 0.0)
        {
          noise(TpcTrackKalmanFitter::QOverPt, TpcTrackKalmanFitter::QOverPt) +=
              square(qop_t * loss_sigma / momentum);
        }
      }
    }

    return noise;
  }

  MeasurementCov measurement_covariance(const TpcTrackPoint &point,
                                        const double var_rphi,
                                        const double var_r,
                                        const double var_z)
  {
    const double radius = std::hypot(point.position.x, point.position.y);
    const double cos_phi = (radius > 0.0) ? point.position.x / radius : 1.0;
    const double sin_phi = (radius > 0.0) ? point.position.y / radius : 0.0;

    MeasurementCov cov = MeasurementCov::Zero();
    cov(0, 0) = var_r * square(cos_phi) + var_rphi * square(sin_phi);
    cov(1, 1) = var_r * square(sin_phi) + var_rphi * square(cos_phi);
    cov(0, 1) = (var_r - var_rphi) * sin_phi * cos_phi;
    cov(1, 0) = cov(0, 1);
    cov(2, 2) = var_z;
    return cov;
  }

  MeasurementCov measurement_local_rotation(const TpcTrackPoint &point)
  {
    const double radius = std::hypot(point.position.x, point.position.y);
    const double cos_phi = (radius > 0.0) ? point.position.x / radius : 1.0;
    const double sin_phi = (radius > 0.0) ? point.position.y / radius : 0.0;

    MeasurementCov rotation = MeasurementCov::Zero();
    rotation(0, 0) = cos_phi;
    rotation(0, 1) = sin_phi;
    rotation(1, 0) = -sin_phi;
    rotation(1, 1) = cos_phi;
    rotation(2, 2) = 1.0;
    return rotation;
  }

  double covariance_sigma(const MeasurementCov &covariance, const int index)
  {
    return std::sqrt(std::max(0.0, covariance(index, index)));
  }

  double covariance_correlation(const MeasurementCov &covariance,
                                const int first,
                                const int second,
                                const double sigma_first,
                                const double sigma_second)
  {
    const double denominator = sigma_first * sigma_second;
    if (!(denominator > 0.0) || !std::isfinite(denominator))
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return std::clamp(covariance(first, second) / denominator, -1.0, 1.0);
  }

  MeasurementVector whiten_innovation(const MeasurementCov &innovation,
                                      const MeasurementVector &residual)
  {
    MeasurementVector whitened = MeasurementVector::Constant(
        std::numeric_limits<double>::quiet_NaN());

    const Eigen::LLT<MeasurementCov> llt(innovation);
    if (llt.info() == Eigen::Success)
    {
      whitened = llt.matrixL().solve(residual);
      return whitened;
    }

    const Eigen::SelfAdjointEigenSolver<MeasurementCov> eigensolver(innovation);
    if (eigensolver.info() != Eigen::Success ||
        eigensolver.eigenvalues().minCoeff() <= 0.0)
    {
      return whitened;
    }

    whitened = eigensolver.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal() *
               eigensolver.eigenvectors().transpose() * residual;
    return whitened;
  }

  bool make_seed(const std::vector<TpcTrackPoint> &points,
                 const TpcKalmanConfig &config,
                 TpcTrackHelix &seed,
                 std::vector<double> &theta_values,
                 std::vector<double> &path_s)
  {
    if (!TpcTrackHelixFitter::fit(points, 0, config.bfield_t, seed))
    {
      return false;
    }

    theta_values.clear();
    theta_values.reserve(points.size());
    for (const auto &point : points)
    {
      double theta = std::atan2(point.position.y - seed.cy, point.position.x - seed.cx);
      if (!theta_values.empty())
      {
        while (theta - theta_values.back() > kPi)
        {
          theta -= 2.0 * kPi;
        }
        while (theta - theta_values.back() < -kPi)
        {
          theta += 2.0 * kPi;
        }
      }
      theta_values.push_back(theta);
    }

    if (theta_values.size() < 2)
    {
      return false;
    }

    path_s.assign(theta_values.size(), 0.0);
    for (std::size_t i = 1; i < theta_values.size(); ++i)
    {
      path_s[i] = path_s[i - 1] + seed.radius * std::abs(theta_values[i] - theta_values[i - 1]);
    }
    return true;
  }

  bool initial_state(const std::vector<TpcTrackPoint> &points,
                     const TpcTrackHelix &seed,
                     const std::vector<double> &theta_values,
                     const TpcKalmanConfig &config,
                     StateVector &state)
  {
    if (points.empty() || theta_values.empty())
    {
      return false;
    }

    const double direction = seed.direction;
    const double theta0 = theta_values.front();
    const TpcTrackVec3 seed_position = TpcTrackHelixFitter::point(seed, theta0);
    if (!TpcTrackHelixFitter::finite(seed_position))
    {
      return false;
    }

    const double denom = kCurvaturePerCm * config.bfield_t * seed.radius;
    if (std::abs(denom) <= 0.0 || !std::isfinite(denom))
    {
      return false;
    }

    state.setZero();
    state(TpcTrackKalmanFitter::X) = seed_position.x;
    state(TpcTrackKalmanFitter::Y) = seed_position.y;
    state(TpcTrackKalmanFitter::Z) = seed_position.z;
    state(TpcTrackKalmanFitter::Phi) = normalize_phi(std::atan2(direction * std::cos(theta0),
                                                                direction * -std::sin(theta0)));
    state(TpcTrackKalmanFitter::QOverPt) = direction / denom;
    state(TpcTrackKalmanFitter::TanLambda) = direction * seed.pitch / seed.radius;

    const double max_abs_qop_t = 1.0 / std::max(config.min_pt_gev, 1.0e-6);
    if (std::abs(state(TpcTrackKalmanFitter::QOverPt)) > max_abs_qop_t)
    {
      state(TpcTrackKalmanFitter::QOverPt) = std::copysign(max_abs_qop_t,
                                                           state(TpcTrackKalmanFitter::QOverPt));
    }
    return true;
  }
}  // namespace

bool TpcTrackKalmanFitter::fit(const std::vector<TpcTrackPoint> &input_points,
                               const int charge,
                               const TpcKalmanConfig &config,
                               TpcKalmanResult &result,
                               const double mass_gev)
{
  result = TpcKalmanResult{};
  result.charge = charge;
  result.bfield_t = config.bfield_t;
  result.magnetic_field = config.magnetic_field;
  result.analytic_uniform_propagation = config.analytic_uniform_propagation;
  result.mass_gev = mass_gev;

  if (input_points.size() < 5)
  {
    result.message = "need at least five TPC points";
    return false;
  }

  std::vector<TpcTrackPoint> points = input_points;
  TpcTrackHelixFitter::order_points(points, config.point_order);

  std::vector<double> theta_values;
  if (!make_seed(points, config, result.seed, theta_values, result.path_s))
  {
    result.message = "helix seed failed";
    return false;
  }

  StateVector state;
  if (!initial_state(points, result.seed, theta_values, config, state))
  {
    result.message = "initial state failed";
    return false;
  }

  const std::size_t npoints = points.size();
  std::vector<StateVector> states_filtered(npoints);
  std::vector<StateVector> states_predicted(npoints);
  std::vector<StateMatrix> covs_filtered(npoints);
  std::vector<StateMatrix> covs_predicted(npoints);
  std::vector<StateMatrix> transport(npoints);

  const double min_meas = std::max(config.min_measurement_sigma_cm, 1.0e-12);
  const double sigma_rphi = std::max(config.meas_sigma_rphi_cm, min_meas);
  const double sigma_r = std::max(config.meas_sigma_r_cm, min_meas);
  const double sigma_z = std::max(config.meas_sigma_z_cm, min_meas);
  const double var_rphi = square(sigma_rphi);
  const double var_r = square(sigma_r);
  const double var_z = square(sigma_z);

  MeasurementMatrix hmat = MeasurementMatrix::Zero();
  hmat(0, X) = 1.0;
  hmat(1, Y) = 1.0;
  hmat(2, Z) = 1.0;

  StateMatrix cov = StateMatrix::Zero();
  cov(X, X) = square(config.initial_sigma_pos_cm);
  cov(Y, Y) = square(config.initial_sigma_pos_cm);
  cov(Z, Z) = square(config.initial_sigma_pos_cm);
  cov(Phi, Phi) = square(config.initial_sigma_phi);
  cov(QOverPt, QOverPt) = square(config.initial_sigma_qop_t);
  cov(TanLambda, TanLambda) = square(config.initial_sigma_tanl);

  double chi2 = 0.0;
  int ndof = 0;
  const StateMatrix eye = StateMatrix::Identity();
  RknDiagnostics rkn_diagnostics;

  result.measurement_chi2.reserve(npoints);
  result.measurement_used.reserve(npoints);
  if (config.collect_innovation_components)
  {
    result.measurement_in_seed.reserve(npoints);
    result.innovation_residual_r.reserve(npoints);
    result.innovation_residual_rphi.reserve(npoints);
    result.innovation_residual_z.reserve(npoints);
    result.prediction_sigma_r.reserve(npoints);
    result.prediction_sigma_rphi.reserve(npoints);
    result.prediction_sigma_z.reserve(npoints);
    result.innovation_sigma_r.reserve(npoints);
    result.innovation_sigma_rphi.reserve(npoints);
    result.innovation_sigma_z.reserve(npoints);
    result.innovation_rho_r_rphi.reserve(npoints);
    result.innovation_rho_r_z.reserve(npoints);
    result.innovation_rho_rphi_z.reserve(npoints);
    result.innovation_whitened_0.reserve(npoints);
    result.innovation_whitened_1.reserve(npoints);
    result.innovation_whitened_2.reserve(npoints);
  }

  const auto copy_rkn_diagnostics = [&]()
  {
    result.rkn_propagations = rkn_diagnostics.propagations;
    result.rkn_accepted_steps = rkn_diagnostics.accepted_steps;
    result.rkn_rejected_trials = rkn_diagnostics.rejected_trials;
    result.rkn_max_trial_accepts = rkn_diagnostics.max_trial_accepts;
    result.rkn_failures = rkn_diagnostics.failures;
    result.rkn_seconds = rkn_diagnostics.seconds;
  };

  for (std::size_t index = 0; index < npoints; ++index)
  {
    StateVector pred_state = state;
    StateMatrix pred_cov = cov;
    StateMatrix fmat = eye;
    if (index > 0)
    {
      const double ds = result.path_s[index] - result.path_s[index - 1];
      fmat = transport_jacobian(state, ds, config, mass_gev, &rkn_diagnostics);
      pred_state = (config.magnetic_field == nullptr && config.analytic_uniform_propagation)
                       ? propagate_uniform_bz(state, ds, config.bfield_t, config, mass_gev)
                       : propagate_rkn4(state, ds, config, mass_gev, &rkn_diagnostics);
      if (!fmat.allFinite() || !pred_state.allFinite())
      {
        copy_rkn_diagnostics();
        result.message = "RK propagation exceeded its step budget or became non-finite";
        return false;
      }
      pred_cov = fmat * cov * fmat.transpose() + process_noise(state, ds, config, mass_gev);
      pred_cov = 0.5 * (pred_cov + pred_cov.transpose()).eval();
    }

    MeasurementVector measurement;
    measurement << points[index].position.x, points[index].position.y, points[index].position.z;
    const MeasurementCov meas_cov = measurement_covariance(points[index], var_rphi, var_r, var_z);
    const MeasurementVector meas_residual = measurement - hmat * pred_state;
    const MeasurementCov innovation = hmat * pred_cov * hmat.transpose() + meas_cov;
    const MeasurementCov innovation_inv =
        innovation.completeOrthogonalDecomposition().solve(MeasurementCov::Identity());
    const double step_chi2 =
        (meas_residual.transpose() * innovation_inv * meas_residual)(0, 0);

    if (config.collect_innovation_components)
    {
      const MeasurementCov local_rotation = measurement_local_rotation(points[index]);
      const MeasurementVector local_residual = local_rotation * meas_residual;
      const MeasurementCov predicted_position_cov = hmat * pred_cov * hmat.transpose();
      const MeasurementCov local_prediction_cov =
          local_rotation * predicted_position_cov * local_rotation.transpose();
      const MeasurementCov local_innovation =
          local_rotation * innovation * local_rotation.transpose();

      const double prediction_sigma_r = covariance_sigma(local_prediction_cov, 0);
      const double prediction_sigma_rphi = covariance_sigma(local_prediction_cov, 1);
      const double prediction_sigma_z = covariance_sigma(local_prediction_cov, 2);
      const double innovation_sigma_r = covariance_sigma(local_innovation, 0);
      const double innovation_sigma_rphi = covariance_sigma(local_innovation, 1);
      const double innovation_sigma_z = covariance_sigma(local_innovation, 2);
      const MeasurementVector whitened = whiten_innovation(local_innovation, local_residual);

      // The current seed is a global helix fit, so every measurement contributes
      // to it. This marker will distinguish bootstrap seed points once the seed
      // is restricted to an initial consecutive subset.
      result.measurement_in_seed.push_back(1U);
      result.innovation_residual_r.push_back(local_residual(0));
      result.innovation_residual_rphi.push_back(local_residual(1));
      result.innovation_residual_z.push_back(local_residual(2));
      result.prediction_sigma_r.push_back(prediction_sigma_r);
      result.prediction_sigma_rphi.push_back(prediction_sigma_rphi);
      result.prediction_sigma_z.push_back(prediction_sigma_z);
      result.innovation_sigma_r.push_back(innovation_sigma_r);
      result.innovation_sigma_rphi.push_back(innovation_sigma_rphi);
      result.innovation_sigma_z.push_back(innovation_sigma_z);
      result.innovation_rho_r_rphi.push_back(
          covariance_correlation(local_innovation, 0, 1,
                                 innovation_sigma_r, innovation_sigma_rphi));
      result.innovation_rho_r_z.push_back(
          covariance_correlation(local_innovation, 0, 2,
                                 innovation_sigma_r, innovation_sigma_z));
      result.innovation_rho_rphi_z.push_back(
          covariance_correlation(local_innovation, 1, 2,
                                 innovation_sigma_rphi, innovation_sigma_z));
      result.innovation_whitened_0.push_back(whitened(0));
      result.innovation_whitened_1.push_back(whitened(1));
      result.innovation_whitened_2.push_back(whitened(2));
    }

    const Eigen::Matrix<double, 6, 3> gain = pred_cov * hmat.transpose() * innovation_inv;
    state = pred_state + gain * meas_residual;
    state(Phi) = normalize_phi(state(Phi));
    cov = (eye - gain * hmat) * pred_cov * (eye - gain * hmat).transpose() +
          gain * meas_cov * gain.transpose();
    cov = 0.5 * (cov + cov.transpose()).eval();

    result.measurement_chi2.push_back(step_chi2);
    result.measurement_used.push_back(1U);
    ++result.naccepted;
    chi2 += step_chi2;
    ndof += 3;

    states_predicted[index] = pred_state;
    covs_predicted[index] = pred_cov;
    transport[index] = fmat;
    states_filtered[index] = state;
    covs_filtered[index] = cov;
  }

  std::vector<StateVector> states_smoothed = states_filtered;
  std::vector<StateMatrix> covs_smoothed = covs_filtered;
  for (std::size_t next_index = npoints - 1; next_index > 0; --next_index)
  {
    const std::size_t index = next_index - 1;
    const StateMatrix pred_inv = covs_predicted[next_index]
                                     .completeOrthogonalDecomposition()
                                     .solve(StateMatrix::Identity());
    const StateMatrix smoother_gain =
        covs_filtered[index] *
        transport[next_index].transpose() *
        pred_inv;
    const StateVector smooth_residual = residual(states_smoothed[next_index],
                                                 states_predicted[next_index]);
    states_smoothed[index] = states_filtered[index] + smoother_gain * smooth_residual;
    states_smoothed[index](Phi) = normalize_phi(states_smoothed[index](Phi));
    covs_smoothed[index] =
        covs_filtered[index] +
        smoother_gain *
            (covs_smoothed[next_index] - covs_predicted[next_index]) *
            smoother_gain.transpose();
    covs_smoothed[index] =
        0.5 * (covs_smoothed[index] +
               covs_smoothed[index].transpose())
                  .eval();
  }

  result.states_filtered.reserve(npoints);
  result.covs_filtered.reserve(npoints);
  result.states_smoothed.reserve(npoints);
  result.covs_smoothed.reserve(npoints);
  for (std::size_t index = 0; index < npoints; ++index)
  {
    result.states_filtered.push_back(to_array(states_filtered[index]));
    result.covs_filtered.push_back(to_array(covs_filtered[index]));
    result.states_smoothed.push_back(to_array(states_smoothed[index]));
    result.covs_smoothed.push_back(to_array(covs_smoothed[index]));
  }

  result.chi2 = chi2;
  result.ndof = ndof - StateDim;
  copy_rkn_diagnostics();
  result.success = true;
  result.message = "ok";
  return true;
}

TpcTrackVec3 TpcTrackKalmanFitter::state_position(const std::array<double, StateDim> &state)
{
  return {state[X], state[Y], state[Z]};
}

TpcTrackVec3 TpcTrackKalmanFitter::state_momentum(const std::array<double, StateDim> &state)
{
  const double qop_t = state[QOverPt];
  const double pt = (std::abs(qop_t) < 1.0e-12) ? 1.0e12 : 1.0 / std::abs(qop_t);
  return {
      pt * std::cos(state[Phi]),
      pt * std::sin(state[Phi]),
      pt * state[TanLambda]};
}

TpcTrackVec3 TpcTrackKalmanFitter::state_tangent(const std::array<double, StateDim> &state)
{
  return {
      std::cos(state[Phi]),
      std::sin(state[Phi]),
      state[TanLambda]};
}

std::array<double, TpcTrackKalmanFitter::StateDim> TpcTrackKalmanFitter::propagation_state(
    const TpcKalmanResult &fit,
    const TpcTrackVec3 & /*reference_vertex*/)
{
  if (!fit.success || fit.states_smoothed.empty())
  {
    return {};
  }

  if (fit.charge == 0)
  {
    return fit.states_smoothed.front();
  }

  const double charge_sign = static_cast<double>((fit.charge > 0) ? 1 : -1);
  const double physical_qop_sign = -charge_sign;
  const double sequence_qop = fit.states_smoothed.front()[QOverPt];
  const bool sequence_runs_against_physical = sequence_qop * physical_qop_sign < 0.0;

  // The fitted states are assumed to be in a continuous along-track sequence.
  // Charge fixes only which end of that sequence is the physical initial point;
  // do not use the reference vertex or any transverse-distance heuristic here.
  auto state = sequence_runs_against_physical ? fit.states_smoothed.back()
                                              : fit.states_smoothed.front();

  const double qop_t = state[QOverPt];
  if (std::abs(qop_t) < 1.0e-12)
  {
    return state;
  }

  const double pt = 1.0 / std::abs(qop_t);
  const double fit_omega = kCurvaturePerCm * fit.bfield_t * qop_t;
  if (std::abs(fit_omega) < 1.0e-12)
  {
    return state;
  }

  const double fit_center_x = state[X] - std::sin(state[Phi]) / fit_omega;
  const double fit_center_y = state[Y] + std::cos(state[Phi]) / fit_omega;

  // Internal convention: omega = 0.003 * B * qop_t.  For a physical charge in
  // a solenoidal field this corresponds to qop_t = -charge / pT.  The Kalman
  // fit itself may have the opposite sign if the input point order was reversed,
  // so choose the tangent direction that preserves the fitted circle center.
  const double physical_qop_t = -charge_sign / pt;
  const double physical_omega = kCurvaturePerCm * fit.bfield_t * physical_qop_t;
  if (std::abs(physical_omega) < 1.0e-12)
  {
    return state;
  }

  auto center_distance2 = [&](const double phi)
  {
    const double cx = state[X] - std::sin(phi) / physical_omega;
    const double cy = state[Y] + std::cos(phi) / physical_omega;
    return square(cx - fit_center_x) + square(cy - fit_center_y);
  };

  const double phi_keep = normalize_phi(state[Phi]);
  const double phi_flip = normalize_phi(state[Phi] + kPi);
  if (center_distance2(phi_flip) < center_distance2(phi_keep))
  {
    state[Phi] = phi_flip;
    state[TanLambda] *= -1.0;
  }
  else
  {
    state[Phi] = phi_keep;
  }
  state[QOverPt] = physical_qop_t;
  return state;
}

std::array<double, TpcTrackKalmanFitter::StateDim> TpcTrackKalmanFitter::propagate_state(
    const std::array<double, StateDim> &state,
    const double ds_cm,
    const TpcKalmanConfig &config,
    const double mass_gev)
{
  const StateVector input = to_eigen(state);
  if (config.magnetic_field == nullptr && config.analytic_uniform_propagation)
  {
    return to_array(propagate_uniform_bz(input, ds_cm, config.bfield_t, config, mass_gev));
  }
  return to_array(propagate_rkn4(input, ds_cm, config, mass_gev));
}

std::pair<double, double> TpcTrackKalmanFitter::dca_to_vertex(const TpcKalmanResult &fit,
                                                              const TpcTrackVec3 &vertex,
                                                              const TpcKalmanConfig *input_config)
{
  if (!fit.success || fit.states_smoothed.empty())
  {
    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  }

  const auto state_array = propagation_state(fit, vertex);
  const StateVector state = to_eigen(state_array);
  const double omega = omega_from_state(state, fit.bfield_t);
  if (std::abs(omega) < 1.0e-10)
  {
    return TpcTrackHelixFitter::line_dca_to_vertex(state_position(state_array),
                                                   state_momentum(state_array),
                                                   vertex);
  }

  const double radius = 1.0 / omega;
  const double center_x = state(X) - std::sin(state(Phi)) / omega;
  const double center_y = state(Y) + std::cos(state(Phi)) / omega;
  const double vx = vertex.x - center_x;
  const double vy = vertex.y - center_y;
  const double distance_to_center = std::sqrt(vx * vx + vy * vy);
  double theta_closest = std::atan2(state(Y) - center_y, state(X) - center_x);
  double dca_xy = std::abs(radius);
  if (distance_to_center > 0.0)
  {
    const double closest_x = center_x + std::abs(radius) * vx / distance_to_center;
    const double closest_y = center_y + std::abs(radius) * vy / distance_to_center;
    const double theta0 = std::atan2(state(Y) - center_y, state(X) - center_x);
    const double theta_raw = std::atan2(closest_y - center_y, closest_x - center_x);
    theta_closest = theta0 + normalize_phi(theta_raw - theta0);
    dca_xy = std::abs(distance_to_center - std::abs(radius));
  }

  const double theta0 = std::atan2(state(Y) - center_y, state(X) - center_x);
  const double dtheta = normalize_phi(theta_closest - theta0);
  const double s_cm = dtheta / omega;
  TpcKalmanConfig config = input_config != nullptr ? *input_config : TpcKalmanConfig{};
  config.bfield_t = fit.bfield_t;
  config.magnetic_field = fit.magnetic_field;
  config.analytic_uniform_propagation = fit.analytic_uniform_propagation;
  const auto closest = propagate_state(state_array, s_cm, config, fit.mass_gev);
  return {dca_xy, std::abs(closest[Z] - vertex.z)};
}
