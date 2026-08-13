#include "TpcTrackHelixFitter.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <tuple>

namespace
{
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kCurvatureRefineMinRadiusCm = 150.0;
  constexpr double kCurvatureRefineMaxThetaSpan = 1.0;
  constexpr double kMinTurnDzForBranchingCm = 0.5;
  constexpr double kAnchorResidualSlackCm = 5.0;
  constexpr double kMaximumSearchTurnFraction = 0.95;
  constexpr int kMaxMeasurementTurns = 32;

  template <class T>
  constexpr T square(const T &value)
  {
    return value * value;
  }

  double quiet_nan()
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double unwrap_to_previous(double theta, const double previous)
  {
    while (theta - previous > kPi)
    {
      theta -= 2.0 * kPi;
    }
    while (theta - previous < -kPi)
    {
      theta += 2.0 * kPi;
    }
    return theta;
  }

  bool solve3x3(double a[3][3], double b[3], double x[3])
  {
    for (int i = 0; i < 3; ++i)
    {
      int pivot = i;
      double pivot_abs = std::abs(a[i][i]);
      for (int row = i + 1; row < 3; ++row)
      {
        const double value_abs = std::abs(a[row][i]);
        if (value_abs > pivot_abs)
        {
          pivot = row;
          pivot_abs = value_abs;
        }
      }

      if (pivot_abs < 1e-16)
      {
        return false;
      }

      if (pivot != i)
      {
        for (int col = i; col < 3; ++col)
        {
          std::swap(a[i][col], a[pivot][col]);
        }
        std::swap(b[i], b[pivot]);
      }

      const double diag = a[i][i];
      for (int col = i; col < 3; ++col)
      {
        a[i][col] /= diag;
      }
      b[i] /= diag;

      for (int row = 0; row < 3; ++row)
      {
        if (row == i)
        {
          continue;
        }
        const double factor = a[row][i];
        if (std::abs(factor) < 1e-20)
        {
          continue;
        }
        for (int col = i; col < 3; ++col)
        {
          a[row][col] -= factor * a[i][col];
        }
        b[row] -= factor * b[i];
      }
    }

    x[0] = b[0];
    x[1] = b[1];
    x[2] = b[2];
    return true;
  }

  bool theta_span_for_circle(const std::vector<TpcTrackPoint> &points,
                             const std::size_t nfit,
                             const double cx,
                             const double cy,
                             double &theta_span)
  {
    theta_span = std::numeric_limits<double>::quiet_NaN();
    if (nfit < 2 || points.size() < nfit)
    {
      return false;
    }

    std::vector<double> theta_values;
    theta_values.reserve(nfit);
    for (std::size_t i = 0; i < nfit; ++i)
    {
      const auto &pos = points[i].position;
      if (!std::isfinite(pos.x) || !std::isfinite(pos.y))
      {
        return false;
      }
      double theta = std::atan2(pos.y - cy, pos.x - cx);
      if (!theta_values.empty())
      {
        theta = unwrap_to_previous(theta, theta_values.back());
      }
      theta_values.push_back(theta);
    }

    const auto minmax = std::minmax_element(theta_values.begin(), theta_values.end());
    if (minmax.first == theta_values.end())
    {
      return false;
    }

    theta_span = *minmax.second - *minmax.first;
    return std::isfinite(theta_span) && theta_span >= 0.0;
  }

  bool prefer_curvature_refined_circle(const std::vector<TpcTrackPoint> &points,
                                       const std::size_t nfit,
                                       const double cx,
                                       const double cy,
                                       const double radius)
  {
    if (nfit < 4 || !std::isfinite(radius) || radius < kCurvatureRefineMinRadiusCm)
    {
      return false;
    }

    double theta_span = std::numeric_limits<double>::quiet_NaN();
    if (!theta_span_for_circle(points, nfit, cx, cy, theta_span))
    {
      return false;
    }

    return theta_span < kCurvatureRefineMaxThetaSpan;
  }
}  // namespace

bool TpcTrackHelixFitter::parse_point_order(const std::string &mode, TpcTrackPointOrder &order)
{
  if (mode == "path")
  {
    order = TpcTrackPointOrder::Path;
    return true;
  }
  if (mode == "input")
  {
    order = TpcTrackPointOrder::Input;
    return true;
  }
  if (mode == "radius" || mode == "r")
  {
    order = TpcTrackPointOrder::Radius;
    return true;
  }
  if (mode == "theta-z" || mode == "thetaz" || mode == "theta_z")
  {
    order = TpcTrackPointOrder::ThetaZ;
    return true;
  }
  if (mode == "auto")
  {
    order = TpcTrackPointOrder::Auto;
    return true;
  }
  return false;
}

void TpcTrackHelixFitter::order_points(std::vector<TpcTrackPoint> &points,
                                       const TpcTrackPointOrder order)
{
  if (points.size() <= 2 || order == TpcTrackPointOrder::Input)
  {
    return;
  }

  auto apply_order = [&points](const std::vector<std::size_t> &indices)
  {
    std::vector<TpcTrackPoint> ordered;
    ordered.reserve(indices.size());
    for (const auto index : indices)
    {
      ordered.push_back(points[index]);
    }
    points = std::move(ordered);
  };

  auto path_indices = [&points]()
  {
    std::vector<std::size_t> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::stable_sort(indices.begin(), indices.end(),
                     [&points](const auto lhs, const auto rhs)
                     { return points[lhs].path < points[rhs].path; });
    return indices;
  };

  auto radius_indices = [&points]()
  {
    std::vector<std::size_t> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::stable_sort(indices.begin(), indices.end(),
                     [&points](const auto lhs, const auto rhs)
                     {
                       const double lhs_r = pt(points[lhs].position);
                       const double rhs_r = pt(points[rhs].position);
                       if (std::abs(lhs_r - rhs_r) > 1.0e-6)
                       {
                         return lhs_r < rhs_r;
                       }
                       return points[lhs].layer < points[rhs].layer;
                     });
    return indices;
  };

  auto inner_first = [&points](std::vector<std::size_t> indices)
  {
    if (indices.size() >= 2)
    {
      const double first_r = pt(points[indices.front()].position);
      const double last_r = pt(points[indices.back()].position);
      if (first_r > last_r)
      {
        std::reverse(indices.begin(), indices.end());
      }
    }
    return indices;
  };

  auto theta_z_indices = [&points, &path_indices, &inner_first]()
  {
    double cx = 0.0;
    double cy = 0.0;
    double radius = 0.0;
    if (!fit_circle_least_squares(points, points.size(), cx, cy, radius))
    {
      return path_indices();
    }

    std::vector<double> raw_theta(points.size(), 0.0);
    double z_min = points.front().position.z;
    double z_max = points.front().position.z;
    for (std::size_t i = 0; i < points.size(); ++i)
    {
      raw_theta[i] = std::atan2(points[i].position.y - cy, points[i].position.x - cx);
      z_min = std::min(z_min, points[i].position.z);
      z_max = std::max(z_max, points[i].position.z);
    }

    std::vector<std::size_t> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0);
    if (z_max - z_min > 1.0e-3)
    {
      std::vector<std::size_t> z_order = indices;
      std::stable_sort(z_order.begin(), z_order.end(),
                       [&points](const auto lhs, const auto rhs)
                       { return points[lhs].position.z < points[rhs].position.z; });

      std::vector<double> unwrapped_theta = raw_theta;
      bool have_previous = false;
      double previous = 0.0;
      for (const auto index : z_order)
      {
        double theta = raw_theta[index];
        if (have_previous)
        {
          theta = unwrap_to_previous(theta, previous);
        }
        unwrapped_theta[index] = theta;
        previous = theta;
        have_previous = true;
      }

      std::stable_sort(indices.begin(), indices.end(),
                       [&unwrapped_theta](const auto lhs, const auto rhs)
                       { return unwrapped_theta[lhs] < unwrapped_theta[rhs]; });
      return inner_first(indices);
    }

    std::stable_sort(indices.begin(), indices.end(),
                     [&raw_theta](const auto lhs, const auto rhs)
                     { return raw_theta[lhs] < raw_theta[rhs]; });
    return inner_first(indices);
  };

  auto looks_like_looper = [&points, &path_indices]()
  {
    std::vector<int> layers;
    layers.reserve(points.size());
    for (const auto &point : points)
    {
      layers.push_back(point.layer);
    }
    std::sort(layers.begin(), layers.end());
    if (std::adjacent_find(layers.begin(), layers.end()) != layers.end())
    {
      return true;
    }

    const auto sequence = path_indices();
    std::vector<double> dr_values;
    dr_values.reserve(sequence.size());
    for (std::size_t i = 1; i < sequence.size(); ++i)
    {
      const double dr = pt(points[sequence[i]].position) - pt(points[sequence[i - 1]].position);
      if (std::abs(dr) > 0.5)
      {
        dr_values.push_back(dr);
      }
    }
    for (std::size_t i = 1; i < dr_values.size(); ++i)
    {
      if (dr_values[i] * dr_values[i - 1] < 0.0)
      {
        return true;
      }
    }

    double cx = 0.0;
    double cy = 0.0;
    double radius = 0.0;
    if (!fit_circle_least_squares(points, points.size(), cx, cy, radius) ||
        points.size() < 4)
    {
      return false;
    }

    double z_min = points.front().position.z;
    double z_max = points.front().position.z;
    std::vector<std::size_t> z_order(points.size());
    std::iota(z_order.begin(), z_order.end(), 0);
    for (const auto &point : points)
    {
      z_min = std::min(z_min, point.position.z);
      z_max = std::max(z_max, point.position.z);
    }
    if (z_max - z_min <= 1.0e-3)
    {
      return false;
    }

    std::stable_sort(z_order.begin(), z_order.end(),
                     [&points](const auto lhs, const auto rhs)
                     { return points[lhs].position.z < points[rhs].position.z; });

    std::vector<double> theta_values;
    theta_values.reserve(points.size());
    for (const auto index : z_order)
    {
      double theta = std::atan2(points[index].position.y - cy, points[index].position.x - cx);
      if (!theta_values.empty())
      {
        theta = unwrap_to_previous(theta, theta_values.back());
      }
      theta_values.push_back(theta);
    }

    const auto minmax_theta = std::minmax_element(theta_values.begin(), theta_values.end());
    return minmax_theta.first != theta_values.end() &&
           *minmax_theta.second - *minmax_theta.first > 1.5 * kPi;
  };

  switch (order)
  {
  case TpcTrackPointOrder::Path:
    apply_order(path_indices());
    break;
  case TpcTrackPointOrder::Radius:
    apply_order(radius_indices());
    break;
  case TpcTrackPointOrder::ThetaZ:
    apply_order(theta_z_indices());
    break;
  case TpcTrackPointOrder::Auto:
    apply_order(looks_like_looper() ? theta_z_indices() : radius_indices());
    break;
  case TpcTrackPointOrder::Input:
    break;
  }
}

bool TpcTrackHelixFitter::fit_circle_least_squares(const std::vector<TpcTrackPoint> &points,
                                                   const std::size_t nfit,
                                                   double &cx,
                                                   double &cy,
                                                   double &radius)
{
  if (nfit < 3 || points.size() < nfit)
  {
    return false;
  }

  Eigen::MatrixXd matrix(nfit, 3);
  Eigen::VectorXd rhs(nfit);
  for (std::size_t i = 0; i < nfit; ++i)
  {
    const auto &pos = points[i].position;
    matrix(static_cast<int>(i), 0) = pos.x;
    matrix(static_cast<int>(i), 1) = pos.y;
    matrix(static_cast<int>(i), 2) = 1.0;
    rhs(static_cast<int>(i)) = -(square(pos.x) + square(pos.y));
  }

  const Eigen::Vector3d solution = matrix.colPivHouseholderQr().solve(rhs);
  cx = -0.5 * solution(0);
  cy = -0.5 * solution(1);
  const double radius2 = square(cx) + square(cy) - solution(2);
  if (radius2 <= 0.0 || !std::isfinite(radius2))
  {
    return false;
  }

  radius = std::sqrt(radius2);
  return radius > 0.0 && std::isfinite(radius);
}

bool TpcTrackHelixFitter::fit_circle_curvature_refined(const std::vector<TpcTrackPoint> &points,
                                                       const std::size_t nfit,
                                                       double &cx,
                                                       double &cy,
                                                       double &radius)
{
  cx = quiet_nan();
  cy = quiet_nan();
  radius = quiet_nan();

  if (nfit < 4 || points.size() < nfit)
  {
    return false;
  }

  std::vector<TpcTrackVec3> positions;
  positions.reserve(nfit);
  for (std::size_t i = 0; i < nfit; ++i)
  {
    const auto &pos = points[i].position;
    if (!finite(pos))
    {
      continue;
    }
    positions.push_back(pos);
  }

  if (positions.size() < 4)
  {
    return false;
  }

  const double npoints = static_cast<double>(positions.size());

  double mean_x = 0.0;
  double mean_y = 0.0;
  for (const auto &pos : positions)
  {
    mean_x += pos.x;
    mean_y += pos.y;
  }
  mean_x /= npoints;
  mean_y /= npoints;

  double cuu = 0.0;
  double cvv = 0.0;
  double cuv = 0.0;
  for (const auto &pos : positions)
  {
    const double u = pos.x - mean_x;
    const double v = pos.y - mean_y;
    cuu += u * u;
    cvv += v * v;
    cuv += u * v;
  }
  const double phi = 0.5 * std::atan2(2.0 * cuv, cuu - cvv);
  const double cphi = std::cos(phi);
  const double sphi = std::sin(phi);

  double sum_s = 0.0;
  double sum_s2 = 0.0;
  double sum_s3 = 0.0;
  double sum_s4 = 0.0;
  double sum_n = 0.0;
  double sum_sn = 0.0;
  double sum_s2n = 0.0;
  for (const auto &pos : positions)
  {
    const double u = pos.x - mean_x;
    const double v = pos.y - mean_y;
    const double s = cphi * u + sphi * v;
    const double n = -sphi * u + cphi * v;
    sum_s += s;
    sum_s2 += s * s;
    sum_s3 += s * s * s;
    sum_s4 += s * s * s * s;
    sum_n += n;
    sum_sn += s * n;
    sum_s2n += s * s * n;
  }

  double seed_matrix[3][3] = {
      {npoints, sum_s, sum_s2},
      {sum_s, sum_s2, sum_s3},
      {sum_s2, sum_s3, sum_s4}};
  double seed_rhs[3] = {sum_n, sum_sn, sum_s2n};
  double seed_coeffs[3] = {0.0, 0.0, 0.0};
  if (!solve3x3(seed_matrix, seed_rhs, seed_coeffs))
  {
    return false;
  }

  const double a = seed_coeffs[0];
  const double b = seed_coeffs[1];
  const double c = seed_coeffs[2];
  const double denom = std::pow(1.0 + b * b, 1.5);
  if (!(denom > 0.0) || !std::isfinite(denom))
  {
    return false;
  }

  const double kappa = 2.0 * c / denom;
  if (!std::isfinite(kappa) || !(std::abs(kappa) > 1e-14))
  {
    return false;
  }

  auto circle_from_local = [&](const double a_par,
                               const double b_par,
                               const double kappa_par,
                               double &cx_out,
                               double &cy_out,
                               double &radius_out) -> bool
  {
    if (!std::isfinite(a_par) || !std::isfinite(b_par) || !std::isfinite(kappa_par))
    {
      return false;
    }
    if (!(std::abs(kappa_par) > 1e-18))
    {
      return false;
    }

    const double slope_norm = std::sqrt(1.0 + b_par * b_par);
    if (!(slope_norm > 0.0) || !std::isfinite(slope_norm))
    {
      return false;
    }

    const double center_s = -b_par / (slope_norm * kappa_par);
    const double center_n = a_par + 1.0 / (slope_norm * kappa_par);
    const double radius_local = std::abs(1.0 / kappa_par);
    if (!(radius_local > 0.0) || !std::isfinite(radius_local))
    {
      return false;
    }

    cx_out = mean_x + center_s * cphi - center_n * sphi;
    cy_out = mean_y + center_s * sphi + center_n * cphi;
    radius_out = radius_local;
    return std::isfinite(cx_out) && std::isfinite(cy_out);
  };

  auto residuals_from_local = [&](const double a_par,
                                  const double b_par,
                                  const double kappa_par,
                                  std::vector<double> &residuals,
                                  double &chi2,
                                  double *cx_out = nullptr,
                                  double *cy_out = nullptr,
                                  double *radius_out = nullptr) -> bool
  {
    double cx_fit = 0.0;
    double cy_fit = 0.0;
    double radius_fit = 0.0;
    if (!circle_from_local(a_par, b_par, kappa_par, cx_fit, cy_fit, radius_fit))
    {
      return false;
    }

    residuals.resize(positions.size());
    chi2 = 0.0;
    for (std::size_t i = 0; i < positions.size(); ++i)
    {
      const double dx = positions[i].x - cx_fit;
      const double dy = positions[i].y - cy_fit;
      const double dist = std::sqrt(dx * dx + dy * dy);
      if (!std::isfinite(dist))
      {
        return false;
      }
      const double resid = dist - radius_fit;
      residuals[i] = resid;
      chi2 += resid * resid;
    }

    if (cx_out)
    {
      *cx_out = cx_fit;
    }
    if (cy_out)
    {
      *cy_out = cy_fit;
    }
    if (radius_out)
    {
      *radius_out = radius_fit;
    }
    return std::isfinite(chi2);
  };

  double a_cur = a;
  double b_cur = b;
  double kappa_cur = kappa;
  double cx_cur = 0.0;
  double cy_cur = 0.0;
  double radius_cur = 0.0;
  std::vector<double> residuals_cur;
  double chi2_cur = 0.0;
  if (!residuals_from_local(a_cur, b_cur, kappa_cur, residuals_cur, chi2_cur,
                            &cx_cur, &cy_cur, &radius_cur))
  {
    return false;
  }

  double lambda = 1e-3;
  for (int iter = 0; iter < 30; ++iter)
  {
    const double step_a = std::max(1e-8, 1e-6 * std::max({1.0, std::abs(a_cur), 0.01 * radius_cur}));
    const double step_b = std::max(1e-8, 1e-6 * std::max(1.0, std::abs(b_cur)));
    const double step_kappa = std::max(1e-12, 1e-6 * std::max(std::abs(kappa_cur), 1.0 / std::max(radius_cur, 1.0)));
    const double steps[3] = {step_a, step_b, step_kappa};

    std::vector<double> jac_cols[3];
    bool jacobian_ok = true;
    for (int ipar = 0; ipar < 3; ++ipar)
    {
      jac_cols[ipar].assign(positions.size(), 0.0);

      double a_plus = a_cur;
      double b_plus = b_cur;
      double kappa_plus = kappa_cur;
      double a_minus = a_cur;
      double b_minus = b_cur;
      double kappa_minus = kappa_cur;
      if (ipar == 0)
      {
        a_plus += steps[ipar];
        a_minus -= steps[ipar];
      }
      else if (ipar == 1)
      {
        b_plus += steps[ipar];
        b_minus -= steps[ipar];
      }
      else
      {
        kappa_plus += steps[ipar];
        kappa_minus -= steps[ipar];
      }

      std::vector<double> residuals_plus;
      std::vector<double> residuals_minus;
      double chi2_dummy = 0.0;
      const bool have_plus = residuals_from_local(a_plus, b_plus, kappa_plus, residuals_plus, chi2_dummy);
      const bool have_minus = residuals_from_local(a_minus, b_minus, kappa_minus, residuals_minus, chi2_dummy);
      if (!have_plus && !have_minus)
      {
        jacobian_ok = false;
        break;
      }

      for (std::size_t ipoint = 0; ipoint < positions.size(); ++ipoint)
      {
        if (have_plus && have_minus)
        {
          jac_cols[ipar][ipoint] = (residuals_plus[ipoint] - residuals_minus[ipoint]) / (2.0 * steps[ipar]);
        }
        else if (have_plus)
        {
          jac_cols[ipar][ipoint] = (residuals_plus[ipoint] - residuals_cur[ipoint]) / steps[ipar];
        }
        else
        {
          jac_cols[ipar][ipoint] = (residuals_cur[ipoint] - residuals_minus[ipoint]) / steps[ipar];
        }
      }
    }
    if (!jacobian_ok)
    {
      break;
    }

    double jtj[3][3] = {};
    double jtr[3] = {};
    for (std::size_t ipoint = 0; ipoint < positions.size(); ++ipoint)
    {
      for (int i = 0; i < 3; ++i)
      {
        const double ji = jac_cols[i][ipoint];
        jtr[i] += ji * residuals_cur[ipoint];
        for (int j = 0; j < 3; ++j)
        {
          jtj[i][j] += ji * jac_cols[j][ipoint];
        }
      }
    }

    bool accepted = false;
    double delta[3] = {0.0, 0.0, 0.0};
    for (int trial = 0; trial < 8; ++trial)
    {
      double system[3][3] = {
          {jtj[0][0], jtj[0][1], jtj[0][2]},
          {jtj[1][0], jtj[1][1], jtj[1][2]},
          {jtj[2][0], jtj[2][1], jtj[2][2]}};
      for (int idiag = 0; idiag < 3; ++idiag)
      {
        system[idiag][idiag] += lambda * std::max(jtj[idiag][idiag], 1.0);
      }

      double neg_jtr[3] = {-jtr[0], -jtr[1], -jtr[2]};
      if (!solve3x3(system, neg_jtr, delta))
      {
        lambda *= 10.0;
        continue;
      }

      const double a_try = a_cur + delta[0];
      const double b_try = b_cur + delta[1];
      const double kappa_try = kappa_cur + delta[2];
      std::vector<double> residuals_try;
      double chi2_try = 0.0;
      double cx_try = 0.0;
      double cy_try = 0.0;
      double radius_try = 0.0;
      if (!residuals_from_local(a_try, b_try, kappa_try, residuals_try, chi2_try,
                                &cx_try, &cy_try, &radius_try) ||
          !(chi2_try < chi2_cur))
      {
        lambda *= 10.0;
        continue;
      }

      a_cur = a_try;
      b_cur = b_try;
      kappa_cur = kappa_try;
      residuals_cur.swap(residuals_try);
      chi2_cur = chi2_try;
      cx_cur = cx_try;
      cy_cur = cy_try;
      radius_cur = radius_try;
      lambda = std::max(1e-12, lambda * 0.3);
      accepted = true;
      break;
    }

    if (!accepted)
    {
      break;
    }

    const double rel_step = std::sqrt(
        std::pow(delta[0] / std::max({1.0, std::abs(a_cur), 0.01 * radius_cur}), 2) +
        std::pow(delta[1] / std::max(1.0, std::abs(b_cur)), 2) +
        std::pow(delta[2] / std::max(std::abs(kappa_cur), 1e-12), 2));
    if (rel_step < 1e-10)
    {
      break;
    }
  }

  cx = cx_cur;
  cy = cy_cur;
  radius = radius_cur;
  return std::isfinite(cx) && std::isfinite(cy) &&
         std::isfinite(radius) && radius > 0.0;
}

bool TpcTrackHelixFitter::fit(const std::vector<TpcTrackPoint> &points,
                              const int fit_first_points,
                              const double bfield_t,
                              TpcTrackHelix &helix)
{
  const std::size_t nfit =
      (fit_first_points > 0) ? std::min(points.size(), static_cast<std::size_t>(fit_first_points)) : points.size();
  if (nfit < 3)
  {
    return false;
  }

  double cx = 0.0;
  double cy = 0.0;
  double radius = 0.0;
  bool have_circle = fit_circle_least_squares(points, nfit, cx, cy, radius);
  if (have_circle && prefer_curvature_refined_circle(points, nfit, cx, cy, radius))
  {
    double refined_cx = 0.0;
    double refined_cy = 0.0;
    double refined_radius = 0.0;
    if (fit_circle_curvature_refined(points, nfit, refined_cx, refined_cy, refined_radius))
    {
      cx = refined_cx;
      cy = refined_cy;
      radius = refined_radius;
    }
  }
  else if (!have_circle)
  {
    have_circle = fit_circle_curvature_refined(points, nfit, cx, cy, radius);
  }

  if (!have_circle)
  {
    return false;
  }

  Eigen::MatrixXd z_matrix(nfit, 2);
  Eigen::VectorXd z_rhs(nfit);
  std::vector<double> theta_values;
  theta_values.reserve(nfit);
  for (std::size_t i = 0; i < nfit; ++i)
  {
    double theta = std::atan2(points[i].position.y - cy, points[i].position.x - cx);
    if (!theta_values.empty())
    {
      theta = unwrap_to_previous(theta, theta_values.back());
    }
    theta_values.push_back(theta);
    z_matrix(static_cast<int>(i), 0) = theta;
    z_matrix(static_cast<int>(i), 1) = 1.0;
    z_rhs(static_cast<int>(i)) = points[i].position.z;
  }

  const auto minmax_theta = std::minmax_element(theta_values.begin(), theta_values.end());
  if (minmax_theta.first == theta_values.end() ||
      *minmax_theta.second - *minmax_theta.first < 1e-4)
  {
    return false;
  }

  const Eigen::Vector2d z_solution = z_matrix.colPivHouseholderQr().solve(z_rhs);
  double direction = (theta_values.back() > theta_values.front()) ? 1.0 : -1.0;
  if (std::abs(theta_values.back() - theta_values.front()) <= 0.0)
  {
    direction = 1.0;
  }

  helix.cx = cx;
  helix.cy = cy;
  helix.radius = radius;
  helix.z0 = z_solution(1);
  helix.pitch = z_solution(0);
  helix.theta_first = theta_values.front();
  helix.theta_last = theta_values.back();
  helix.theta_min = *minmax_theta.first;
  helix.theta_max = *minmax_theta.second;
  helix.direction = direction;
  helix.bfield_t = bfield_t;
  return true;
}

bool TpcTrackHelixFitter::from_state(const TpcTrackVec3 &position,
                                     const TpcTrackVec3 &momentum_value,
                                     const int charge,
                                     const double bfield_t,
                                     TpcTrackHelix &helix)
{
  if (charge == 0 || bfield_t == 0.0 || !finite(position) || !finite(momentum_value))
  {
    return false;
  }

  const double p_t = pt(momentum_value);
  if (p_t <= 0.0 || !std::isfinite(p_t))
  {
    return false;
  }

  const double radius = p_t / (0.003 * std::abs(bfield_t));
  if (radius <= 0.0 || !std::isfinite(radius))
  {
    return false;
  }

  const double direction = -static_cast<double>(charge) * ((bfield_t > 0.0) ? 1.0 : -1.0);
  const double radial_x = direction * momentum_value.y / p_t;
  const double radial_y = -direction * momentum_value.x / p_t;
  const double cx = position.x - radius * radial_x;
  const double cy = position.y - radius * radial_y;
  const double theta = std::atan2(position.y - cy, position.x - cx);
  const double pitch = momentum_value.z * radius / (direction * p_t);
  const double z0 = position.z - pitch * theta;

  helix.cx = cx;
  helix.cy = cy;
  helix.radius = radius;
  helix.z0 = z0;
  helix.pitch = pitch;
  helix.theta_first = theta;
  helix.theta_last = theta;
  helix.theta_min = theta;
  helix.theta_max = theta;
  helix.direction = direction;
  helix.bfield_t = bfield_t;
  return finite(point(helix, theta)) && finite(momentum(helix, theta));
}

bool TpcTrackHelixFitter::orient_to_charge(TpcTrackHelix &helix, const int charge)
{
  if (charge == 0 || helix.bfield_t == 0.0 || !std::isfinite(helix.bfield_t))
  {
    return false;
  }

  const double charge_sign = static_cast<double>((charge > 0) ? 1 : -1);
  const double bfield_sign = (helix.bfield_t > 0.0) ? 1.0 : -1.0;
  const double physical_direction = -charge_sign * bfield_sign;

  // The geometric fit may be ordered opposite to physical time.  Keep the
  // fitted helix, but make theta_first the physical initial branch.
  if (helix.direction * physical_direction < 0.0)
  {
    std::swap(helix.theta_first, helix.theta_last);
  }
  helix.direction = physical_direction;
  return true;
}

TpcTrackVec3 TpcTrackHelixFitter::point(const TpcTrackHelix &helix, const double theta)
{
  return {
      helix.cx + helix.radius * std::cos(theta),
      helix.cy + helix.radius * std::sin(theta),
      helix.z0 + helix.pitch * theta};
}

TpcTrackVec3 TpcTrackHelixFitter::tangent(const TpcTrackHelix &helix, const double theta)
{
  return {
      -helix.radius * std::sin(theta),
      helix.radius * std::cos(theta),
      helix.pitch};
}

TpcTrackVec3 TpcTrackHelixFitter::momentum(const TpcTrackHelix &helix, const double theta)
{
  double p_t = 0.3 * std::abs(helix.bfield_t) * (helix.radius / 100.0);
  if (p_t <= 0.0 || !std::isfinite(p_t))
  {
    p_t = 1.0;
  }

  return {
      helix.direction * p_t * (-std::sin(theta)),
      helix.direction * p_t * std::cos(theta),
      helix.direction * p_t * helix.pitch / helix.radius};
}

std::pair<double, double> TpcTrackHelixFitter::theta_search_range(const TpcTrackHelix &helix,
                                                                  const double theta_extension,
                                                                  const double downstream_margin)
{
  const double upstream = helix.theta_first - helix.direction * theta_extension;
  const double downstream = helix.theta_first + helix.direction * downstream_margin;
  return {std::min(upstream, downstream), std::max(upstream, downstream)};
}

bool TpcTrackHelixFitter::measurement_anchored_search_range(
    const TpcTrackHelix &helix,
    const std::vector<TpcTrackPoint> &points,
    const double max_upstream_cm,
    const double downstream_margin_cm,
    TpcTrackHelixSearchRange &range)
{
  range = {};
  if (points.empty() || !(helix.radius > 0.0) ||
      !std::isfinite(helix.radius) || !std::isfinite(helix.theta_first) ||
      !std::isfinite(helix.direction) || std::abs(helix.direction) < 0.5)
  {
    return false;
  }

  const double direction = helix.direction > 0.0 ? 1.0 : -1.0;
  const double reference_theta = helix.theta_first;
  const double circumference_cm = 2.0 * kPi * helix.radius;
  if (!(circumference_cm > 0.0) || !std::isfinite(circumference_cm))
  {
    return false;
  }

  struct Projection
  {
    int point_index{-1};
    double theta{0.0};
    double path_cm{0.0};
    double residual_cm{0.0};
  };

  std::vector<Projection> projections;
  projections.reserve(points.size());
  double minimum_residual_cm = std::numeric_limits<double>::infinity();
  const double turn_dz_cm = 2.0 * kPi * helix.pitch * direction;

  for (std::size_t index = 0; index < points.size(); ++index)
  {
    const auto &measurement = points[index].position;
    if (!finite(measurement))
    {
      continue;
    }

    const double raw_theta = std::atan2(measurement.y - helix.cy,
                                        measurement.x - helix.cx);
    // Path is signed relative to the transverse perigee. Measurements may be
    // either before or after that point in physical time, so start from the
    // nearest angular branch and let z resolve any additional full turns.
    const double base_phase = std::remainder(
        direction * (raw_theta - reference_theta), 2.0 * kPi);
    const double base_theta = reference_theta + direction * base_phase;
    const TpcTrackVec3 base_position = point(helix, base_theta);

    int turn_guess = 0;
    if (std::abs(turn_dz_cm) >= kMinTurnDzForBranchingCm)
    {
      turn_guess = static_cast<int>(std::llround(
          (measurement.z - base_position.z) / turn_dz_cm));
      turn_guess = std::clamp(turn_guess,
                              -kMaxMeasurementTurns,
                              kMaxMeasurementTurns);
    }

    const int first_turn = std::max(-kMaxMeasurementTurns, turn_guess - 2);
    const int last_turn = std::min(kMaxMeasurementTurns, turn_guess + 2);
    Projection best;
    best.point_index = static_cast<int>(index);
    best.residual_cm = std::numeric_limits<double>::infinity();
    int best_turn = -1;
    for (int turn = first_turn; turn <= last_turn; ++turn)
    {
      const double theta = base_theta + direction * 2.0 * kPi * turn;
      const TpcTrackVec3 projected = point(helix, theta);
      const double residual_cm = distance(projected, measurement);
      if (!std::isfinite(residual_cm))
      {
        continue;
      }

      if (residual_cm < best.residual_cm - 1.0e-12 ||
          (std::abs(residual_cm - best.residual_cm) <= 1.0e-12 &&
           (best_turn < 0 || turn < best_turn)))
      {
        best.theta = theta;
        best.path_cm = helix.radius * (base_phase + 2.0 * kPi * turn);
        best.residual_cm = residual_cm;
        best_turn = turn;
      }
    }

    if (!std::isfinite(best.residual_cm))
    {
      continue;
    }
    projections.push_back(best);
    minimum_residual_cm = std::min(minimum_residual_cm, best.residual_cm);
  }

  if (projections.empty() || !std::isfinite(minimum_residual_cm))
  {
    return false;
  }

  const double residual_limit_cm = minimum_residual_cm + kAnchorResidualSlackCm;
  auto anchor_iter = projections.end();
  for (auto iter = projections.begin(); iter != projections.end(); ++iter)
  {
    if (iter->residual_cm > residual_limit_cm)
    {
      continue;
    }
    if (anchor_iter == projections.end() ||
        iter->path_cm < anchor_iter->path_cm - 1.0e-9 ||
        (std::abs(iter->path_cm - anchor_iter->path_cm) <= 1.0e-9 &&
         iter->residual_cm < anchor_iter->residual_cm))
    {
      anchor_iter = iter;
    }
  }
  if (anchor_iter == projections.end())
  {
    return false;
  }

  const double requested_upstream_cm = std::max(0.0, max_upstream_cm);
  const double requested_downstream_cm = std::max(0.0, downstream_margin_cm);
  // Leave a finite gap between the two ends of the search domain. A cap that
  // differs from one turn only by machine epsilon still permits duplicate
  // geometric intersections and rounds back to 2*pi in float QA branches.
  const double maximum_span_cm =
      kMaximumSearchTurnFraction * circumference_cm;
  const double downstream_cm =
      std::min(requested_downstream_cm, maximum_span_cm);
  const double upstream_cm = std::min(
      requested_upstream_cm, std::max(0.0, maximum_span_cm - downstream_cm));
  if (!(upstream_cm + downstream_cm > 0.0))
  {
    return false;
  }

  const double lower_path_cm = anchor_iter->path_cm - upstream_cm;
  const double upper_path_cm = anchor_iter->path_cm + downstream_cm;
  const double theta_a = reference_theta + direction * lower_path_cm / helix.radius;
  const double theta_b = reference_theta + direction * upper_path_cm / helix.radius;

  range.valid = true;
  range.anchor_point_index = anchor_iter->point_index;
  range.anchor_theta = anchor_iter->theta;
  range.anchor_path_cm = anchor_iter->path_cm;
  range.anchor_residual_cm = anchor_iter->residual_cm;
  range.theta_min = std::min(theta_a, theta_b);
  range.theta_max = std::max(theta_a, theta_b);
  range.upstream_cm = upstream_cm;
  range.downstream_cm = downstream_cm;
  return std::isfinite(range.theta_min) && std::isfinite(range.theta_max) &&
         range.theta_max > range.theta_min;
}

bool TpcTrackHelixFitter::line_line_pca(const TpcTrackVec3 &pos1,
                                        const TpcTrackVec3 &dir1,
                                        const TpcTrackVec3 &pos2,
                                        const TpcTrackVec3 &dir2,
                                        TpcTrackLinePca &pca,
                                        const bool normalize_dirs)
{
  TpcTrackVec3 u1 = normalize_dirs ? unit(dir1) : dir1;
  TpcTrackVec3 u2 = normalize_dirs ? unit(dir2) : dir2;
  if (!finite(u1) || !finite(u2))
  {
    return false;
  }

  const TpcTrackVec3 w0 = subtract(pos1, pos2);
  const double a = dot(u1, u1);
  const double b = dot(u1, u2);
  const double c = dot(u2, u2);
  const double d = dot(u1, w0);
  const double e = dot(u2, w0);
  const double denom = a * c - b * b;
  if (std::abs(denom) < 1e-12)
  {
    return false;
  }

  const double s = (b * e - c * d) / denom;
  const double t = (a * e - b * d) / denom;
  pca.pca1 = add(pos1, scale(u1, s));
  pca.pca2 = add(pos2, scale(u2, t));
  pca.dca = distance(pca.pca1, pca.pca2);
  pca.step1 = s;
  pca.step2 = t;
  return true;
}

TpcTrackHelixPca TpcTrackHelixFitter::refine_pair(const TpcTrackHelix &helix1,
                                                  const TpcTrackHelix &helix2,
                                                  double theta1,
                                                  double theta2,
                                                  const double min1,
                                                  const double max1,
                                                  const double min2,
                                                  const double max2,
                                                  double max_step)
{
  double best_dca2 = square(distance(point(helix1, theta1), point(helix2, theta2)));

  for (int iter = 0; iter < 30; ++iter)
  {
    TpcTrackLinePca line_pca;
    if (!line_line_pca(point(helix1, theta1), tangent(helix1, theta1),
                       point(helix2, theta2), tangent(helix2, theta2),
                       line_pca, false))
    {
      break;
    }

    const double step1 = std::clamp(line_pca.step1, -max_step, max_step);
    const double step2 = std::clamp(line_pca.step2, -max_step, max_step);
    if (std::abs(step1) < 1e-5 && std::abs(step2) < 1e-5)
    {
      break;
    }

    const double candidate_theta1 = std::clamp(theta1 + step1, min1, max1);
    const double candidate_theta2 = std::clamp(theta2 + step2, min2, max2);
    const double candidate_dca2 = square(distance(point(helix1, candidate_theta1),
                                                  point(helix2, candidate_theta2)));
    if (candidate_dca2 < best_dca2)
    {
      theta1 = candidate_theta1;
      theta2 = candidate_theta2;
      best_dca2 = candidate_dca2;
    }
    else
    {
      max_step *= 0.5;
      if (max_step < 1e-4)
      {
        break;
      }
    }
  }

  TpcTrackHelixPca output;
  output.theta1 = theta1;
  output.theta2 = theta2;
  output.pca1 = point(helix1, theta1);
  output.pca2 = point(helix2, theta2);
  output.dca = distance(output.pca1, output.pca2);
  return output;
}

std::vector<TpcTrackHelixPca> TpcTrackHelixFitter::pca_candidates(
    const TpcTrackHelix &helix1,
    const TpcTrackHelix &helix2,
    const double theta_extension,
    const int coarse_steps,
    const double downstream_margin,
    const int max_candidates)
{
  TpcTrackHelixSearchRange range1;
  const auto legacy_range1 = theta_search_range(helix1, theta_extension, downstream_margin);
  range1.valid = true;
  range1.theta_min = legacy_range1.first;
  range1.theta_max = legacy_range1.second;
  TpcTrackHelixSearchRange range2;
  const auto legacy_range2 = theta_search_range(helix2, theta_extension, downstream_margin);
  range2.valid = true;
  range2.theta_min = legacy_range2.first;
  range2.theta_max = legacy_range2.second;
  return pca_candidates_in_ranges(helix1, helix2, range1, range2,
                                  coarse_steps, max_candidates);
}

std::vector<TpcTrackHelixPca> TpcTrackHelixFitter::pca_candidates_in_ranges(
    const TpcTrackHelix &helix1,
    const TpcTrackHelix &helix2,
    const TpcTrackHelixSearchRange &range1,
    const TpcTrackHelixSearchRange &range2,
    const int coarse_steps,
    const int max_candidates)
{
  if (!range1.valid || !range2.valid ||
      !std::isfinite(range1.theta_min) || !std::isfinite(range1.theta_max) ||
      !std::isfinite(range2.theta_min) || !std::isfinite(range2.theta_max) ||
      !(range1.theta_max > range1.theta_min) ||
      !(range2.theta_max > range2.theta_min))
  {
    return {};
  }

  const int n_steps = std::max(8, coarse_steps);

  std::vector<double> theta1_values;
  std::vector<double> theta2_values;
  theta1_values.reserve(n_steps);
  theta2_values.reserve(n_steps);
  for (int i = 0; i < n_steps; ++i)
  {
    const double fraction = (n_steps == 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(n_steps - 1);
    theta1_values.push_back(range1.theta_min + fraction * (range1.theta_max - range1.theta_min));
    theta2_values.push_back(range2.theta_min + fraction * (range2.theta_max - range2.theta_min));
  }

  std::vector<std::tuple<double, int, int>> coarse;
  coarse.reserve(static_cast<std::size_t>(n_steps) * static_cast<std::size_t>(n_steps));
  for (int i = 0; i < n_steps; ++i)
  {
    const TpcTrackVec3 point1 = point(helix1, theta1_values[i]);
    for (int j = 0; j < n_steps; ++j)
    {
      const TpcTrackVec3 point2 = point(helix2, theta2_values[j]);
      coarse.emplace_back(square(distance(point1, point2)), i, j);
    }
  }

  std::sort(coarse.begin(), coarse.end(),
            [](const auto &lhs, const auto &rhs)
            { return std::get<0>(lhs) < std::get<0>(rhs); });

  const double max_step = std::max(range1.theta_max - range1.theta_min,
                                   range2.theta_max - range2.theta_min) /
                          static_cast<double>(n_steps);
  const int n_candidates = std::min({std::max(1, max_candidates), static_cast<int>(coarse.size())});

  std::vector<TpcTrackHelixPca> candidates;
  candidates.reserve(n_candidates);
  for (int index = 0; index < n_candidates; ++index)
  {
    const int i = std::get<1>(coarse[index]);
    const int j = std::get<2>(coarse[index]);
    candidates.push_back(refine_pair(
        helix1, helix2, theta1_values[i], theta2_values[j],
        range1.theta_min, range1.theta_max,
        range2.theta_min, range2.theta_max, max_step));
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const TpcTrackHelixPca &lhs, const TpcTrackHelixPca &rhs)
            { return lhs.dca < rhs.dca; });
  return candidates;
}

std::pair<double, double> TpcTrackHelixFitter::helix_dca_to_vertex(const TpcTrackHelix &helix,
                                                                   const TpcTrackVec3 &vertex)
{
  const double vx = vertex.x - helix.cx;
  const double vy = vertex.y - helix.cy;
  const double distance_to_center = std::sqrt(square(vx) + square(vy));

  double theta_raw = helix.theta_first;
  double dca_xy = helix.radius;
  if (distance_to_center > 0.0)
  {
    theta_raw = std::atan2(vy, vx);
    dca_xy = std::abs(distance_to_center - helix.radius);
  }

  const double wraps = std::round((helix.theta_first - theta_raw) / (2.0 * kPi));
  const double theta = theta_raw + 2.0 * kPi * wraps;
  const TpcTrackVec3 closest = point(helix, theta);
  return {dca_xy, std::abs(closest.z - vertex.z)};
}

std::pair<double, double> TpcTrackHelixFitter::line_dca_to_vertex(const TpcTrackVec3 &pos,
                                                                  const TpcTrackVec3 &mom,
                                                                  const TpcTrackVec3 &vertex)
{
  const TpcTrackVec3 rel = subtract(pos, vertex);
  const double pt2 = square(mom.x) + square(mom.y);
  if (pt2 <= 0.0)
  {
    return {quiet_nan(), quiet_nan()};
  }
  const double dca_xy = std::abs(rel.x * mom.y - rel.y * mom.x) / std::sqrt(pt2);
  const double sxy = -(rel.x * mom.x + rel.y * mom.y) / pt2;
  const TpcTrackVec3 closest = add(pos, scale(mom, sxy));
  return {dca_xy, std::abs(closest.z - vertex.z)};
}

bool TpcTrackHelixFitter::finite(const TpcTrackVec3 &value)
{
  return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

TpcTrackVec3 TpcTrackHelixFitter::add(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs)
{
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

TpcTrackVec3 TpcTrackHelixFitter::subtract(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs)
{
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

TpcTrackVec3 TpcTrackHelixFitter::scale(const TpcTrackVec3 &value, const double factor)
{
  return {value.x * factor, value.y * factor, value.z * factor};
}

double TpcTrackHelixFitter::dot(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs)
{
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

double TpcTrackHelixFitter::norm(const TpcTrackVec3 &value)
{
  return std::sqrt(dot(value, value));
}

TpcTrackVec3 TpcTrackHelixFitter::unit(const TpcTrackVec3 &value)
{
  const double length = norm(value);
  if (length <= 0.0 || !std::isfinite(length))
  {
    return {quiet_nan(), quiet_nan(), quiet_nan()};
  }
  return scale(value, 1.0 / length);
}

double TpcTrackHelixFitter::pt(const TpcTrackVec3 &value)
{
  return std::sqrt(square(value.x) + square(value.y));
}

double TpcTrackHelixFitter::distance(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs)
{
  return norm(subtract(lhs, rhs));
}

double TpcTrackHelixFitter::vector_cosine(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs)
{
  const double denom = norm(lhs) * norm(rhs);
  if (denom <= 0.0 || !std::isfinite(denom))
  {
    return quiet_nan();
  }
  return dot(lhs, rhs) / denom;
}
