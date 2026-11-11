#ifndef PHGARFIELD_ADAPTIVEMAP_H
#define PHGARFIELD_ADAPTIVEMAP_H

#include "PHGarfieldForwardMap.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <limits>
#include <string>
#include <vector>

// Temporary adaptive-builder state. Final LUTs are written directly in compact
// form. Coordinates are(phi,r,signed z); corner order is(phi_bit<<2)|(r_bit<<1)|z_bit.
class PHGarfieldAdaptiveMap
{
 public:
  static constexpr const char* ProbeScheme = "tensor_3x3x3_19_noncorner_actual_axis_curvature_v1";
  static constexpr const char* SplitPolicy = "largest_actual_midpoint_error_binary_axis_v1";
  static constexpr const char* AllLostPolicy = "direct_query_no_finite_probe_loss_certificate_v1";
  static constexpr const char* LeafEncoding = "0=Uncertified,1=Certified,2=Mixed,3=SampledAllLost,4=NumericallyUnresolved";
  static constexpr const char* FallbackEncoding = "0=None,1=PhysicalMixed,2=Numerical,3=DepthBudget,4=NodeBudget,5=QueryBudget,6=AngularBranch,7=SampledAllLost";
  // Engineering ceiling only. Actual vector growth remains governed by the
  // explicit per-file policy and checked actual point/cell counts.
  static constexpr std::uint64_t HardMaximumEstimatedBytes = 1024ULL * 1024 * 1024 * 1024;
  // Matches PHGarfieldForwardTransport::ForwardDrift's radial_tolerance_cm. This is a numerical
  // cage-boundary convention (1e-5um), not the micrometre interpolation guard.
  static constexpr double EndpointRadialToleranceCm = 1.e-9;

  enum class LeafClass : int
  {
    Uncertified = 0,
    Certified = 1,
    Mixed = 2,
    AllLost = 3,
    NumericallyUnresolved = 4
  };
  enum class FallbackReason : int
  {
    None = 0,
    PhysicalMixed = 1,
    Numerical = 2,
    DepthBudget = 3,
    NodeBudget = 4,
    QueryBudget = 5,
    AngularBranch = 6,
    SampledAllLost = 7
  };
  struct Geometry
  {
    unsigned int phi_nodes{72};
    unsigned int radial_nodes{32};
    unsigned int longitudinal_nodes{28};
    double inner_radius_cm{21.6};
    double outer_radius_cm{76.4};
    double central_membrane_half_width_cm{0.28};
    double readout_abs_z_cm{102.605};
    double nominal_velocity_cm_ns{0.00755};
  };
  struct Policy
  {
    double guard_um{10.0};
    unsigned int max_depth{12};
    unsigned int max_nodes_per_root{127};
    unsigned int max_queries_per_root{2048};
    std::uint64_t max_total_points{5000000};
    std::uint64_t max_total_cells{2000000};
    std::uint64_t max_estimated_bytes{2ULL * 1024 * 1024 * 1024};
  };
  struct Bounds
  {
    std::array<double, 3> lo{};
    std::array<double, 3> hi{};
  };
  struct Point
  {
    std::uint64_t id{0};
    // Phi is always canonical[0,2pi), including a final-cell upper corner.
    std::array<double, 3> coordinate{};
    int reaches{0};
    int forward_status{0};
    std::array<double, 3> value{};  // dr_cm,r0*dphi_cm,absolute time_ns
  };
  struct Cell
  {
    std::uint64_t id{0};
    std::uint64_t root_id{0};
    std::array<double, 3> lo{};
    std::array<double, 3> hi{};
    std::array<std::int64_t, 8> corner_id{{-1, -1, -1, -1, -1, -1, -1, -1}};
    std::array<std::int64_t, 2> child_id{{-1, -1}};
    int split_axis{-1};  // midpoint; equality belongs to upper child
    unsigned int depth{0};
    int leaf_class{static_cast<int>(LeafClass::Uncertified)};
    int fallback_reason{static_cast<int>(FallbackReason::QueryBudget)};
    unsigned int probe_count{0};  // noncorner probes; complete stencil=19
    unsigned int probe_reaching{0};
    int angular_valid{0};
    double max_transverse_error_um{-1.0};
    double max_time_error_um{-1.0};
    std::array<double, 3> directional_error_um{{-1.0, -1.0, -1.0}};
    std::array<unsigned int, 3> directional_reach_changes{{0, 0, 0}};
  };
  struct Root
  {
    std::uint64_t root_id{0};
    std::uint64_t cell_id{0};
  };
  struct Data
  {
    Geometry geometry;
    Policy policy;
    std::vector<Point> points;
    std::vector<Cell> cells;
    std::vector<Root> roots;
    // Runtime is complete iff[first_root,end_root)==[0,RootCount(geometry)).
    std::uint64_t first_root{0};
    std::uint64_t end_root{0};
    double maximum_drift_time_ns{0.0};
  };
  struct BuildSummary
  {
    std::uint64_t actual_sampler_queries{0};
    std::uint64_t reused_points{0};
    std::uint64_t certified_leaves{0};
    std::uint64_t direct_leaves{0};
    std::uint64_t budget_leaves{0};
  };
  using Sampler = std::function<bool(const std::array<double, 3>&, Point&, std::ostream&)>;

  static Geometry FromConfig(const PHGarfieldForwardMap::Config& c);
  static double CanonicalPhi(double phi);
  static std::uint64_t RootCount(const Geometry& g);
  // Caller checks the gas domain and canonicalizes phi before indexing.
  static std::uint64_t RootIndex(const Geometry& g, unsigned int side, double phi, double r, double z);
  static Bounds RootBounds(const Geometry& g, std::uint64_t id);
  static Bounds ChildBounds(const Bounds& parent, int axis, unsigned int child);
  static std::array<double, 3> CornerCoordinate(const Bounds& b, unsigned int corner);
  static bool Validate(const Data& data, std::ostream& error, bool complete = true);
  // Shared by sampling and builder validation. Invalid samples fail before
  // further refinement; diagnostics preserve full-precision coordinates/values.
  static bool ValidatePoint(const Point& p, const Geometry& g, double maximum_drift_time_ns,
                            std::ostream& error);
  static bool HasUsableCorners(const Data& data, const Cell& cell);
  // Physical/finite/branch integrity is independent of interpolation accuracy.
  template <class CornerValue>
  static bool UsableCornerValues(const Geometry& geometry, const Bounds& bounds,
                                 double maximum_time_ns, const CornerValue& corner_value)
  {
    if (!std::isfinite(maximum_time_ns) || maximum_time_ns < 0.)
    {
      return false;
    }
    double minimum_angle = std::numeric_limits<double>::infinity();
    double maximum_angle = -minimum_angle;
    for (unsigned int corner = 0; corner < 8; ++corner)
    {
      const auto* value = corner_value(corner);
      const auto coordinate = CornerCoordinate(bounds, corner);
      const double radius = coordinate[1];
      if (!value || !std::isfinite(radius) || radius <= 0. ||
          !std::all_of(value->begin(), value->end(), [](double v)
                       { return std::isfinite(v); }) ||
          (*value)[2] < 0. || (*value)[2] > maximum_time_ns ||
          radius + (*value)[0] < geometry.inner_radius_cm - EndpointRadialToleranceCm ||
          radius + (*value)[0] > geometry.outer_radius_cm + EndpointRadialToleranceCm ||
          std::abs((*value)[1] / radius) > std::acos(-1.) ||
          (std::abs(coordinate[2]) == geometry.readout_abs_z_cm && *value != std::array<double, 3>{{0., 0., 0.}}))
      {
        return false;
      }
      const double angle = (*value)[1] / radius;
      minimum_angle = std::min(minimum_angle, angle);
      maximum_angle = std::max(maximum_angle, angle);
    }
    return maximum_angle - minimum_angle <= std::acos(-1.);
  }
  // Both runtime layouts use the same arithmetic and corner order. The
  // accessor returns null for an invalid corner reference.
  template <class CornerValue>
  static bool InterpolateCorners(const Bounds& bounds, const std::array<double, 3>& coordinate,
                                 const CornerValue& corner_value, PHGarfieldForwardMap::NodeValue& result)
  {
    std::array<double, 3> fraction{};
    std::array<double, 3> value{};
    for (unsigned int axis = 0; axis < 3; ++axis)
    {
      fraction[axis] = (coordinate[axis] - bounds.lo[axis]) / (bounds.hi[axis] - bounds.lo[axis]);
      if (!std::isfinite(fraction[axis]) || fraction[axis] < 0. || fraction[axis] > 1.)
      {
        return false;
      }
    }
    for (unsigned int corner = 0; corner < 8; ++corner)
    {
      const auto* sample = corner_value(corner);
      if (!sample)
      {
        return false;
      }
      const double weight = ((corner & 4U) ? fraction[0] : 1. - fraction[0]) * ((corner & 2U) ? fraction[1] : 1. - fraction[1]) * ((corner & 1U) ? fraction[2] : 1. - fraction[2]);
      for (unsigned int component = 0; component < 3; ++component)
      {
        value[component] += weight * (*sample)[component];
      }
    }
    result = {true, value[0], value[1], value[2]};
    return true;
  }
  static std::string BuildCanonicalPolicy(const Data& d);
  static bool EstimatedBytes(std::uint64_t points, std::uint64_t cells,
                             std::uint64_t roots, std::uint64_t& bytes);
  // All budgets are charged per root independently of cross-root cache hits.
  // Samplerfalse is fatal; ordinary physical loss is a successful Point with
  // reaches0. Incomplete/budget/Mixed/AllLost leaves remain explicit direct.
  // Physical generation packs the completed state and writes the final compact
  // table, then independently reopens it before publication.
  static bool Generate(const PHGarfieldForwardTransport& garfield, const PHGarfieldForwardMap::Config& config,
                       const Policy& policy, std::uint64_t first,
                       std::uint64_t end, BuildSummary& build_summary,
                       std::ostream& log, std::ostream& error);
  static bool Build(const Geometry& geometry, const Policy& policy, std::uint64_t first,
                    std::uint64_t end, const Sampler& sampler, Data& output,
                    BuildSummary& summary, std::ostream& log, std::ostream& error);
};

#endif
