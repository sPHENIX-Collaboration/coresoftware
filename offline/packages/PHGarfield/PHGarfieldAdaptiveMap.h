#ifndef PHGARFIELD_ADAPTIVEMAP_H
#define PHGARFIELD_ADAPTIVEMAP_H

#include "PHGarfieldForwardMap.h"

#include <array>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

class TDirectory;

// Schema4 deliberately owns primitive arrays, not live ROOT objects. Coordinates
// are(phi,r,signed z); corner order is(phi_bit<<2)|(r_bit<<1)|z_bit.
class PHGarfieldAdaptiveMap
{
 public:
  static constexpr int SchemaVersion = 4;
  static constexpr const char* ProbeScheme = "tensor_3x3x3_19_noncorner_actual_axis_curvature_v1";
  static constexpr const char* SplitPolicy = "largest_actual_midpoint_error_binary_axis_v1";
  static constexpr const char* AllLostPolicy = "direct_query_no_finite_probe_loss_certificate_v1";
  static constexpr const char* LeafEncoding = "0=Uncertified,1=Certified,2=Mixed,3=SampledAllLost,4=NumericallyUnresolved";
  static constexpr const char* FallbackEncoding = "0=None,1=PhysicalMixed,2=Numerical,3=DepthBudget,4=NodeBudget,5=QueryBudget,6=AngularBranch,7=SampledAllLost";
  // Engineering ceiling only. Actual vector growth remains governed by the
  // explicit per-file policy and checked actual point/cell counts.
  static constexpr std::uint64_t HardMaximumEstimatedBytes = 1024ULL * 1024 * 1024 * 1024;
  // Matches PHGarfield::ForwardDrift's radial_tolerance_cm. This is a numerical
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
    // Filled once after validation. It is derived, never persisted/trusted.
    std::vector<std::int64_t> root_lookup;
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

  static Geometry FromConfig(const PHGarfieldForwardMap::Config&);
  static double CanonicalPhi(double phi);
  static std::uint64_t RootCount(const Geometry&);
  static Bounds RootBounds(const Geometry&, std::uint64_t root_id);
  static Bounds ChildBounds(const Bounds&, int axis, unsigned int child);
  static std::array<double, 3> CornerCoordinate(const Bounds&, unsigned int corner);
  static bool Validate(Data&, std::ostream&, bool require_complete = true);
  // Shared by sampling and closed-file validation. Invalid samples fail before
  // further refinement; diagnostics preserve full-precision coordinates/values.
  static bool ValidatePoint(const Point&, const Geometry&, double maximum_drift_time_ns,
                            std::ostream&);
  static const Cell* FindLeaf(const Data&, unsigned int side, double phi, double r, double z);
  static bool Interpolate(const Data&, const Cell&, double phi, double r, double z,
                          PHGarfieldForwardMap::NodeValue&);
  // directory is GarfieldFullTransport/Adaptive. Read checks primitive branch
  // types, complete rows, content fingerprint, topology and declared policy;
  // it copies data and resets branch addresses before returning.
  static bool Read(TDirectory*, const Geometry&, Data&, std::ostream&,
                   bool require_complete = true);
  static bool Write(TDirectory*, const Data&, std::ostream&);
  static std::string ContentFingerprint(const Data&);
  static std::string BuildCanonicalPolicy(const Data&);
  static bool EstimatedBytes(std::uint64_t points, std::uint64_t cells,
                             std::uint64_t roots, std::uint64_t& bytes);
  // All budgets are charged per root independently of cross-root cache hits.
  // Samplerfalse is fatal; ordinary physical loss is a successful Point with
  // reaches0. Incomplete/budget/Mixed/AllLost leaves remain explicit direct.
  // Physical generator uses the identical preflight/canonical metadata helper
  // as legacy maps, then reopens all arrays before atomic publication.
  static bool Generate(const PHGarfield&, const PHGarfieldForwardMap::Config&,
                       const Policy&, std::uint64_t first_root,
                       std::uint64_t end_root, BuildSummary&,
                       std::ostream& log, std::ostream& error);
  static bool Build(const Geometry&, const Policy&, std::uint64_t first_root,
                    std::uint64_t end_root, const Sampler&, Data&,
                    BuildSummary&, std::ostream& log, std::ostream& error);
};

#endif
