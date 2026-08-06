// Behavioral test for PhotonClusterBuilder::topocluster_raw_iso and the
// exact-first feature resolver used by calculate_bdt_score.
//
// Run manually after building the package:
//   ./testPhotonClusterBuilderTopoIso
// Exit code 0 means every check passed.

#include "PhotonClusterBuilder.h"

#include <calobase/PhotonClusterv1.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterv1.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace
{
  int failures = 0;

  void check(bool passed, const std::string& label)
  {
    if (!passed)
    {
      ++failures;
      std::cerr << "FAIL: " << label << std::endl;
    }
    else
    {
      std::cout << "ok:   " << label << std::endl;
    }
  }

  void check_close(float value, float expected, const std::string& label, float tolerance = 1e-5F)
  {
    const bool passed = std::isfinite(value) && (std::fabs(value - expected) <= tolerance);
    if (!passed)
    {
      ++failures;
      std::cerr << "FAIL: " << label << " (value=" << value << " expected=" << expected << ")" << std::endl;
    }
    else
    {
      std::cout << "ok:   " << label << std::endl;
    }
  }

  // Place a topo cluster on a cylinder of radius r_cm so that, measured from a
  // vertex at (0, 0, vertex_z=0), its pseudorapidity is exactly eta.
  void add_topo(RawClusterContainer* container, float eta, float phi, float energy, float r_cm = 100.0F)
  {
    auto* topo = new RawClusterv1();
    topo->set_r(r_cm);
    topo->set_phi(phi);
    topo->set_z(r_cm * std::sinh(eta));
    topo->set_energy(energy);
    container->AddCluster(topo);
  }

  float et_of(float energy, float eta)
  {
    return energy / std::cosh(eta);
  }
}  // namespace

int main()
{
  const float nan_value = std::numeric_limits<float>::quiet_NaN();

  // --- Signed sums, exclusion, nesting, single subtraction ------------------
  {
    RawClusterContainer topoclusters;
    const float seed_eta = 0.10F;
    const float seed_phi = 0.50F;
    const float candidate_et = 3.0F;

    // In-cone at the axis (the candidate's own topo cluster stays in the sum;
    // the single -candidate_et subtraction is the only removal).
    add_topo(&topoclusters, 0.10F, 0.50F, 2.0F);
    // In-cone (dR = 0.20) with negative signed energy: must contribute with sign.
    add_topo(&topoclusters, 0.30F, 0.50F, -0.5F);
    // Out of every cone (dR = 0.70): must never contribute.
    add_topo(&topoclusters, 0.10F, 1.20F, 5.0F);
    // Between the nested radii (dR = 0.35): only the R=0.4 cone sees it.
    add_topo(&topoclusters, 0.45F, 0.50F, 1.0F);

    const float sum_03 = et_of(2.0F, 0.10F) + et_of(-0.5F, 0.30F);
    const float sum_04 = sum_03 + et_of(1.0F, 0.45F);

    const float iso_03 = PhotonClusterBuilder::topocluster_raw_iso(&topoclusters, seed_eta, seed_phi, 0.3F, 0.0F, candidate_et);
    const float iso_04 = PhotonClusterBuilder::topocluster_raw_iso(&topoclusters, seed_eta, seed_phi, 0.4F, 0.0F, candidate_et);

    check_close(iso_03, sum_03 - candidate_et, "R=0.3 signed in-cone sum minus candidate ET exactly once");
    check_close(iso_04, sum_04 - candidate_et, "R=0.4 includes the dR=0.35 cluster the R=0.3 cone excludes");
    const float positive_only_sum_03 = et_of(2.0F, 0.10F);
    check(iso_03 < positive_only_sum_03 - candidate_et - 1e-4F,
          "negative topo-cluster ET is retained with sign, not clipped");
    check_close(iso_04 - iso_03, et_of(1.0F, 0.45F), "nested-cone difference is exactly the shell cluster ET");
  }

  // --- Phi wrapping across +/- pi -------------------------------------------
  {
    RawClusterContainer topoclusters;
    const float seed_phi = 3.10F;
    // Across the boundary: delta-phi wraps to 2*pi - 6.20 ~= 0.083, inside R=0.3.
    add_topo(&topoclusters, 0.0F, -3.10F, 1.0F);
    // Same side, |delta-phi| = 0.35: outside R=0.3, inside R=0.4.
    add_topo(&topoclusters, 0.0F, 2.75F, 1.0F);

    const float iso_03 = PhotonClusterBuilder::topocluster_raw_iso(&topoclusters, 0.0F, seed_phi, 0.3F, 0.0F, 1.0F);
    const float iso_04 = PhotonClusterBuilder::topocluster_raw_iso(&topoclusters, 0.0F, seed_phi, 0.4F, 0.0F, 1.0F);

    check_close(iso_03, et_of(1.0F, 0.0F) - 1.0F, "delta-phi wraps across +/- pi into the cone");
    check_close(iso_04, 2.0F * et_of(1.0F, 0.0F) - 1.0F, "unwrapped same-side cluster obeys the plain cone edge");
  }

  // --- Mirrored wrap direction (seed below -pi side) ------------------------
  {
    RawClusterContainer topoclusters;
    // Seed at -3.10, cluster at +3.10: delta-phi wraps to ~ -0.083, inside R=0.3.
    add_topo(&topoclusters, 0.0F, 3.10F, 1.0F);
    const float iso_03 = PhotonClusterBuilder::topocluster_raw_iso(&topoclusters, 0.0F, -3.10F, 0.3F, 0.0F, 1.0F);
    check_close(iso_03, et_of(1.0F, 0.0F) - 1.0F, "delta-phi wraps across +/- pi in the mirrored direction");
  }

  // --- Vertex-relative kinematics -------------------------------------------
  {
    RawClusterContainer topoclusters;
    const float vertex_z = 30.0F;
    const float z_topo = 100.0F * std::sinh(0.10F);
    add_topo(&topoclusters, 0.10F, 0.50F, 2.0F);

    const float eta_seen = std::asinh((z_topo - vertex_z) / 100.0F);
    const float iso = PhotonClusterBuilder::topocluster_raw_iso(&topoclusters, eta_seen, 0.50F, 0.3F, vertex_z, 1.0F);
    check_close(iso, et_of(2.0F, eta_seen) - 1.0F, "topo-cluster eta/ET are measured from the event vertex");
  }

  // --- Empty container is a real (zero-activity) measurement ----------------
  {
    RawClusterContainer topoclusters;
    const float iso = PhotonClusterBuilder::topocluster_raw_iso(&topoclusters, 0.0F, 0.0F, 0.4F, 0.0F, 2.5F);
    check_close(iso, -2.5F, "empty topo container yields -candidate_et, not an error");
  }

  // --- Unusable inputs must be NaN, never a plausible physics value ---------
  {
    RawClusterContainer topoclusters;
    add_topo(&topoclusters, 0.0F, 0.0F, 1.0F);

    check(std::isnan(PhotonClusterBuilder::topocluster_raw_iso(nullptr, 0.0F, 0.0F, 0.4F, 0.0F, 1.0F)),
          "missing topo container yields NaN");
    check(std::isnan(PhotonClusterBuilder::topocluster_raw_iso(&topoclusters, nan_value, 0.0F, 0.4F, 0.0F, 1.0F)),
          "non-finite seed eta yields NaN");
    check(std::isnan(PhotonClusterBuilder::topocluster_raw_iso(&topoclusters, 0.0F, 0.0F, 0.0F, 0.0F, 1.0F)),
          "non-positive radius yields NaN");
    check(std::isnan(PhotonClusterBuilder::topocluster_raw_iso(&topoclusters, 0.0F, 0.0F, 0.4F, nan_value, 1.0F)),
          "non-finite vertex z yields NaN");
    check(std::isnan(PhotonClusterBuilder::topocluster_raw_iso(&topoclusters, 0.0F, 0.0F, 0.4F, 0.0F, nan_value)),
          "non-finite candidate ET yields NaN");
  }

  // --- Exact-first BDT feature resolution -----------------------------------
  {
    PhotonClusterv1 photon;
    photon.set_energy(9.0F);
    photon.set_shower_shape_parameter("cluster_eta", 0.42F);
    photon.set_shower_shape_parameter("cluster_phi", -0.75F);
    photon.set_shower_shape_parameter("cluster_pt", 7.5F);
    photon.set_shower_shape_parameter("cluster_ietacent", 41.0F);
    photon.set_shower_shape_parameter("cluster_iphicent", 173.0F);
    photon.set_shower_shape_parameter("eta", -9.0F);
    photon.set_shower_shape_parameter("phi", -9.0F);
    photon.set_shower_shape_parameter("pt", -9.0F);
    photon.set_shower_shape_parameter("ietacent", -9.0F);
    photon.set_shower_shape_parameter("iphicent", -9.0F);
    photon.set_shower_shape_parameter("weta", 0.18F);
    photon.set_shower_shape_parameter("e33", 6.0F);

    const float vertex_z = 12.5F;

    check_close(PhotonClusterBuilder::resolve_bdt_feature(&photon, "cluster_eta", vertex_z), 0.42F,
                "exact cluster_eta wins before cluster_ prefix fallback");
    check_close(PhotonClusterBuilder::resolve_bdt_feature(&photon, "cluster_phi", vertex_z), -0.75F,
                "exact cluster_phi wins before cluster_ prefix fallback");
    check_close(PhotonClusterBuilder::resolve_bdt_feature(&photon, "cluster_pt", vertex_z), 7.5F,
                "exact cluster_pt wins before cluster_ prefix fallback");
    check_close(PhotonClusterBuilder::resolve_bdt_feature(&photon, "cluster_ietacent", vertex_z), 41.0F,
                "exact cluster_ietacent wins before cluster_ prefix fallback");
    check_close(PhotonClusterBuilder::resolve_bdt_feature(&photon, "cluster_iphicent", vertex_z), 173.0F,
                "exact cluster_iphicent wins before cluster_ prefix fallback");

    check_close(PhotonClusterBuilder::resolve_bdt_feature(&photon, "cluster_Eta", vertex_z), 0.42F,
                "PPG12/TMVA cluster_Eta alias resolves after exact miss");
    check_close(PhotonClusterBuilder::resolve_bdt_feature(&photon, "cluster_Phi", vertex_z), -0.75F,
                "PPG12/TMVA cluster_Phi alias resolves after exact miss");
    check_close(PhotonClusterBuilder::resolve_bdt_feature(&photon, "cluster_Et", vertex_z),
                9.0F / std::cosh(0.42F), "PPG12/TMVA cluster_Et alias resolves after exact miss");
    check_close(PhotonClusterBuilder::resolve_bdt_feature(&photon, "vertexz", vertex_z), vertex_z,
                "vertexz alias resolves after exact miss");
    check_close(PhotonClusterBuilder::resolve_bdt_feature(&photon, "cluster_weta", vertex_z), 0.18F,
                "cluster_ prefix fallback resolves an unprefixed shower shape");
    check_close(PhotonClusterBuilder::resolve_bdt_feature(&photon, "cluster_weta_over_e33", vertex_z), 0.03F,
                "ratio inputs use the same alias and prefix resolution contract");

    photon.set_shower_shape_parameter("cluster_Eta", 0.99F);
    check_close(PhotonClusterBuilder::resolve_bdt_feature(&photon, "cluster_Eta", vertex_z), 0.99F,
                "an exact stored key wins even when its spelling is also a known alias");
    check(std::isnan(PhotonClusterBuilder::resolve_bdt_feature(&photon, "cluster_missing", vertex_z)),
          "a missing exact key and missing fallback remain NaN");
    check(std::isnan(PhotonClusterBuilder::resolve_bdt_feature(nullptr, "cluster_eta", vertex_z)),
          "a null photon returns NaN");
  }

  if (failures)
  {
    std::cerr << failures << " check(s) failed" << std::endl;
    return 1;
  }
  std::cout << "all checks passed" << std::endl;
  return 0;
}
