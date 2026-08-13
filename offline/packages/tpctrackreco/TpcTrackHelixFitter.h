// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef TPCTRACKRECO_TPCTRACKHELIXFITTER_H
#define TPCTRACKRECO_TPCTRACKHELIXFITTER_H

#include "TpcTrackFit.h"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

class TpcTrackHelixFitter
{
 public:
  static bool parse_point_order(const std::string &mode, TpcTrackPointOrder &order);
  static void order_points(std::vector<TpcTrackPoint> &points, TpcTrackPointOrder order);

  static bool fit_circle_least_squares(const std::vector<TpcTrackPoint> &points,
                                       std::size_t nfit,
                                       double &cx,
                                       double &cy,
                                       double &radius);
  static bool fit_circle_curvature_refined(const std::vector<TpcTrackPoint> &points,
                                           std::size_t nfit,
                                           double &cx,
                                           double &cy,
                                           double &radius);
  static bool fit(const std::vector<TpcTrackPoint> &points,
                  int fit_first_points,
                  double bfield_t,
                  TpcTrackHelix &helix);
  static bool from_state(const TpcTrackVec3 &position,
                         const TpcTrackVec3 &momentum,
                         int charge,
                         double bfield_t,
                         TpcTrackHelix &helix);
  static bool orient_to_charge(TpcTrackHelix &helix, int charge);

  static TpcTrackVec3 point(const TpcTrackHelix &helix, double theta);
  static TpcTrackVec3 tangent(const TpcTrackHelix &helix, double theta);
  static TpcTrackVec3 momentum(const TpcTrackHelix &helix, double theta);

  static std::pair<double, double> theta_search_range(const TpcTrackHelix &helix,
                                                      double theta_extension,
                                                      double downstream_margin);
  static bool measurement_anchored_search_range(
      const TpcTrackHelix &helix,
      const std::vector<TpcTrackPoint> &points,
      double max_upstream_cm,
      double downstream_margin_cm,
      TpcTrackHelixSearchRange &range);
  static bool line_line_pca(const TpcTrackVec3 &pos1,
                            const TpcTrackVec3 &dir1,
                            const TpcTrackVec3 &pos2,
                            const TpcTrackVec3 &dir2,
                            TpcTrackLinePca &pca,
                            bool normalize_dirs);
  static TpcTrackHelixPca refine_pair(const TpcTrackHelix &helix1,
                                      const TpcTrackHelix &helix2,
                                      double theta1,
                                      double theta2,
                                      double min1,
                                      double max1,
                                      double min2,
                                      double max2,
                                      double max_step);
  static std::vector<TpcTrackHelixPca> pca_candidates(const TpcTrackHelix &helix1,
                                                      const TpcTrackHelix &helix2,
                                                      double theta_extension,
                                                      int coarse_steps,
                                                      double downstream_margin,
                                                      int max_candidates);
  static std::vector<TpcTrackHelixPca> pca_candidates_in_ranges(
      const TpcTrackHelix &helix1,
      const TpcTrackHelix &helix2,
      const TpcTrackHelixSearchRange &range1,
      const TpcTrackHelixSearchRange &range2,
      int coarse_steps,
      int max_candidates);

  static std::pair<double, double> helix_dca_to_vertex(const TpcTrackHelix &helix,
                                                       const TpcTrackVec3 &vertex);
  static std::pair<double, double> line_dca_to_vertex(const TpcTrackVec3 &pos,
                                                      const TpcTrackVec3 &mom,
                                                      const TpcTrackVec3 &vertex);

  static bool finite(const TpcTrackVec3 &value);
  static TpcTrackVec3 add(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs);
  static TpcTrackVec3 subtract(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs);
  static TpcTrackVec3 scale(const TpcTrackVec3 &value, double factor);
  static double dot(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs);
  static double norm(const TpcTrackVec3 &value);
  static TpcTrackVec3 unit(const TpcTrackVec3 &value);
  static double pt(const TpcTrackVec3 &value);
  static double distance(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs);
  static double vector_cosine(const TpcTrackVec3 &lhs, const TpcTrackVec3 &rhs);
};

#endif
