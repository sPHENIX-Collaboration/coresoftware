#ifndef TRACKBASE_TRACKFITUTILS_H
#define TRACKBASE_TRACKFITUTILS_H

#include "ActsSurfaceMaps.h"
#include "TrkrDefs.h"

#include <Acts/Definitions/Algebra.hpp>

#include <tuple>
#include <utility>
#include <vector>

class ActsGeometry;
class TrkrClusterContainer;

namespace TrackFitUtils
{
  using position_t = std::pair<double, double>;
  using position_vector_t = std::vector<position_t>;

  /// circle fit output [R, x0, y0]
  using circle_fit_output_t = std::tuple<double, double, double>;

  /**
   * Circle fit to a given set of data points (in 2D)
   * This is an algebraic fit, due to Taubin, based on the journal article
   * G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
   * Space Curves Defined By Implicit Equations, With
   * Applications To Edge And Range Image Segmentation",
   * IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
   * It works well whether data points are sampled along an entire circle or along a small arc.
   * It still has a small bias and its statistical accuracy is slightly lower than that of the geometric fit (minimizing geometric distances),
   * It provides a very good initial guess for a subsequent geometric fit.
   * Nikolai Chernov  (September 2012)
   */
  circle_fit_output_t circle_fit_by_taubin(const position_vector_t&);

  /// convenient overload
  circle_fit_output_t circle_fit_by_taubin(const std::vector<Acts::Vector3>&);

  /// line fit output [slope, intercept]
  using line_fit_output_t = std::tuple<double, double>;

  /// line_fit
  /**
   * copied from: https://www.bragitoff.com
   * typically used to fit we want to fit z vs radius
   * 
   * Updated to use the "Deming model" by default:
   * minimizing by orthoncal distance to line in x and y
   * (instead of the y-distance). For details, see:
   * http://staff.pubhealth.ku.dk/~bxc/MethComp/Deming.pdf
   */
  line_fit_output_t line_fit(const position_vector_t&);

  /*
   * Need to make a metric for distance from points to lines origin (pca).
   *  - project point "global" to the line.
   *  - return distance on line to the pca (the point of closest approach to origin)
   */
  double line_dist_to_pca (const double slope, const double intercept, 
      const Acts::Vector2& pca, const Acts::Vector3& global);

  /// convenient overload
  line_fit_output_t line_fit(const std::vector<Acts::Vector3>&);

  line_fit_output_t line_fit_xy(const std::vector<Acts::Vector3>& positions);
  line_fit_output_t line_fit_xz(const std::vector<Acts::Vector3>& positions);

  /// line-circle intersection output. (xplus, yplus, xminus, yminus)
  using line_circle_intersection_output_t = std::tuple<double, double, double, double>;
  /**
  * r is radius of sPHENIX layer
  * m and b are parameters of line fitted to TPC clusters in the x-y plane (slope and intersection)
  * the solutions are xplus, xminus, yplus, yminus
  * The intersection of the line-circle occurs when
  * y = m*x + b
  * Here we assume that circle is an sPHENIX layer centered on x1=y1=0
  * x^2 +y^2 = r^2
  * substitute for y in equation of circle, solve for x, and then for y.
  **/
  line_circle_intersection_output_t line_circle_intersection(double r, double m, double b);

  /// circle-circle intersection output. (xplus, yplus, xminus, yminus)
  using circle_circle_intersection_output_t = std::tuple<double, double, double, double>;

  /**
   * r1 is radius of sPHENIX layer
   * r2, x2 and y2 are parameters of circle fitted to TPC clusters
   * the solutions are xplus, xminus, yplus, yminus

   * The intersection of two circles occurs when
   * (x-x1)^2 + (y-y1)^2 = r1^2,  / (x-x2)^2 + (y-y2)^2 = r2^2
   * Here we assume that circle 1 is an sPHENIX layer centered on x1=y1=0, and circle 2 is arbitrary
   * x^2 +y^2 = r1^2,   (x-x2)^2 + (y-y2)^2 = r2^2
   * expand the equations and subtract to eliminate the x^2 and y^2 terms, gives the radical line connecting the intersection points
   * iy = - (2*x2*x - D) / 2*y2,
   * then substitute for y in equation of circle 1
   */
  circle_circle_intersection_output_t circle_circle_intersection(double r1, double r2, double x2, double y2);

  unsigned int addClusters(std::vector<float>& fitpars,
                                  double dca_cut,
                                  ActsGeometry* _tGeometry,
                                  TrkrClusterContainer* _cluster_map,
                                  std::vector<Acts::Vector3>& global_vec,
                                  std::vector<TrkrDefs::cluskey>& cluskey_vec,
                                  unsigned int startLayer,
                                  unsigned int endLayer);

  unsigned int addClustersOnLine(TrackFitUtils::line_fit_output_t& fitpars,
                                        const bool& isXY,
                                        const double& dca_cut,
                                        ActsGeometry* tGeometry,
                                        TrkrClusterContainer* clusterContainer,
                                        std::vector<Acts::Vector3>& global_vec,
                                        std::vector<TrkrDefs::cluskey>& cluskey_vec,
                                        const unsigned int& startLayer,
                                        const unsigned int& endLayer);

  std::pair<Acts::Vector3, Acts::Vector3> get_helix_tangent(const std::vector<float>& fitpars, Acts::Vector3& global);

  Acts::Vector3 get_helix_pca(std::vector<float>& fitpars, const Acts::Vector3& global);

  Acts::Vector2 get_circle_point_pca(float radius, float x0, float y0, Acts::Vector3 global);

  std::vector<float> fitClusters(std::vector<Acts::Vector3>& global_vec,
                                        std::vector<TrkrDefs::cluskey> cluskey_vec,
                                        bool use_intt = false);
  void getTrackletClusters(ActsGeometry* _tGeometry,
                                  TrkrClusterContainer* _cluster_map,
                                  std::vector<Acts::Vector3>& global_vec,
                                  const std::vector<TrkrDefs::cluskey>& cluskey_vec);
  Acts::Vector3 surface_3Dline_intersection(const TrkrDefs::cluskey& key,
                                                   TrkrCluster* cluster, ActsGeometry* geometry, float& xyslope, float& xyint, float& yzslope, float& yzint);

  Acts::Vector3 getPCALinePoint(const Acts::Vector3& global, const Acts::Vector3& tangent, const Acts::Vector3& posref);
  Acts::Vector3 get_helix_surface_intersection(const Surface& surf,
                                                      std::vector<float>& fitpars, Acts::Vector3& global, ActsGeometry* _tGeometry);
  Acts::Vector3 get_line_plane_intersection(const Acts::Vector3& pca, const Acts::Vector3& tangent,
                                                   const Acts::Vector3& sensorCenter, const Acts::Vector3& sensorNormal);

  std::vector<double> getLineClusterResiduals(position_vector_t& rz_pts, float slope, float intercept);
  std::vector<double> getCircleClusterResiduals(position_vector_t& xy_pts, float R, float X0, float Y0);

  Acts::Vector2 get_line_point_pca(double slope, double intercept, Acts::Vector3 global);
  std::vector<float> fitClustersZeroField(std::vector<Acts::Vector3>& global_vec,
						       std::vector<TrkrDefs::cluskey> cluskey_vec, bool use_intt);

  float get_helix_pathlength(std::vector<float>& fitpars, const Acts::Vector3& start_point, const Acts::Vector3& end_point);
  float get_helix_surface_pathlength(const Surface& surf, std::vector<float>& fitpars, const Acts::Vector3& start_point, ActsGeometry* tGeometry);

  /* std::tuple<double,double> dca_on_line2D_to_point(const double x0, const double y0, const double xy_m, const double xy_b); */
  // get track fit parameters for track matching:: is-good-fit, phi, eta, pt==1, pos_vec3, mom_vec3
  std::tuple<bool, double, double, double, Acts::Vector3, Acts::Vector3> 
    zero_field_track_params(
      ActsGeometry* _tGeometry, 
      TrkrClusterContainer* _cluster_map, 
      const std::vector<TrkrDefs::cluskey>& clusters
    );

   double z_fit_to_pca(const double slope, const double intercept, 
    const std::vector<Acts::Vector3>& glob_pts);
};

#endif
