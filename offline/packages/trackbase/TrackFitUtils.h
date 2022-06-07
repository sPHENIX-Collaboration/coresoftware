#ifndef TRACKBASE_TRACKFITUTILS_H
#define TRACKBASE_TRACKFITUTILS_H

#include <Acts/Definitions/Algebra.hpp>

#include <tuple>
#include <utility>
#include <vector>

class TrackFitUtils
{

  public:

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
  static circle_fit_output_t circle_fit_by_taubin(const position_vector_t&);

  /// convenient overload
  static circle_fit_output_t circle_fit_by_taubin(const std::vector<Acts::Vector3>& );

  /// line fit output [slope, intercept]
  using line_fit_output_t = std::tuple<double,double>;

  /// line_fit
  /**
   * copied from: https://www.bragitoff.com
   * typically used to fit we want to fit z vs radius
   */
  static line_fit_output_t line_fit(const position_vector_t&);

  /// convenient overload
  static line_fit_output_t line_fit(const std::vector<Acts::Vector3>& );

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
  static circle_circle_intersection_output_t circle_circle_intersection( double r1, double r2, double x2, double y2 );  
  
};

#endif
