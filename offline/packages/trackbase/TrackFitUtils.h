#ifndef TRACKBASE_TRACKFITUTILS_H
#define TRACKBASE_TRACKFITUTILS_H

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

  /// line fit output [slope, intercept]
  using line_fit_output_t = std::tuple<double,double>;
  
  /// line_fit
  /**
   * copied from: https://www.bragitoff.com
   * typically used to fit we want to fit z vs radius
   */
  static line_fit_output_t line_fit(const position_vector_t&);

};

#endif
