#include "TrackFitUtils.h"

#include <cmath>

namespace
{

  //! convenience square method
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }
}

//_________________________________________________________________________________
TrackFitUtils::circle_fit_output_t TrackFitUtils::circle_fit_by_taubin( const TrackFitUtils::position_vector_t& positions )
{

  // Compute x- and y- sample means
  double meanX = 0;
  double meanY = 0;
  double weight = 0;

  for( const auto& [x,y] : positions)
  {
    meanX += x;
    meanY += y;
    ++weight;
  }
  meanX /= weight;
  meanY /= weight;

  //     computing moments

  double Mxx = 0;
  double Myy = 0;
  double Mxy = 0;
  double Mxz = 0;
  double Myz = 0;
  double Mzz = 0;

  for (auto& [x,y] : positions)
  {
    double Xi = x - meanX;   //  centered x-coordinates
    double Yi = y - meanY;   //  centered y-coordinates
    double Zi = square(Xi) + square(Yi);

    Mxy += Xi*Yi;
    Mxx += Xi*Xi;
    Myy += Yi*Yi;
    Mxz += Xi*Zi;
    Myz += Yi*Zi;
    Mzz += Zi*Zi;
  }
  Mxx /= weight;
  Myy /= weight;
  Mxy /= weight;
  Mxz /= weight;
  Myz /= weight;
  Mzz /= weight;

  //  computing coefficients of the characteristic polynomial

  const double Mz = Mxx + Myy;
  const double Cov_xy = Mxx*Myy - Mxy*Mxy;
  const double Var_z = Mzz - Mz*Mz;
  const double A3 = 4*Mz;
  const double A2 = -3*Mz*Mz - Mzz;
  const double A1 = Var_z*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
  const double A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
  const double A22 = A2 + A2;
  const double A33 = A3 + A3 + A3;

  //    finding the root of the characteristic polynomial
  //    using Newton's method starting at x=0
  //    (it is guaranteed to converge to the right root)
  static constexpr int iter_max = 99;
  double x = 0;
  double y = A0;

  // usually, 4-6 iterations are enough
  for( int iter=0; iter<iter_max; ++iter)
  {
    const double Dy = A1 + x*(A22 + A33*x);
    const double xnew = x - y/Dy;
    if ((xnew == x)||(!std::isfinite(xnew))) break;

    const double ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
    if (std::abs(ynew)>=std::abs(y))  break;

    x = xnew;  y = ynew;

  }

  //  computing parameters of the fitting circle
  const double DET = square(x) - x*Mz + Cov_xy;
  const double Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2;
  const double Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2;

  //  assembling the output
  double X0 = Xcenter + meanX;
  double Y0 = Ycenter + meanY;
  double R = std::sqrt( square(Xcenter) + square(Ycenter) + Mz);
  return std::make_tuple( R, X0, Y0 );

}

//_________________________________________________________________________________
TrackFitUtils::circle_fit_output_t TrackFitUtils::circle_fit_by_taubin( const std::vector<Acts::Vector3>& positions )
{
  position_vector_t positions_2d;
  for( const auto& position:positions )
  { positions_2d.emplace_back( position.x(), position.y() ); }

  return circle_fit_by_taubin( positions_2d );
}

//_________________________________________________________________________________
TrackFitUtils::line_fit_output_t TrackFitUtils::line_fit(const TrackFitUtils::position_vector_t& positions)
{
  double xsum=0;
  double x2sum=0;
  double ysum=0;
  double xysum=0;
  for( const auto& [r,z]:positions )
  {
    xsum=xsum+r;                        //calculate sigma(xi)
    ysum=ysum+z;                        //calculate sigma(yi)
    x2sum=x2sum+square(r);              //calculate sigma(x^2i)
    xysum=xysum+r*z;                    //calculate sigma(xi*yi)
  }

  const auto npts = positions.size();
  const double denominator = (x2sum*npts-square(xsum));
  const double a= (xysum*npts-xsum*ysum)/denominator;            //calculate slope
  const double b= (x2sum*ysum-xsum*xysum)/denominator;           //calculate intercept
  return std::make_tuple( a, b );
}

//_________________________________________________________________________________
TrackFitUtils::line_fit_output_t TrackFitUtils::line_fit( const std::vector<Acts::Vector3>& positions )
{
  position_vector_t positions_2d;
  for( const auto& position:positions )
  { positions_2d.emplace_back( std::sqrt( square(position.x())+square(position.y())), position.z() ); }

  return line_fit( positions_2d );
}

//_________________________________________________________________________________
TrackFitUtils::circle_circle_intersection_output_t TrackFitUtils::circle_circle_intersection(double r1, double r2, double x2, double y2 )
{
  
  const double D = square(r1) - square(r2) + square(x2) + square(y2);
  const double a = 1.0 + square(x2/y2);
  const double b = - D * x2/square(y2);
  const double c = square(D/(2.*y2)) - square(r1);
  const double delta = square(b)-4.*a*c;
  
  const double sqdelta = std::sqrt( delta );
  
  const double xplus = (-b + sqdelta ) / (2. * a);
  const double xminus = (-b - sqdelta ) / (2. * a);
  
  const double yplus  = -(2*x2*xplus - D) / (2.*y2);
  const double yminus = -(2*x2*xminus - D) / (2.*y2);

  
  return std::make_tuple( xplus, yplus, xminus, yminus );

}
