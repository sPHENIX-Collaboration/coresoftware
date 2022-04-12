#include "TrackSeed_v1.h"
#include "ActsTransformations.h"

TrackSeed_v1::TrackSeed_v1()
{}

TrackSeed_v1::TrackSeed_v1(const TrackSeed& seed)
{ TrackSeed_v1::CopyFrom( seed ); }

// have to suppress uninitMenberVar from cppcheck since it triggers many false positive
// cppcheck-suppress uninitMemberVar
TrackSeed_v1::TrackSeed_v1(const TrackSeed_v1& seed)
{ TrackSeed_v1::CopyFrom( seed ); }

TrackSeed_v1& TrackSeed_v1::operator=(const TrackSeed_v1& seed)
{ if( this != &seed ) CopyFrom( seed ); return *this; }

TrackSeed_v1::~TrackSeed_v1()
{}

void TrackSeed_v1::CopyFrom( const TrackSeed& seed )
{
  if( this == &seed ) return;
  TrackSeed::CopyFrom( seed );

  m_qOverR = seed.get_qOverR();
  m_X0 = seed.get_X0();
  m_Y0 = seed.get_Y0();
  m_slope = seed.get_slope();
  m_Z0 = seed.get_Z0();
  m_crossing = seed.get_crossing();

  m_cluster_keys.clear();
  std::copy(seed.begin_cluster_keys(), seed.end_cluster_keys(), 
	    std::inserter(m_cluster_keys, m_cluster_keys.begin() ) );

}


void TrackSeed_v1::identify(std::ostream& os) const
{
  os << "TrackSeed_v1 object ";
  os << "charge " << get_charge() << std::endl;
  os << "(px,py,pz) = (" << get_px() << ", " << get_py()
     << ", " << get_pz() << ")" << std::endl;
  os << "(x,y,z) = (" << get_x() << ", " << get_y() << ", " << get_z() 
     << ")" << std::endl;
  os << "list of cluster keys ";
  if(m_cluster_keys.size() > 0) 
    {
      for (TrackSeed::ConstClusterKeyIter iter = begin_cluster_keys();
	   iter != end_cluster_keys();
	   ++iter)
	{
	  TrkrDefs::cluskey cluster_key = *iter;
	  os << cluster_key << ", ";
	}
    }
  
  os << std::endl;
  return;
}



void TrackSeed_v1::circleFitByTaubin(TrkrClusterContainer *clusters,
					 ActsSurfaceMaps* surfmaps,
					 ActsTrackingGeometry *tGeometry)
{
  /**  
   *   Circle fit to a given set of data points (in 2D)
   *   This is an algebraic fit, due to Taubin, based on the journal article
   *   G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
   *               Space Curves Defined By Implicit Equations, With 
   *               Applications To Edge And Range Image Segmentation",
   *               IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)
   *  It works well whether data points are sampled along an entire circle 
   *  or along a small arc. 
   *  It still has a small bias and its statistical accuracy is slightly lower 
   *  than that of the geometric fit (minimizing geometric distances),
   *  It provides a very good initial guess for a subsequent geometric fit. 
   *    Nikolai Chernov  (September 2012)
   */
  
  int iter, IterMAX = 99;
  
  double Mz, Mxy, Mxx, Myy, Mxz, Myz, Mzz, Cov_xy, Var_z;
  double A0, A1, A2, A22, A3, A33;
  double x, y;
  double DET, Xcenter, Ycenter;
  
  // Compute x- and y- sample means   
  double meanX = 0;
  double meanY = 0;
  double weight = 0;

  ActsTransformations transformer;
  std::vector<Acts::Vector3> globalPositions;
  for(const auto& key: m_cluster_keys)
    {
      Acts::Vector3 globalPos = transformer.getGlobalPosition(
			          clusters->findCluster(key),
				  surfmaps, tGeometry);
      globalPositions.push_back(globalPos);
      meanX += globalPos(0);
      meanY += globalPos(1);
      weight++;
    }

  meanX /= weight;
  meanY /= weight;

  Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;

  for(const auto& pos : globalPositions)
    {
      
      double Xi = pos(0) - meanX;
      double Yi = pos(1) - meanY;
      double Zi = Xi * Xi + Yi * Yi;

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

  Mz = Mxx + Myy;
  Cov_xy = Mxx * Myy - Mxy * Mxy;
  Var_z = Mzz - Mz * Mz;
  A3 = 4 * Mz;
  A2 = -3 * Mz * Mz - Mzz;
  A1 = Var_z * Mz + 4 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz;
  A0 = Mxz * (Mxz * Myy - Myz * Mxy) + 
    Myz * (Myz * Mxx - Mxz * Mxy) - Var_z * Cov_xy;
  A22 = A2 + A2;
  A33 = A3 + A3 + A3;

  for (x=0., y=A0, iter=0; iter<IterMAX; iter++)  // usually, 4-6 iterations are enough
    {
      double Dy = A1 + x * (A22 + A33 * x);
      double xnew = x - y / Dy;
      if ((xnew == x)||(!std::isfinite(xnew))) break;
      double ynew = A0 + xnew * (A1 + xnew * (A2 + xnew * A3));
      if (fabs(ynew)>=fabs(y))  break;
      x = xnew;  y = ynew;
    }
  
  //  computing parameters of the fitting circle
  
  DET = x*x - x*Mz + Cov_xy;
  Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/2;
  Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2;
  
  //  Update the circle fit parameters
  m_X0 = Xcenter + meanX;
  m_Y0 = Ycenter + meanY;
  m_qOverR = 1. / sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz);

  /// Set the charge
  Acts::Vector3 firstpos = globalPositions.at(0);
  Acts::Vector3 secondpos = globalPositions.at(1);
  
  float firstphi = atan2(firstpos(1), firstpos(0));
  float secondphi = atan2(secondpos(1), secondpos(0));
  float dphi = secondphi - firstphi;
  if(dphi > M_PI) dphi = 2.*M_PI - dphi;
  if(dphi < -M_PI) dphi = 2*M_PI + dphi;
  if(dphi<0) m_qOverR *= -1;

}

void TrackSeed_v1::lineFit(TrkrClusterContainer *clusters,
			       ActsSurfaceMaps *surfMaps, 
			       ActsTrackingGeometry *tGeometry)
{
  // copied from: https://www.bragitoff.com
  // we want to fit z vs radius
  double xsum = 0,x2sum = 0,ysum = 0,xysum = 0;    
  ActsTransformations transformer;
  for(const auto& key : m_cluster_keys)
    {
      Acts::Vector3 pos = transformer.getGlobalPosition(clusters->findCluster(key),
							surfMaps, tGeometry);
      double z = pos(2);
      double r = sqrt(pow(pos(0),2) + pow(pos(1), 2));
      
      xsum=xsum+r;               // calculate sigma(xi)
      ysum=ysum+z;               // calculate sigma(yi)
      x2sum=x2sum+pow(r,2);      // calculate sigma(x^2i)
      xysum=xysum+r*z;           // calculate sigma(xi*yi)
    
    }
  
  /// calculate slope
  m_slope = (m_cluster_keys.size()*xysum-xsum*ysum) / (m_cluster_keys.size()*x2sum-xsum*xsum);

  /// calculate intercept
  m_Z0 = (x2sum*ysum-xsum*xysum) / (x2sum*m_cluster_keys.size()-xsum*xsum);
  
}   

void TrackSeed_v1::findRoot(float& x, float& y) const
{
  /**
   * We need to determine the closest point on the circle to the origin
   * since we can't assume that the track originates from the origin
   * The eqn for the circle is (x-X0)^2+(y-Y0)^2=R^2 and we want to 
   * minimize d = sqrt((0-x)^2+(0-y)^2), the distance between the 
   * origin and some (currently, unknown) point on the circle x,y.
   * 
   * Solving the circle eqn for x and substituting into d gives an eqn for
   * y. Taking the derivative and setting equal to 0 gives the following 
   * two solutions. We take the smaller solution as the correct one, as 
   * usually one solution is wildly incorrect (e.g. 1000 cm)
   */
  float R = fabs(1./m_qOverR);
  double miny = (sqrt(pow(m_X0, 2) * pow(R, 2) * pow(m_Y0, 2) + pow(R, 2) 
		      * pow(m_Y0,4)) + pow(m_X0,2) * m_Y0 + pow(m_Y0, 3)) 
    / (pow(m_X0, 2) + pow(m_Y0, 2));

  double miny2 = (-sqrt(pow(m_X0, 2) * pow(R, 2) * pow(m_Y0, 2) + pow(R, 2) 
		      * pow(m_Y0,4)) + pow(m_X0,2) * m_Y0 + pow(m_Y0, 3)) 
    / (pow(m_X0, 2) + pow(m_Y0, 2));

  double minx = sqrt(pow(R, 2) - pow(miny - m_Y0, 2)) + m_X0;
  double minx2 = -sqrt(pow(R, 2) - pow(miny2 - m_Y0, 2)) + m_X0;
  

  /// Figure out which of the two roots is actually closer to the origin
  if(fabs(minx) < fabs(minx2))
    { x = minx; }
  else
    { x = minx2; }

  if(fabs(miny) < fabs(miny2))
    { y = miny; }
  else
    { y = miny2; }
  
}
float TrackSeed_v1::findRoot(bool findX) const
{
  float x=NAN, y=NAN;
  findRoot(x,y);
  if(findX) 
    { return x; }

  return y;
}

float TrackSeed_v1::get_x() const
{
  return findRoot(true);
}

float TrackSeed_v1::get_y() const
{
  return findRoot(false);
}

float TrackSeed_v1::get_z() const
{
  return get_Z0();
}

float TrackSeed_v1::get_pt() const
{
  /// Scaling factor for radius in 1.4T field
  return 0.3 * 1.4 / 100. * fabs(1./m_qOverR);
}

float TrackSeed_v1::get_phi() const
{
  float x=NAN, y=NAN;
  findRoot(x,y);
  return atan2(-1 * (m_X0-x), m_Y0-y);
}

float TrackSeed_v1::get_theta() const
{
  float theta = atan(1./m_slope);
  /// Normalize to 0<theta<pi
  if(theta < 0) 
    { theta += M_PI; }
  return theta;
}

float TrackSeed_v1::get_eta() const
{
  return -log(tan(get_theta() / 2.));
}

float TrackSeed_v1::get_p() const
{
  return get_pt() * cosh(get_eta());
}

float TrackSeed_v1::get_px() const
{
  return get_p() * sin(get_theta()) * cos(get_phi());
}

float TrackSeed_v1::get_py() const
{
  return get_p() * sin(get_theta()) * sin(get_phi());
}

float TrackSeed_v1::get_pz() const
{
  return get_p() * cos(get_theta());
}

int TrackSeed_v1::get_charge() const
{
  return ( m_qOverR < 0 ) ? -1 : 1;
}
