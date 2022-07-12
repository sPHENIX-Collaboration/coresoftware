#include "TrackSeed_v1.h"

#include <trackbase/TrackFitUtils.h>

namespace
{

  //! convenience square method
  template<class T>
    inline constexpr T square( const T& x ) { return x*x; }
}

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
  os << "(pt,pz) = (" << get_pt()
     << ", " << get_pz() << ")" << std::endl;
  os << "(x,y,z) = (" << get_x() << ", " << get_y() << ", " << get_z() 
     << ")" << std::endl;
  os << "(X0,Y0,Z0) = (" << m_X0 << ", " << m_Y0 << ", " << m_Z0
     << ")" << std::endl;
  os << "R and slope " << fabs(1./m_qOverR) << ", " << m_slope << std::endl;
  os << "list of cluster keys size: " << m_cluster_keys.size() << std::endl;;
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
				     ActsGeometry *tGeometry,
				     uint8_t startLayer,
				     uint8_t endLayer)
{
  std::map<TrkrDefs::cluskey, Acts::Vector3> positions;

  for(const auto& key: m_cluster_keys)
    {
      auto layer = TrkrDefs::getLayer(key);
      if(layer < startLayer or layer > endLayer)
	{ continue; }
      
      Acts::Vector3 pos = tGeometry->getGlobalPosition(
           key, clusters->findCluster(key));

      positions.insert(std::make_pair(key, pos));
    }

  circleFitByTaubin(positions, startLayer, endLayer);

}

void TrackSeed_v1::lineFit(TrkrClusterContainer *clusters,
			   ActsGeometry *tGeometry,
			   uint8_t startLayer,
			   uint8_t endLayer)
{
  std::map<TrkrDefs::cluskey, Acts::Vector3> positions;

  for(const auto& key: m_cluster_keys)
    {
      auto layer = TrkrDefs::getLayer(key);
      if(layer < startLayer or layer > endLayer)
	{ continue; }
      
      Acts::Vector3 pos = tGeometry->getGlobalPosition(
           key, clusters->findCluster(key));

      positions.insert(std::make_pair(key, pos));
    }

  lineFit(positions, startLayer, endLayer);

}

void TrackSeed_v1::circleFitByTaubin(std::map<TrkrDefs::cluskey, Acts::Vector3>& positions,
				       uint8_t startLayer,
				       uint8_t endLayer)
{

  TrackFitUtils::position_vector_t positions_2d;
  for(const auto& key: m_cluster_keys)
  {
    const auto layer = TrkrDefs::getLayer(key);
    if(layer < startLayer or layer > endLayer) { continue; }

    const auto iter = positions.find(key);
    
    /// you supplied the wrong key...
    if(iter == positions.end())
    { continue; }

    // add to 2d position list
    const Acts::Vector3& pos = iter->second;
    positions_2d.emplace_back( pos.x(), pos.y() );
  }
  
  // do the fit
  const auto [r, x0, y0] = TrackFitUtils::circle_fit_by_taubin( positions_2d );
  
  // assign 
  m_X0 = x0;
  m_Y0 = y0;
  m_qOverR = 1. / r;
 
  /// Set the charge
  const auto& firstpos = positions_2d.at(0);
  const auto& secondpos = positions_2d.at(1);
  
  const auto firstphi = atan2(firstpos.second, firstpos.first);
  const auto secondphi = atan2(secondpos.second, secondpos.first);
  auto dphi = secondphi - firstphi;
  if(dphi > M_PI) dphi = 2.*M_PI - dphi;
  if(dphi < -M_PI) dphi = 2*M_PI + dphi;
  if(dphi > 0) m_qOverR *= -1;

}

void TrackSeed_v1::lineFit(std::map<TrkrDefs::cluskey, Acts::Vector3>& positions,
			   uint8_t startLayer,
			   uint8_t endLayer)
{
  
  TrackFitUtils::position_vector_t positions_2d;
  for(const auto& key : m_cluster_keys)
  {
    const auto layer = TrkrDefs::getLayer(key);
    if(layer < startLayer or layer > endLayer)
    { continue; }

    const auto iter = positions.find(key);

    /// The wrong key was supplied...
    if(iter == positions.end())
    { continue; }

    // store (r,z)
    const Acts::Vector3& pos = iter->second;
    positions_2d.emplace_back( std::sqrt( square(pos.x()) + square(pos.y())), pos.z() );
  }
  
  // do the fit
  const auto [slope, intercept] = TrackFitUtils::line_fit( positions_2d );
  
  // assign
  m_slope = slope;
  m_Z0 = intercept;
  
}   

std::pair<float,float> TrackSeed_v1::findRoot() const
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
  const float R = std::abs(1./m_qOverR);
  const double miny = (std::sqrt(square(m_X0) * square(R) * square(m_Y0) + square(R) 
		      * pow(m_Y0,4)) + square(m_X0) * m_Y0 + pow(m_Y0, 3)) 
    / (square(m_X0) + square(m_Y0));

  const double miny2 = (-std::sqrt(square(m_X0) * square(R) * square(m_Y0) + square(R) 
		      * pow(m_Y0,4)) + square(m_X0) * m_Y0 + pow(m_Y0, 3)) 
    / (square(m_X0) + square(m_Y0));

  const double minx = std::sqrt(square(R) - square(miny - m_Y0)) + m_X0;
  const double minx2 = -std::sqrt(square(R) - square(miny2 - m_Y0)) + m_X0;
  
  /// Figure out which of the two roots is actually closer to the origin
  const float x = ( std::abs(minx) < std::abs(minx2)) ? minx:minx2;
  const float y = ( std::abs(miny) < std::abs(miny2)) ? miny:miny2;
  return std::make_pair(x,y);
  
}

float TrackSeed_v1::get_x() const
{ return findRoot().first; }

float TrackSeed_v1::get_y() const
{ return findRoot().second; }

float TrackSeed_v1::get_z() const
{ return get_Z0(); }

float TrackSeed_v1::get_pt() const
{
  /// Scaling factor for radius in 1.4T field
  return 0.3 * 1.4 / 100. * fabs(1./m_qOverR);
}
float TrackSeed_v1::get_phi(std::map<TrkrDefs::cluskey, Acts::Vector3>& positions) const
{
  const auto [x,y] = findRoot();
  float phi = atan2(-1* (m_X0-x), (m_Y0-y));
  Acts::Vector3 pos0 = positions.find(*(m_cluster_keys.begin()))->second;
  Acts::Vector3 pos1 = positions.find(*(std::next(m_cluster_keys.begin(), 1)))->second;
  /// convert to the angle of the tangent to the circle
  // we need to know if the track proceeds clockwise or CCW around the circle
  double dx0 = pos0(0) - m_X0;
  double dy0 = pos0(1) - m_Y0;
  double phi0 = atan2(dy0, dx0);
  double dx1 = pos1(0) - m_X0;
  double dy1 = pos1(1) - m_Y0;
  double phi1 = atan2(dy1, dx1);
  double dphi = phi1 - phi0;

  // need to deal with the switch from -pi to +pi at phi = 180 degrees
  // final phi - initial phi must be < 180 degrees for it to be a valid track
  if(dphi > M_PI) dphi -= 2.0 * M_PI;
  if(dphi < - M_PI) dphi += M_PI;

  // whether we add 180 degrees depends on the angle of the bend
  if(dphi < 0)
    { 
      phi += M_PI; 
      if(phi > M_PI)
	{ phi -= 2. * M_PI; }
    }
  
  return phi;
}
float TrackSeed_v1::get_phi(TrkrClusterContainer *clusters,
			    ActsGeometry *tGeometry) const
{
  Acts::Vector3 pos0 = tGeometry->getGlobalPosition(
		       *(m_cluster_keys.begin()),
		       clusters->findCluster(*(m_cluster_keys.begin())));

  auto key = *std::next(m_cluster_keys.begin(), 1);
  Acts::Vector3 pos1 = tGeometry->getGlobalPosition(key,
					     clusters->findCluster(key));
  std::map<TrkrDefs::cluskey, Acts::Vector3> positions;
  positions.insert(std::make_pair(*(m_cluster_keys.begin()), pos0));
  positions.insert(std::make_pair(key, pos1));
  return get_phi(positions);
 
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

float TrackSeed_v1::get_px(TrkrClusterContainer *clusters,
			   ActsGeometry *tGeometry) const
{
  return get_pt() * cos(get_phi(clusters, tGeometry));
}

float TrackSeed_v1::get_py(TrkrClusterContainer *clusters,
			   ActsGeometry *tGeometry) const
{
  return get_pt() * sin(get_phi(clusters, tGeometry));
}

float TrackSeed_v1::get_pz() const
{
  return get_p() * cos(get_theta());
}

int TrackSeed_v1::get_charge() const
{
  return ( m_qOverR < 0 ) ? -1 : 1;
}
