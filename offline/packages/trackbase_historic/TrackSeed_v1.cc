#include "TrackSeed_v1.h"
#include <trackbase/TrkrCluster.h>

#include <trackbase/TrackFitUtils.h>

#include <cmath>

namespace
{

  //! convenience square method
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }
}  // namespace

TrackSeed_v1::TrackSeed_v1(const TrackSeed& seed)
{
  TrackSeed_v1::CopyFrom(seed);
}

// have to suppress missingMemberCopy from cppcheck, it does not
// go down to the CopyFrom method where things are done correctly
// cppcheck-suppress missingMemberCopy
TrackSeed_v1::TrackSeed_v1(const TrackSeed_v1& seed)
  : TrackSeed(seed)
{
  TrackSeed_v1::CopyFrom(seed);
}

TrackSeed_v1& TrackSeed_v1::operator=(const TrackSeed_v1& seed)
{
  if (this != &seed)
  {
    CopyFrom(seed);
  }
  return *this;
}

void TrackSeed_v1::CopyFrom(const TrackSeed& seed)
{
  if (this == &seed)
  {
    return;
  }
  TrackSeed::CopyFrom(seed);

  m_qOverR = seed.get_qOverR();
  m_X0 = seed.get_X0();
  m_Y0 = seed.get_Y0();
  m_slope = seed.get_slope();
  m_Z0 = seed.get_Z0();
  m_crossing = seed.get_crossing();

  m_cluster_keys.clear();
  std::copy(seed.begin_cluster_keys(), seed.end_cluster_keys(),
            std::inserter(m_cluster_keys, m_cluster_keys.begin()));
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
  os << "R and slope " << fabs(1. / m_qOverR) << ", " << m_slope << std::endl;
  os << "list of cluster keys size: " << m_cluster_keys.size() << std::endl;
  ;
  if (m_cluster_keys.size() > 0)
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

float TrackSeed_v1::get_pt() const
{
  /// Scaling factor for radius in 1.4T field
  return 0.3 * 1.4 / 100. * fabs(1. / m_qOverR);
}

float TrackSeed_v1::get_theta() const
{
  float theta = atan(1. / m_slope);
  /// Normalize to 0<theta<pi
  if (theta < 0)
  {
    theta += M_PI;
  }
  return theta;
}

float TrackSeed_v1::get_eta() const
{
  return -log(tan(get_theta() / 2.));
}

float TrackSeed_v1::get_p() const
{
  return get_pt() * std::cosh(get_eta());
}

float TrackSeed_v1::get_pz() const
{
  return get_p() * std::cos(get_theta());
}

int TrackSeed_v1::get_charge() const
{
  return (m_qOverR < 0) ? -1 : 1;
}
