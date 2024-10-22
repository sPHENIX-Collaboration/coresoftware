#include "TrackSeed_v2.h"
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

TrackSeed_v2::TrackSeed_v2(const TrackSeed& seed)
{
  TrackSeed_v2::CopyFrom(seed);
}

// have to suppress missingMemberCopy from cppcheck, it does not
// go down to the CopyFrom method where things are done correctly
// cppcheck-suppress missingMemberCopy
TrackSeed_v2::TrackSeed_v2(const TrackSeed_v2& seed)
  : TrackSeed(seed)
{
  TrackSeed_v2::CopyFrom(seed);
}

TrackSeed_v2& TrackSeed_v2::operator=(const TrackSeed_v2& seed)
{
  if (this != &seed)
  {
    CopyFrom(seed);
  }
  return *this;
}

void TrackSeed_v2::CopyFrom(const TrackSeed& seed)
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
  m_phi = seed.get_phi();
  m_cluster_keys.clear();
  std::copy(seed.begin_cluster_keys(), seed.end_cluster_keys(),
            std::inserter(m_cluster_keys, m_cluster_keys.begin()));
}

void TrackSeed_v2::identify(std::ostream& os) const
{
  os << "TrackSeed_v2 object ";
  os << "charge " << get_charge() << std::endl;
  os << "beam crossing " << get_crossing() << std::endl;
  os << "(pt,pz) = (" << get_pt()
     << ", " << get_pz() << ")" << std::endl;
  os << " phi " << m_phi << " eta " << get_eta() << std::endl;
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

float TrackSeed_v2::get_pt() const
{
  /// Scaling factor for radius in 1.4T field
  return 0.3 * 1.4 / 100. * fabs(1. / m_qOverR);
}

float TrackSeed_v2::get_theta() const
{
  float theta = atan(1. / m_slope);
  /// Normalize to 0<theta<pi
  if (theta < 0)
  {
    theta += M_PI;
  }
  return theta;
}

float TrackSeed_v2::get_eta() const
{
  return -log(tan(get_theta() / 2.));
}

float TrackSeed_v2::get_p() const
{
  return get_pt() * std::cosh(get_eta());
}

float TrackSeed_v2::get_px() const
{
  return get_pt() * std::cos(m_phi);
}

float TrackSeed_v2::get_py() const
{
  return get_pt() * std::sin(m_phi);
}

float TrackSeed_v2::get_pz() const
{
  return get_p() * std::cos(get_theta());
}

int TrackSeed_v2::get_charge() const
{
  return (m_qOverR < 0) ? -1 : 1;
}
