#include <Acts/EventData/Measurement.hpp>
#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/EventData/SourceLinkConcept.hpp>
#include <Acts/EventData/detail/fittable_type_generator.hpp>


/**
 * This class creates an Acts::SourceLink that relates TrkrClusters to the
 * surface they were measured on. The source link is needed for the fitting
 */
class TrkrClusterSourceLink
{
public:

  /// Instantiate with a clusid, associated surface, and values that actually
  /// make the measurement. Acts requires the surface be available in this class
  TrkrClusterSourceLink(unsigned int clusid,
			std::shared_ptr<const Acts::Surface> surface,
			Acts::BoundVector loc,
			Acts::ActsSymMatrixD<3> cov)
    : m_clusid(clusid)
    , m_surface(surface)
    , m_loc(loc)
    , m_cov(cov)
{
}

  /// Must be default constructible to satisfy SourceLinkConcept
  TrkrClusterSourceLink()    = default;
  TrkrClusterSourceLink(TrkrClusterSourceLink&&)      = default;
  TrkrClusterSourceLink(const TrkrClusterSourceLink&) = default;

  /// Needs equality operators defined to satisfy SourceLinkConcept
  TrkrClusterSourceLink& operator=(TrkrClusterSourceLink&&)      = default;
  TrkrClusterSourceLink& operator=(const TrkrClusterSourceLink&) = default;
  
  /// Needs referenceSurface function to satisfy SourceLinkConcept
  const Acts::Surface& referenceSurface() const 
  {
    return *m_surface;
  }
  
  /// Create Acts::FittableMeasurement from information in SourceLink
  Acts::FittableMeasurement<TrkrClusterSourceLink> operator*() const
  {

    return Acts::Measurement<TrkrClusterSourceLink, 
			     Acts::ParDef::eLOC_0,
			     Acts::ParDef::eLOC_1>
      {
	m_surface,
	  *this,
	  m_cov.topLeftCorner<2, 2>(),
	  m_loc[0],
	  m_loc[1]};
  }

private:
  /// Hitindex corresponding to clusID and the corresponding 
  /// surface to which it belongs to
  unsigned int m_clusid;
  std::shared_ptr<const Acts::Surface> m_surface;
  /// Local x and y position for cluster
  Acts::BoundVector m_loc;
  /// Cluster covariance matrix
  Acts::ActsSymMatrixD<3> m_cov;

  /// Needs equality operator defined to satisfy SourceLinkConcept
  /// Equate the cluster keys
  friend constexpr bool
  operator==(const TrkrClusterSourceLink& lhs, const TrkrClusterSourceLink& rhs)
  {
    return lhs.m_clusid == rhs.m_clusid;
  }

};


/// Ensure that the SourceLink class satisfies SourceLinkConcept conditions
static_assert(Acts::SourceLinkConcept<TrkrClusterSourceLink>, 
	      "TrkrClusterSourceLink does not fulfill SourceLinkConcept");
