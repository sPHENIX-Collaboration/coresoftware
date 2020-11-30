#ifndef TRACKRECO_PHACTSSILICONSEEDING_H
#define TRACKRECO_PHACTSSILICONSEEDING_H

#include "ActsTrackingGeometry.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Seeding/Seedfinder.hpp>
#include <Acts/Utilities/Units.hpp>

#include <Acts/Seeding/BinFinder.hpp>
#include <Acts/Seeding/SpacePointGrid.hpp>

#include <ActsExamples/EventData/TrkrClusterSourceLink.hpp>

#include <string>
#include <map>
class PHCompositeNode;

using SourceLink = ActsExamples::TrkrClusterSourceLink;

/**
 * A struct for Acts to take cluster information for seeding
 */
struct SpacePoint {
  unsigned int m_hitId;
  float m_x;
  float m_y;
  float m_z;
  float m_r;
  Acts::GeometryIdentifier m_geoId;
  float m_varianceRphi;
  float m_varianceZ;
  
  unsigned int Id() const { return m_hitId; }

  /// These are needed by Acts
  float x() const { return m_x; }
  float y() const { return m_y; }
  float z() const { return m_z; }
  float r() const { return m_r; }

};

inline bool operator==(SpacePoint a, SpacePoint b) {
  return (a.m_hitId == b.m_hitId);
}

using SpacePointPtr = std::unique_ptr<SpacePoint>;

class PHActsSiliconSeeding : public SubsysReco
{
 public:
  PHActsSiliconSeeding(const std::string& name = "PHActsSiliconSeeding");
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:

  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);

  Acts::SeedfinderConfig<SpacePoint> configureSeeder();
  Acts::SpacePointGridConfig configureSPGrid();
  
  SpacePointPtr makeSpacePoint(const unsigned int& hitId,
			    const SourceLink& sl);

  std::map<unsigned int, SourceLink> *m_sourceLinks;
  ActsTrackingGeometry *m_tGeometry;
  
  Acts::SeedfinderConfig<SpacePoint> m_seedFinderCfg;
  Acts::SpacePointGridConfig m_gridCfg;

  /// Configurable parameters
  /// seed pt has to be in MeV
  float m_minSeedPt = 100;

  /// How many seeds a given hit can be the middle hit of the seed
  int m_maxSeedsPerSpM = 1;

  /// Limiting location of measurements (e.g. detector constraints)
  float m_rMax = 250.;
  float m_rMin = 20.;
  float m_zMax = 500.;
  float m_zMin = -500.;
 
  /// max distance between two measurements in one seed
  float m_deltaRMax = 50;
  
  /// Cot of maximum theta angle. Equivalent to eta=1.1 here
  float m_cotThetaMax = 1.335647;
  
  /// B field value in z direction
  /// bfield for space point grid neds to be in kiloTesla
  float m_bField = 1.4 / 1000.;

  std::shared_ptr<Acts::BinFinder<SpacePoint>> 
    m_bottomBinFinder, m_topBinFinder;

  int m_event = 0;

};


#endif
