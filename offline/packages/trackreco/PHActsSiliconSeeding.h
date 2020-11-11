#ifndef TRACKRECO_PHACTSSILICONSEEDING_H
#define TRACKRECO_PHACTSILICONSEEDING_H

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

using SpacePointPtr = std::unique_ptr<SpacePoint>;

class PHActsSiliconSeeding : public SubsysReco
{
 public:
  PHActsSiliconSeeding(std::string& name);
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
  float m_minSeedPt = 0.1 * Acts::UnitConstants::GeV;

  /// How many seeds a given hit can be the middle hit of the seed
  int m_maxSeedsPerSpM = 3;

  /// Limiting location of measurements (e.g. detector constraints)
  float m_rMax = 15. * Acts::UnitConstants::cm;
  float m_zMax = 20. * Acts::UnitConstants::cm;
  float m_zMin = -20. * Acts::UnitConstants::cm;
 
  /// max distance between two measurements in one seed
  float m_deltaRMax = 15. * Acts::UnitConstants::cm;
  
  /// Cot of maximum theta angle. Equivalent to eta=1.1 here
  float m_cotThetaMax = 1.335647;
  
  /// B field value in z direction
  float m_bField = 1.4 * Acts::UnitConstants::T;

  std::shared_ptr<Acts::BinFinder<SpacePoint>> m_bottomBinFinder, m_topBinFinder;

};


#endif
