#ifndef TRACKRECO_PHACTSTRACKPROJECTION_H
#define TRACKRECO_PHACTSTRACKPROJECTION_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>

#include <trackbase/ActsGeometry.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/TrackParameters.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Acts/Propagator/Propagator.hpp>
#pragma GCC diagnostic pop

#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Utilities/Result.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>

class PHCompositeNode;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;

#include <map>
#include <memory>
#include <string>

using BoundTrackParam =
    const Acts::BoundTrackParameters;
using BoundTrackParamResult = Acts::Result<BoundTrackParam>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;
using Trajectory = ActsExamples::Trajectories;

/**
 * This class takes final fitted tracks from the Acts track fitting
 * and projects them out to cylinders with radius at the same radius
 * as the three calorimeters. Cluster matching is performed with the
 * projections and the SvtxTrack object is updated.
 */

class PHActsTrackProjection : public SubsysReco
{
 public:
  PHActsTrackProjection(const std::string &name = "PHActsTrackProjection");

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void useConstField(bool field) { m_constField = field; }

  /// Set an arbitrary radius to project to, in cm
  void setLayerRadius(SvtxTrack::CAL_LAYER layer,
                      const float rad)
  {
    if (m_caloRadii.find(layer) != m_caloRadii.end())
    {
      m_caloRadii[layer] = rad;
    }
    else
    {
      m_caloRadii.insert(std::make_pair(layer, rad));
    }
  }

 private:
  int getNodes(PHCompositeNode *topNode);
  int projectTracks(int caloLayer);

  /// Propagate the fitted track parameters to a surface with Acts
  BoundTrackParamResult propagateTrack(
      const Acts::BoundTrackParameters &params,
      const SurfacePtr &targetSurf);

  /// Set the particular calo nodes depending on which layer
  int setCaloContainerNodes(PHCompositeNode *topNode,
                            const int caloLayer);

  /// Make Acts::CylinderSurface objects corresponding to the calos
  int makeCaloSurfacePtrs(PHCompositeNode *topNode);

  /// Update the SvtxTrack object with the track-cluster match
  void updateSvtxTrack(const Acts::BoundTrackParameters &params,
                       SvtxTrack *svtxTrack,
                       const int caloLayer);

  /// Get 3x3 and 5x5 tower sums matched to a track
  void getSquareTowerEnergies(int phiBin, int etaBin,
                              double &energy3x3,
                              double &energy5x5);

  /// Get the cluster values for a particular matched track
  void getClusterProperties(double phi, double eta,
                            double &minIndex, double &minDphi,
                            double &minDeta, double &minE);
  Acts::BoundTrackParameters makeTrackParams(SvtxTrack *track);
  double deltaPhi(const double &phi);
  Acts::Vector3 getVertex(SvtxTrack *track);

  /// Objects containing the Acts track fit results
  ActsGeometry *m_tGeometry = nullptr;
  SvtxTrackMap *m_trackMap = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;

  /// Objects to hold calorimeter information. There are
  /// only 3 calo layers
  const static int m_nCaloLayers = 3;
  std::vector<std::string> m_caloNames;
  std::vector<SvtxTrack::CAL_LAYER> m_caloTypes;
  std::map<std::string, SurfacePtr> m_caloSurfaces;
  /// An optional map that allows projection to an arbitrary radius
  /// Results are written to the SvtxTrack based on the provided CAL_LAYER
  std::map<SvtxTrack::CAL_LAYER, float> m_caloRadii;

  RawTowerGeomContainer *m_towerGeomContainer = nullptr;
  RawTowerContainer *m_towerContainer = nullptr;
  RawClusterContainer *m_clusterContainer = nullptr;

  bool m_constField = true;
  bool m_useCemcPosRecalib = false;
  bool m_calosAvailable = true;

  int m_event = 0;
};

#endif
