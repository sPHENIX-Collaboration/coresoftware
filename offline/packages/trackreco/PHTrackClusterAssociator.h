// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHACTSTRACKCLUSTERASSOCIATOR_H
#define PHACTSTRACKCLUSTERASSOCIATOR_H

#include <fun4all/SubsysReco.h>
#include <trackbase_historic/SvtxTrack.h>

#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxVertexMap;
class SvtxVertex;
class RawCluster;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class SvtxTrackCaloClusterMap;

class PHTrackClusterAssociator : public SubsysReco
{
 public:
  PHTrackClusterAssociator(const std::string& name = "PHTrackClusterAssociator");

  ~PHTrackClusterAssociator() override;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int ResetEvent(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void setCaloRadius(const std::string& name, double rad)
  {
    m_caloRadii.insert(std::make_pair(name, rad));
  }

 private:
  int getNodes(PHCompositeNode* topNode);
  int createNodes(PHCompositeNode* topNode);
  int getCaloNodes(PHCompositeNode* topNode, const int caloLayer);
  int matchTracks(PHCompositeNode* topNode, const int caloLayer);
  RawCluster* getCluster(double phi, double eta, SvtxVertex* vertex);
  SvtxTrackMap* m_trackMap = nullptr;
  SvtxVertexMap* m_vertexMap = nullptr;

  /// Objects to hold calorimeter information. There are
  /// only 3 calo layers
  const static int m_nCaloLayers = 3;
  std::vector<std::string> m_caloNames;
  /// An optional map that allows projection to an arbitrary radius
  /// Results are written to the SvtxTrack based on the provided CAL_LAYER
  std::map<std::string, float> m_caloRadii;

  RawTowerGeomContainer* m_towerGeomContainer = nullptr;
  RawTowerContainer* m_towerContainer = nullptr;
  RawClusterContainer* m_clusterContainer = nullptr;

  SvtxTrackCaloClusterMap* m_trackClusterMap = nullptr;

  bool m_useCemcPosRecalib = false;
  bool m_calosAvailable = true;
};

#endif  // PHACTSTRACKCLUSTERASSOCIATOR_H
