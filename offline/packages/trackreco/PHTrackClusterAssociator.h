// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHACTSTRACKCLUSTERASSOCIATOR_H
#define PHACTSTRACKCLUSTERASSOCIATOR_H

#include <fun4all/SubsysReco.h>
#include <trackbase_historic/SvtxTrack.h>

#include <string>
#include <map>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class RawCluster;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;

class PHTrackClusterAssociator : public SubsysReco
{
 public:
  typedef std::map<SvtxTrack*, RawCluster*> TrackClusterMap;

  PHTrackClusterAssociator(const std::string &name = "PHTrackClusterAssociator");

  ~PHTrackClusterAssociator() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:

  int getNodes(PHCompositeNode* topNode);
  int createNodes(PHCompositeNode* topNode);
  int getCaloNodes(PHCompositeNode* topNode, const int caloLayer);
  int matchTracks(PHCompositeNode* topNode, const int caloLayer);

  SvtxTrackMap* m_trackMap = nullptr;

  /// Objects to hold calorimeter information. There are 
  /// only 3 calo layers
  const static int m_nCaloLayers = 3;
  std::vector<std::string> m_caloNames;
  std::vector<SvtxTrack::CAL_LAYER> m_caloTypes;
  /// An optional map that allows projection to an arbitrary radius
  /// Results are written to the SvtxTrack based on the provided CAL_LAYER
  std::map<SvtxTrack::CAL_LAYER, float> m_caloRadii;

  RawTowerGeomContainer *m_towerGeomContainer = nullptr;
  RawTowerContainer *m_towerContainer = nullptr;
  RawClusterContainer *m_clusterContainer = nullptr;

  TrackClusterMap* m_trackClusterMap = nullptr;
  
  bool m_useCemcPosRecalib = false;
  bool m_calosAvailable = true;
  
};

#endif // PHACTSTRACKCLUSTERASSOCIATOR_H
