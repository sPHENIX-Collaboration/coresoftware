#ifndef TPC_POLYCLUSTERRESIDUALS_H
#define TPC_POLYCLUSTERRESIDUALS_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class Tpc_PolyTrackContainer;
class Tpc_PolyTrackVertexContainer;
class PHCompositeNode;
class TFile;
class TTree;
class Tpc_PolyClusterContainer;

class Tpc_PolyClusterResiduals : public SubsysReco
{
 public:
  explicit Tpc_PolyClusterResiduals(const std::string& name = "Tpc_PolyClusterResiduals",
                                   const std::string& outfilename = "tpc_polycluster_residuals.root");
  ~Tpc_PolyClusterResiduals() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setClusterNodeName(const std::string& n) { m_clusterNodeName = n; }
  void setTpc_PolyTrackNodeName(const std::string& n) { m_finalTrackNodeName = n; }
  void setTpc_PolyTrackVertexNodeName(const std::string& n) { m_finalTrackVertexNodeName = n; }
  void setOutputFileName(const std::string& n) { m_outfilename = n; }
  void setMagneticFieldTesla(double b) { m_magneticFieldTesla = b; }
  void setMinPt(double v) { m_minPt = v; }
  void setMaxPt(double v) { m_maxPt = v; }
  void setMinTpcClusters(unsigned int v) { m_minTpcClusters = v; }
  void setMaxTpcClusters(unsigned int v) { m_maxTpcClusters = v; }
  void setUseStraightLineTracks(bool v) { m_useStraightLineTracks = v; }

 private:
  bool get_nodes(PHCompositeNode* topNode);
  void reset_tree_values();

  std::string m_outfilename;
  std::string m_clusterNodeName;
  std::string m_finalTrackNodeName;
  std::string m_finalTrackVertexNodeName;

  double m_magneticFieldTesla {1.4};
  double m_minPt {0.0};
  double m_maxPt {1.0e30};
  unsigned int m_minTpcClusters {0};
  unsigned int m_maxTpcClusters {0xffffffffu};
  bool m_useStraightLineTracks {false};

  unsigned int m_evt {0};
  TFile* m_outfile {nullptr};
  TTree* m_tree {nullptr};
  Tpc_PolyClusterContainer* m_clusters {nullptr};
  Tpc_PolyTrackContainer* m_finalTracks {nullptr};
  Tpc_PolyTrackVertexContainer* m_finalTrackVertices {nullptr};

  unsigned int m_event {0};
  unsigned int m_finalTrackId {0};
  unsigned int m_sourceClusterId {0};
  unsigned int m_sourceAssembledTrackId {0};
  int m_side {0};
  unsigned int m_ntpcClusters {0};
  int m_fitStatus {0};
  double m_pt {0.0};
  double m_px {0.0};
  double m_py {0.0};
  double m_pz {0.0};
  double m_eta {0.0};
  double m_theta {0.0};
  double m_charge {0.0};
  double m_chi2 {0.0};
  double m_ndf {0.0};
  double m_quality {0.0};
  double m_dedx {0.0};
  double m_vertexX {0.0};
  double m_vertexY {0.0};
  double m_vertexZ {0.0};
  double m_vertexR {0.0};
  double m_collisionXSide0 {0.0};
  double m_collisionYSide0 {0.0};
  double m_collisionZSide0 {0.0};
  double m_collisionXSide1 {0.0};
  double m_collisionYSide1 {0.0};
  double m_collisionZSide1 {0.0};
  double m_collisionXMid {0.0};
  double m_collisionYMid {0.0};
  double m_collisionZMid {0.0};
  double m_pcaX {0.0};
  double m_pcaY {0.0};
  double m_pcaZ {0.0};
  double m_zDCA {0.0};
  double m_zDCAMid {0.0};
  double m_rDCA {0.0};
  double m_rDCAZero {0.0};
  double m_R {0.0};
  double m_rzSlope {0.0};
  std::vector<unsigned int> m_clusterIndex;
  std::vector<unsigned int> m_sector;
  std::vector<unsigned int> m_layer;
  std::vector<double> m_clusterX;
  std::vector<double> m_clusterY;
  std::vector<double> m_clusterZ;
  std::vector<double> m_clusterR;
  std::vector<double> m_clusterPhi;
  std::vector<double> m_clusterAdc;
  std::vector<unsigned int> m_clusterPadSize;
  std::vector<double> m_stateX;
  std::vector<double> m_stateY;
  std::vector<double> m_stateZ;
  std::vector<double> m_stateZDca;
  std::vector<double> m_stateR;
  std::vector<double> m_statePhi;
  std::vector<double> m_deltaPhi;
  std::vector<double> m_residualRPhi;
  std::vector<double> m_residualZ;
};

#endif
