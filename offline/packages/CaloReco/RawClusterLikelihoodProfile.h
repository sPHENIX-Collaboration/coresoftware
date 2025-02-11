#ifndef RAWCLUSTERLIKELIHOODPROFILE_H
#define RAWCLUSTERLIKELIHOODPROFILE_H

#include "ClusterCDFCalculator.h"

#include <fun4all/SubsysReco.h>

class PHCompositeNode;
class RawClusterContainer;

class RawClusterLikelihoodProfile : public SubsysReco
{
 public:
  RawClusterLikelihoodProfile(const std::string &name = "RawClusterLikelihoodProfile");

  ~RawClusterLikelihoodProfile() override = default;

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;

  void set_profile_filepath(const std::string &filepath) {m_profile_path = filepath;}

  void set_inputNodeName(const std::string &inputNodeName) { m_inputNodeName = inputNodeName; }

  void set_outputNodeName(const std::string &outputNodeName) { m_outputNodeName = outputNodeName; }

  void set_towerNodeName(const std::string &towerNodeName) { m_towerNodeName = towerNodeName; }
  
  void set_globalVertexNodeName(const std::string &globalvertexNodeName) { m_globalvtxNodeName = globalvertexNodeName; }

  void set_min_cluster_e(const float min_cluster_e) { m_min_cluster_e = min_cluster_e; }
  
  void set_tower_threshold(const float tower_thres_e) { m_tower_thres_e = tower_thres_e; }
  
  void set_profile_dimension(const int nProfDim) { m_profile_dimension = nProfDim; }

 private:
  ClusterCDFCalculator *cdfcalc{nullptr};
  
  const int inputDimx{7};
  const int inputDimy{7};

  std::string m_profile_path;

  std::string m_inputNodeName{"CLUSTERINFO_CEMC"};
  std::string m_outputNodeName{"CLUSTERINFO_CEMC_PROFILE"};
  std::string m_towerNodeName{"TOWERINFO_CALIB_CEMC"};
  std::string m_globalvtxNodeName{"GlobalVertexMap"};

  bool inplace{false};

  RawClusterContainer *_clusters{nullptr};

  float m_min_cluster_e{0.5};
  float m_tower_thres_e{0.070};
  int m_profile_dimension{3};

  float vx{0};
  float vy{0};
  float vz{0};

  void CreateNodes(PHCompositeNode* topNode);

};

#endif
