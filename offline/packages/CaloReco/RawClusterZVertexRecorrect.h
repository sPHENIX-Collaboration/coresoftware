#ifndef CALORECO_RAWCLUSTERZVTXRECORRECT_H
#define CALORECO_RAWCLUSTERZVTXRECORRECT_H

#include "CaloRecoUtility.h"

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
//class RawClusterContainer;

class RawClusterZVertexRecorrect : public SubsysReco
{
 public:
  explicit RawClusterZVertexRecorrect(const std::string &name);

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  //  void CreateNodeTree(PHCompositeNode *topNode);


  void set_UseTowerInfo(const int useMode)
  {  // 0 only old tower, 1 only new (TowerInfo based),
    m_UseTowerInfo = useMode;
  }


  void set_UseBbcZVtx(const bool useBbc)
  { 
    // this should (could?) be replaced by flags to GlobalVertexMap
    m_UseBbcZVtx = useBbc;
  }


 private:

  //  RawClusterContainer *_recalib_clusters{};

  std::string _det_name;

  CaloRecoUtility m_calrecoUtilInstance;
  int m_UseTowerInfo = 0;  // 0 only old tower, 1 only new (TowerInfo based),
  bool m_UseBbcZVtx = false;

};

#endif
