#ifndef __PHG4DSTCOMPRESSRECO__
#define __PHG4DSTCOMPRESSRECO__

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4HitContainer.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4cemc/RawTowerContainer.h>

#include <set>
#include <string>

class PHG4DstCompressReco : public SubsysReco {
  
public:

  PHG4DstCompressReco(const std::string &name = "PHG4DstCompressReco");
  virtual ~PHG4DstCompressReco(){}
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
  //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode) {return Fun4AllReturnCodes::EVENT_OK;}
  
  void AddHitContainer(const std::string name) {_compress_g4hit_names.insert(name);}
  void AddCellContainer(const std::string name) {_compress_g4cell_names.insert(name);}
  void AddTowerContainer(const std::string name) {_compress_tower_names.insert(name);}

private:

  void SearchG4HitNodes(PHCompositeNode *topNode);
  
  PHG4TruthInfoContainer* _truth_info;
  std::set<std::string> _compress_g4hit_names;
  std::set<std::string> _compress_g4cell_names;
  std::set<std::string> _compress_tower_names;

  std::set<PHG4CylinderCellContainer*> _g4cells;
  std::set<PHG4HitContainer*> _g4hits;  
  std::set<PHG4HitContainer*> _keep_g4hits;

  std::set<RawTowerContainer*> _towers;
};

#endif
