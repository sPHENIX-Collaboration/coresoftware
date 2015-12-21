#ifndef __PHG4DSTCOMPRESSRECO__
#define __PHG4DSTCOMPRESSRECO__

#include <fun4all/SubsysReco.h>

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
  int End(PHCompositeNode *topNode) {}
  
  void AddHitContainer(std::string name) {_compress_g4hit_names.insert(name);}
  void AddCellContainer(std::string name) {_compress_g4cell_names.insert(name);}
  void AddTowerContainer(std::string name) {_compress_tower_names.insert(name);}

private:

  void SearchG4HitNodes(PHCompositeNode *topNode);
  
  PHG4TruthInfoContainer *_truth_info;
  std::set<std::string> _compress_g4hit_names;
  std::set<std::string> _compress_g4cell_names;
  std::set<std::string> _compress_tower_names;

  std::set<PHG4CylinderCellContainer*> _g4cells;
  std::set<PHG4HitContainer*> _g4hits;  
  std::set<PHG4HitContainer*> _keep_g4hits;

  std::set<RawTowerContainer*> _towers;
};

#endif
