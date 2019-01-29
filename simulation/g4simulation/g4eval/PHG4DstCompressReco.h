#ifndef G4EVAL_PHG4DSTCOMPRESSRECO_H
#define G4EVAL_PHG4DSTCOMPRESSRECO_H

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class PHCompositeNode;

class PHG4CellContainer;
class PHG4HitContainer;
class PHG4TruthInfoContainer;

class RawTowerContainer;

class PHG4DstCompressReco : public SubsysReco {
  
public:

  PHG4DstCompressReco(const std::string &name = "PHG4DstCompressReco");
  virtual ~PHG4DstCompressReco(){}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
  //! event processing
  int process_event(PHCompositeNode *topNode);
  
  void AddHitContainer(const std::string &name) {_compress_g4hit_names.insert(name);}
  void AddCellContainer(const std::string &name) {_compress_g4cell_names.insert(name);}
  void AddTowerContainer(const std::string &name) {_compress_tower_names.insert(name);}

private:

  void SearchG4HitNodes(PHCompositeNode *topNode);
  
  PHG4TruthInfoContainer* _truth_info;
  std::set<std::string> _compress_g4hit_names;
  std::set<std::string> _compress_g4cell_names;
  std::set<std::string> _compress_tower_names;

  std::set<PHG4CellContainer*> _g4cells;
  std::set<PHG4HitContainer*> _g4hits;  
  std::set<PHG4HitContainer*> _keep_g4hits;

  std::set<RawTowerContainer*> _towers;
};

#endif
