#ifndef TPC_TPCMISALIGNER_HP_H
#define G4EVAL_TRACKINGEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <trackbase/TrkrDefs.h>

class TrkrCluster;
class TrkrClusterContainer;

class TpcMisaligner_hp : public SubsysReco
{
  public:

  /// constructor
  TpcMisaligner_hp(  const std::string& = "TPCMISALIGNER_HP" );

  /// global initialization
  virtual int Init(PHCompositeNode*);

  /// event processing
  virtual int process_event(PHCompositeNode*);

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// transform clusters
  void transform_clusters();

  /// transform clusters
  void transform_cluster( TrkrCluster* );

  /// cluster container
  TrkrClusterContainer* _cluster_map = nullptr;

};

#endif
