#ifndef BCOLUMICOUNT_BCOLUMIRECO_H
#define BCOLUMICOUNT_BCOLUMIRECO_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>

class SyncObject;

class BcoLumiReco : public SubsysReco
{
 public:
  BcoLumiReco(const std::string &name = "BCOLUMIRECO");
  ~BcoLumiReco() override = default;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void push(uint64_t value);
  uint64_t get_previous_bco() {return  bco[0];}
  uint64_t get_current_bco() const {return bco[1];}
  uint64_t get_future_bco() const {return bco[2];}
 private:
  static int CreateNodeTree(PHCompositeNode *topNode);
  SyncObject *synccopy {nullptr};
  SyncObject *tmpsync {nullptr};
  std::array<uint64_t,3> bco {0};
};

#endif // BCOLUMICOUNT_BCOLUMIRECO_H
