#ifndef TPCCONDITIONSRECO__H
#define TPCCONDITIONSRECO__H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <map>
#include <string>
#include <vector>

class PHCompositeNode;

class CDBInterface;
class CDBTTree;
class TpcConditions;

class TpcConditionsReco : public SubsysReco
{
 public:
  TpcConditionsReco(const std::string &name = "TpcConditionsReco");
  ~TpcConditionsReco() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

 private:
  CDBInterface *m_cdb{nullptr};
  CDBTTree *m_tree{nullptr};
  TpcConditions *m_conditions{nullptr};
  float get_MedianCurrent(int channel, const std::vector<std::string> &channels);

  std::map<uint64_t, int> m_bco_to_channel;
};

#endif
