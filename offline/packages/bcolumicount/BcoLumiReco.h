#ifndef BCOLUMICOUNT_BCOLUMIRECO_H
#define BCOLUMICOUNT_BCOLUMIRECO_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <cstdint>
#include <string>

class PHCompositeNode;
class SyncObject;

class BcoLumiReco : public SubsysReco
{
 public:
  BcoLumiReco(const std::string &name = "BCOLUMIRECO");
  ~BcoLumiReco() override;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void push_bco(uint64_t value);
  uint64_t get_previous_bco() { return m_bco[0]; }
  uint64_t get_current_bco() const { return m_bco[1]; }
  uint64_t get_future_bco() const { return m_bco[2]; }

  void push_evtno(int value);
  int get_previous_evtno() { return m_evtno[0]; }
  int get_current_evtno() const { return m_evtno[1]; }
  int get_future_evtno() const { return m_evtno[2]; }

 private:
  static int CreateNodeTree(PHCompositeNode *topNode);
  SyncObject *m_synccopy{nullptr};
  SyncObject *m_tmpsync{nullptr};
  std::array<uint64_t, 3> m_bco{0};
  std::array<int, 3> m_evtno{0};
};

#endif  // BCOLUMICOUNT_BCOLUMIRECO_H
