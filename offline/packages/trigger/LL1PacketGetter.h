// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef LL1PACKETGETTER_H
#define LL1PACKETGETTER_H

#include <fun4all/SubsysReco.h>
#include "LL1Outv2.h"
#include "TriggerPrimitivev1.h"
#include "TriggerPrimitiveContainerv1.h"
#include "TriggerDefs.h"
#include <climits>
#include <string>


class PHCompositeNode;
class LL1Outv2;
class TriggerPrimitive;
class TriggerPrimitiveContainerv1;

class LL1PacketGetter : public SubsysReco
{
 public:

  explicit LL1PacketGetter(const std::string &name = "LL1PacketGetter", const std::string &trigger = "NONE", const std::string &ll1 = "NONE");
  ~LL1PacketGetter() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int CreateNodeTree(PHCompositeNode *topNode);

  void SetVerbosity(int v){_verbose = v;}
  void set_nsamples(int _nsamples)
  {
    m_nsamples = _nsamples;
    return;
  }

  void setTriggerType(const std::string &name);

  void set_dataflag(bool flag)
  {
    m_isdata = flag;
    return;
  }

 private:
  std::string m_trigger = "NONE";
  std::string m_ll1 = "NONE";

  LL1Outv2 *m_ll1out = nullptr;
  TriggerPrimitiveContainer *_trigger_primitives = nullptr;
  TriggerPrimitivev1 *_trigger_primitive = nullptr;
  std::map<unsigned int, std::vector<unsigned int>> *_trigger_words = nullptr;

  std::map<TriggerDefs::DetectorId, int> m_prim_map;
  TriggerDefs::TriggerId m_triggerid;
  TriggerDefs::DetectorId m_detectorid;
  TriggerDefs::PrimitiveId m_primitiveid;
  std::vector<unsigned int> *_sum = nullptr;

  int m_packet_low{0};
  int m_packet_high{0};

  int m_nsamples{16};
  int m_nchannels{256};
  int m_nchannels_per_primitive{0};
  int m_nprimitives{0};
  int m_ntriggerwords{0};
  int _verbose{0};
  bool m_isdata{0};
};

#endif  // LL1TOWERBUILDER_H
