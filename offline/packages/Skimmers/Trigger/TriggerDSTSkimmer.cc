
#include "TriggerDSTSkimmer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <ffarawobjects/Gl1Packet.h>

#include <phool/getClass.h>

#include <iostream>
#include <map>
#include <utility>

//____________________________________________________________________________..
TriggerDSTSkimmer::TriggerDSTSkimmer(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int TriggerDSTSkimmer::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    if (ievent % 1000 == 0)
    {
      std::cout << "Processing event " << ievent << std::endl;
    }
    ievent++;
  }

  // trigger cut
  int n_trigger_index = m_trigger_index.size();
  for (int it = 0; it < n_trigger_index; ++it)
  {
    if (m_trigger_index[it] < 0 || m_trigger_index[it] > 63)
    {
      std::cerr << "ERROR in TriggerDSTSkimmer: trigger index " 
                << m_trigger_index[it] << " is out of range [0..63]. "
                << "Aborting run." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  
  if (n_trigger_index != 0)
  {
    bool trigger_fired = false;
    Gl1Packet *_gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    int gl1_trigger_vector_scaled[64] = {0};
    if (_gl1PacketInfo)
    {
      uint64_t scaledtriggervec = _gl1PacketInfo->lValue(0, "ScaledVector");
      for (auto &val : gl1_trigger_vector_scaled)
      {
        val = ((scaledtriggervec & 0x1U) == 0x1U);
        scaledtriggervec = (scaledtriggervec >> 1U) & 0xffffffffU;
      }
    }
    else
    {
      std::cout << "TriggerDSTSkimmer::process_event - Error - Can't find Trigger Node Gl1Packet therefore no selection can be made" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    for (int it = 0; it < n_trigger_index; ++it)
    {
      if (gl1_trigger_vector_scaled[m_trigger_index[it]] == 1)
      {
        trigger_fired = true;
        break;
      }
    }
    if (!trigger_fired)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
