#ifndef INTT_COMBINEDRAWDATADECODER_H
#define INTT_COMBINEDRAWDATADECODER_H

#include "InttMapping.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class PHCompositeNode;

class InttCombinedRawDataDecoder : public SubsysReco
{
 public:
  InttCombinedRawDataDecoder(std::string const& name = "InttCombinedRawDataDecoder");

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  int LoadHotChannelMapLocal(std::string const& = "INTT_HotChannelMap.root");
  int LoadHotChannelMapRemote(std::string const& = "INTT_HotChannelMap");

 private:
  std::string m_InttRawNodeName = "INTTRAWHIT";
  typedef std::set<InttNameSpace::RawData_s, InttNameSpace::RawDataComparator> Set_t;
  Set_t m_HotChannelSet;
};

#endif  // INTT_COMBINEDRAWDATADECODER_H
