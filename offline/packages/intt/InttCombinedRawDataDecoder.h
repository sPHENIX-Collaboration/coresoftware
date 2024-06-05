#ifndef INTT_COMBINEDRAWDATADECODER_H
#define INTT_COMBINEDRAWDATADECODER_H

#include "InttBadChannelMap.h"
#include "InttBCOMap.h"
#include "InttDacMap.h"
#include "InttFeeMap.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class PHCompositeNode;
class InttEventInfo;

class InttCombinedRawDataDecoder : public SubsysReco
{
 public:
  InttCombinedRawDataDecoder(std::string const& name = "InttCombinedRawDataDecoder");

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void runInttStandalone(bool runAlone) { m_runStandAlone = runAlone; }
  void writeInttEventHeader(bool write) { m_writeInttEventHeader = write; }

  int LoadBadMap(std::string const& name = "INTT_HotChannelMap") {return m_badmap.Load(name);}
  int LoadBcoMap(std::string const& name = "INTT_BCOMAP") {return m_bcomap.Load(name);}
  int LoadDacMap(std::string const& name = "INTT_DACMAP") {return m_dacmap.Load(name);}
  int LoadFeeMap(std::string const& name = "INTT_FEEMAP") {return m_feemap.Load(name);}

 private:
  InttEventInfo* intt_event_header{nullptr};
  std::string m_InttRawNodeName = "INTTRAWHIT";

  bool m_runStandAlone{false};
  bool m_writeInttEventHeader{false};

  InttBadChannelMap m_badmap;
  InttBCOMap m_bcomap;
  InttDacMap m_dacmap;
  InttFeeMap m_feemap;
};

#endif  // INTT_COMBINEDRAWDATADECODER_H
