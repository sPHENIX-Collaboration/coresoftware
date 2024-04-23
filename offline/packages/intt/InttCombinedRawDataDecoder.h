#ifndef INTT_COMBINEDRAWDATADECODER_H
#define INTT_COMBINEDRAWDATADECODER_H

#include "InttBadMap.h"
#include "InttBcoMap.h"
#include "InttDacMap.h"
#include "InttFeeMap.h"
#include "InttMapping.h"

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

  int LoadBadMap(std::string const& s = "") { return m_badmap.Load(s); }
  int LoadBcoMap(std::string const& s = "") { return m_bcomap.Load(s); }
  int LoadDacMap(std::string const& s = "") { return m_dacmap.Load(s); }
  int LoadFeeMap(std::string const& s = "") { return m_feemap.Load(s); }

  void runInttStandalone(bool runAlone) { m_runStandAlone = runAlone; }
  void writeInttEventHeader(bool write) { m_writeInttEventHeader = write; }

 private:
  InttEventInfo* intt_event_header = nullptr;
  std::string m_InttRawNodeName = "INTTRAWHIT";

  bool m_runStandAlone = false;
  bool m_writeInttEventHeader = false;

  InttBadMap m_badmap;
  InttBcoMap m_bcomap;
  InttDacMap m_dacmap;
  InttFeeMap m_feemap;
};

#endif  // INTT_COMBINEDRAWDATADECODER_H
