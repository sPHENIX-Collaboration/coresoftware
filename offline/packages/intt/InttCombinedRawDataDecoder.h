#ifndef INTT_COMBINEDRAWDATADECODER_H
#define INTT_COMBINEDRAWDATADECODER_H

#include "InttMapping.h"
#include "InttDacMap.h"
#include "InttBCOMap.h"

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
  enum CalibRef {
    CDB  = 0,
    FILE = 1,
  };

  InttCombinedRawDataDecoder(std::string const& name = "InttCombinedRawDataDecoder");

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  int LoadHotChannelMapLocal(std::string const& = "INTT_HotChannelMap.root");
  int LoadHotChannelMapRemote(std::string const& = "INTT_HotChannelMap");

  void SetCalibDAC(std::string const& calibname= "INTT_DACMAP", const CalibRef& calibref=CDB) 
               { m_calibinfoDAC = std::pair< std::string, CalibRef>(calibname, calibref); }

  void SetCalibBCO(std::string const& calibname= "INTT_BCOMAP", const CalibRef& calibref=CDB) 
               { m_calibinfoBCO = std::pair< std::string, CalibRef>(calibname, calibref); }


  void runInttStandalone(bool runAlone) { m_runStandAlone = runAlone; }

  void writeInttEventHeader(bool write) { m_writeInttEventHeader = write; }

 private:
  InttEventInfo* intt_event_header = nullptr;
  std::string m_InttRawNodeName = "INTTRAWHIT";
  typedef std::set<InttNameSpace::RawData_s, InttNameSpace::RawDataComparator> Set_t;
  Set_t m_HotChannelSet;
  bool m_runStandAlone = false;
  bool m_writeInttEventHeader = false;

  std::pair<std::string, CalibRef> m_calibinfoDAC;
  std::pair<std::string, CalibRef> m_calibinfoBCO;

  InttDacMap          m_dacmap;
  InttBCOMap          m_bcomap;
};

#endif  // INTT_COMBINEDRAWDATADECODER_H
