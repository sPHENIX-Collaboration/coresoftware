#ifndef INTT_COMBINEDRAWDATADECODER_H
#define INTT_COMBINEDRAWDATADECODER_H

#include "InttBadChannelMap.h"
#include "InttBCOMap.h"
#include "InttDacMap.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/SubsysReco.h>

#include <set>
#include <string>
#include <vector>
#include <map>

class PHCompositeNode;
class InttEventInfo;

class InttCombinedRawDataDecoder : public SubsysReco
{
 public:
  enum CalibRef
  {
    CDB = 0,
    FILE = 1,
  };

  InttCombinedRawDataDecoder(std::string const& name = "InttCombinedRawDataDecoder");

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  /// Overloaded; no arguments loads with default tag
  int LoadBadChannelMap() {return m_badmap.Load();}
  int LoadBadChannelMap(std::string const& s) {return m_badmap.Load(s);}

  /// Depreciated; use LoadHotChannelMap(const std::string&);
  int LoadHotChannelMapLocal(std::string const& s = "INTT_HotChannelMap.root") {return LoadBadChannelMap(s);}
  /// Depreciated; use LoadHotChannelMap(const std::string&);
  int LoadHotChannelMapRemote(std::string const& s = "INTT_HotChannelMap") {return LoadBadChannelMap(s);}

  void SetCalibDAC(std::string const& calibname = "INTT_DACMAP", const CalibRef& calibref = CDB)
  {
    m_calibinfoDAC = std::pair<std::string, CalibRef>(calibname, calibref);
  }

  void SetCalibBCO(std::string const& calibname = "INTT_BCOMAP", const CalibRef& calibref = CDB)
  {
    m_calibinfoBCO = std::pair<std::string, CalibRef>(calibname, calibref);
  }
  void useRawHitNodeName(const std::string& name) { m_InttRawNodeName = name; }
  void runInttStandalone(bool runAlone) { m_runStandAlone = runAlone; }

  void writeInttEventHeader(bool write) { m_writeInttEventHeader = write; }

  void set_inttFeeOffset(int offset) { m_inttFeeOffset = offset; }
  void set_outputBcoDiff(bool flag) {m_outputBcoDiff = flag; }
  void set_triggeredMode(bool flag) {m_triggeredMode = flag; }
  void set_bcoFilter(bool flag) {m_bcoFilter = flag; }
  void set_SaturatedChipRejection(bool flag){m_SaturatedChipRejection = flag;} // note : this is for removing a fraction of the saturated chips
  void set_HighChipMultiplicityCut(int cut){HighChipMultiplicityCut = cut;}

 private:
  InttEventInfo* intt_event_header = nullptr;
  std::string m_InttRawNodeName = "INTTRAWHIT";
  bool m_runStandAlone = false;
  bool m_writeInttEventHeader = false;
  bool m_bcoFilter = false;
  bool m_SaturatedChipRejection = true; // note : true as default
  std::pair<std::string, CalibRef> m_calibinfoDAC;
  std::pair<std::string, CalibRef> m_calibinfoBCO;

  InttBadChannelMap m_badmap;
  InttDacMap m_dacmap;
  InttBCOMap m_bcomap;

  int m_inttFeeOffset = 23;   //23 is the offset for INTT in streaming mode
  bool m_outputBcoDiff = false;
  bool m_triggeredMode = false;

  std::vector<std::string> evt_inttHits_vec;
  std::map<std::string, int> evt_ChipHit_count_map;
  int HighChipMultiplicityCut = 71;

};

#endif  // INTT_COMBINEDRAWDATADECODER_H
