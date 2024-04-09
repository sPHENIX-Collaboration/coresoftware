#ifndef INTT_COMBINEDRAWDATADECODER_H
#define INTT_COMBINEDRAWDATADECODER_H

#include "InttMap.h"
#include "InttBadMap.h"
#include "InttBcoMap.h"
#include "InttDacMap.h"
#include "InttFeeMap.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/SubsysReco.h>

#include <set>
#include <map>
#include <string>

class PHCompositeNode;
class InttEventInfo;

class InttCombinedRawDataDecoder : public SubsysReco
{
 public:
  enum class calib_load_e : char
  {
    SKIP = 'S',
    FROM_CDB = 'C',
    FROM_FILE = 'F',
  };
  static constexpr auto SKIP =      calib_load_e::SKIP;
  static constexpr auto FROM_CDB =  calib_load_e::FROM_CDB;
  static constexpr auto FROM_FILE = calib_load_e::FROM_FILE;

  InttCombinedRawDataDecoder(std::string const& name = "InttCombinedRawDataDecoder");

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  int LoadHotChannelMapLocal(std::string const& = "INTT_HotChannelMap.root");
  int LoadHotChannelMapRemote(std::string const& = "INTT_HotChannelMap");

  int SetCalib(std::string const&, calib_load_e const&, std::string const& = "");
  int ClearCalib(std::string const&);
  int ClearCalibs();

  void runInttStandalone(bool runAlone) { m_runStandAlone = runAlone; }
  void writeInttEventHeader(bool write) { m_writeInttEventHeader = write; }

 private:
  std::string m_InttRawNodeName = "INTTRAWHIT";

  bool m_writeInttEventHeader = false;
  bool m_runStandAlone = false;

  struct calib_load_s {
    calib_load_e method = SKIP;
    std::string filename = "";
	InttLoadable* const ptr = nullptr;
  };
  typedef std::map<std::string, struct calib_load_s> calib_map_t;
  calib_map_t m_calibs;
  InttBadMap m_badmap;
  InttBcoMap m_bcomap;
  InttDacMap m_dacmap;
  InttFeeMap m_feemap;

  int LoadCalibs();
};

#endif  // INTT_COMBINEDRAWDATADECODER_H
