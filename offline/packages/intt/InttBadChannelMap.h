#ifndef INTT_BAD_CHANNEL_MAP_H
#define INTT_BAD_CHANNEL_MAP_H

#include "InttMapping.h"

#include <iostream>
#include <string>
#include <set>

class CDBTTree;

class InttBadChannelMap
{
 public:
  InttBadChannelMap() = default;
  virtual ~InttBadChannelMap() = default;

  int Load(std::string const& = "INTT_HotMap"); // Should match CDB tag

  /// Depreciated; use int Load(const& std::string)
  int LoadFromFile(std::string const& s = "InttBadChannelMap.root") {return Load(s);}
  /// Depreciated; use int Load(const& std::string)
  int LoadFromCDB(std::string const& s =  "INTT_HotMap") {return Load(s);}

  void identify(std::ostream& = std::cout) const;
  void Print(std::ostream& = std::cout) const;
  int size() const {return m_size;}

  bool IsBad(InttNameSpace::Online_s const&) const;
  bool IsBad(InttNameSpace::Offline_s const&) const;
  bool IsBad(InttNameSpace::RawData_s const&) const;

  /// If the calibration has been loaded using offline indexing conventions
  bool OfflineLoaded() {return m_offline_loaded;}

  /// If the calibration has been loaded using rawdata indexing conventions
  bool RawDataLoaded() {return m_rawdata_loaded;}

 protected:
  virtual int v_LoadFromCDBTTree(CDBTTree&);

 private:
  typedef std::set<InttNameSpace::Offline_s> OfflineSet_t;
  typedef std::set<InttNameSpace::RawData_s> RawDataSet_t;
  OfflineSet_t m_offline_set;
  RawDataSet_t m_rawdata_set;

  bool m_offline_loaded{false};
  bool m_rawdata_loaded{false};

  int m_size{0};
};

#endif  // INTT_BAD_CHANNEL_MAP_H

