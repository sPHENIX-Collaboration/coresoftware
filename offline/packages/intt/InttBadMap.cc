#include "InttBadMap.h"

#include <cdbobjects/CDBTTree.h>

#include <iostream>
#include <limits>
#include <map>
#include <string>

InttBadMap::~InttBadMap()
{
  delete m_bad_channel_set;
}

bool InttBadMap::IsBad(
    InttMap::RawData_s const& raw) const
{
  if (!m_bad_channel_set)
  {
    return false;
  }
  return m_bad_channel_set->find(raw) != m_bad_channel_set->end();
}

int InttBadMap::LoadFromCDBTTree(
    CDBTTree& cdbttree)
{
  delete m_bad_channel_set;
  m_bad_channel_set = new Set_t;
  m_bad_channel_set->clear();

  struct InttMap::RawData_s raw;
  std::map<std::string, InttMap::field_t*> raw_fields = {
    {"felix_server",  &raw.pid},
    {"felix_channel", &raw.fee},
    {"chip",          &raw.chp},
    {"channel",       &raw.chn},
  };

  for(int n = 0, N = cdbttree.GetSingleIntValue("size"); n < N; ++n)
  {
    for(auto& p : raw_fields)
    {
      if((*p.second = cdbttree.GetIntValue(n, p.first)) != std::numeric_limits<int>::min())
      {
        continue;
      }

      std::cerr << __PRETTY_FUNCTION__ << "\n"
                << "\tCDBTTree::GetIntValue returned std::numeric_limits<int>::min()\n"
                << "\t" << p.first << "[" << n << "] combination does not exist in \"cdbttree\"\n" << std::endl;
      return 1;
    }

	raw.pid += 3001;
    m_bad_channel_set->insert(raw);
  }

  return 0;
}

