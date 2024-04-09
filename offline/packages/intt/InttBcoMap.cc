#include "InttBcoMap.h"

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <filesystem>
#include <iostream>
#include <map>

InttBcoMap::InttBcoMap()
{
  for(auto& felix_server : m_bco)
  {
    for(auto& felix_channel : felix_server)
	{
      felix_channel = -1;
	}
  }
}

int InttBcoMap::LoadFromCDBTTree(CDBTTree &cdbttree)
{
  int felix_server = 0;
  int felix_channel = 0;
  int bco_diff = 0;

  std::map<std::string, int*> fields = {
    {"felix_server",  &felix_server},
    {"felix_channel", &felix_channel},
    {"bco_diff",      &bco_diff},
  };

  for(int n = 0, N = cdbttree.GetSingleIntValue("size"); n < N; ++n)
  {
    for(auto& p : fields)
    {
      if((*p.second = cdbttree.GetIntValue(n, p.first)) != std::numeric_limits<int>::min())
      {
        continue;
      }

      std::cerr << __PRETTY_FUNCTION__ << "\n"
                << "\tCDBTTree::GetIntValue returned std::numeric_limits<int>::min()\n"
                << "\tname, channel combination does not exist in \"cdbttree\"\n" << std::endl;
      return 1;
    }

	try
	{
      m_bco.at(felix_server).at(felix_channel) = bco_diff;
	}
	catch (std::out_of_range const& e)
	{
      std::cerr << __PRETTY_FUNCTION__ << "\n"
                << "\t" << e.what() << std::endl;
	  return 1;
	}
  }

  return 0;
}

bool InttBcoMap::IsBad(
  int const& felix_server,
  int const& felix_channel,
  uint64_t const& bco_full,
  int const& bco)
{
  int bco_peak = 0;
  try
  {
    bco_peak = m_bco.at(felix_server).at(felix_channel);
  }
  catch (std::out_of_range const& e)
  {
    std::cerr << __PRETTY_FUNCTION__ << "\n"
              << "\t" << e.what() << std::endl;
	return false;
  }
  if(bco_peak == -1)
  {
    return false;
  }

  int lower = (bco_peak - WIDTH + 128) % 128;
  int upper = (bco_peak + WIDTH + 128) % 128;
  int bco_diff = ((bco_full & 0x7fU) - (bco % 128) + 128) % 128;

  if(lower < upper)
  {
  	return lower <= bco_diff && bco_diff <= upper;
  }

  return lower <= bco_diff || bco_diff <= upper;
}

bool InttBcoMap::IsBad(InttNameSpace::RawData_s const &raw, uint64_t const &bco_full, const int &bco)
{
  return IsBad(raw.felix_server, raw.felix_channel, bco_full, bco);
}

bool InttBcoMap::IsBad(InttNameSpace::Offline_s const &off, uint64_t const &bco_full, const int &bco)
{
  return IsBad(InttNameSpace::ToRawData(off), bco_full, bco);
}
