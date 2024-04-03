#include "InttBCOMap.h"

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <filesystem>
#include <iostream>

InttBCOMap::InttBCOMap()
{
  for(int felix_server = 0; felix_server<8;felix_server++)
  {
    for(int felix_channel = 0;felix_channel<14;felix_channel++)
    {
      m_bco[felix_server][felix_channel] = -1;
    }
  }
}

int InttBCOMap::LoadFromCDB(std::string const &calibname)
{
  if (calibname.empty())
  {
    std::cout << "InttBCOMap::LoadFromCDB(std::string const& name)" << std::endl;
    std::cout << "\tArgument 'name' is empty string" << std::endl;
    return -1;
  }

  std::string database = CDBInterface::instance()->getUrl(calibname);
  if (database.empty())
  {
    std::cout << "InttBCOMap::LoadFromCDB(std::string const& name)" << std::endl;
    std::cout << "\tArgument 'database' is empty string. calibname invalid :" << calibname << std::endl;
    return -1;
  }

  return LoadFromFile(database);
}

int InttBCOMap::LoadFromFile(std::string const &filename)
{
  if (filename.empty())
  {
    std::cout << "InttBCOMap::LoadFromFile(std::string const& name)" << std::endl;
    std::cout << "\tArgument 'filename' is empty string" << std::endl;
    return 1;
  }

  if (!std::filesystem::exists(filename))
  {
    std::cout << "int InttBCOMap::LoadFromFile(std::string const& filename)" << std::endl;
    std::cout << "\tFile '" << filename << "' does not exist" << std::endl;
    return 1;
  }

  std::cout << "CDBFile: " << filename << std::endl;

  CDBTTree cdbttree = CDBTTree(filename);
  cdbttree.LoadCalibrations();

  return LoadFromCDBTTree(cdbttree);
}

int InttBCOMap::LoadFromCDBTTree(CDBTTree &cdbttree)
{
  ///////////////
  std::cout << "LoadFromCDBTTree::LoadFromCDBTTree" << std::endl;

  uint64_t N = cdbttree.GetSingleIntValue("size");
  for (uint64_t n = 0; n < N; ++n)
  {
    int felix_server = cdbttree.GetIntValue(n, "felix_server");
    int felix_channel = cdbttree.GetIntValue(n, "felix_channel");
    int bco_diff = cdbttree.GetIntValue(n, "bco_diff");
    m_bco[felix_server][felix_channel] = bco_diff;

    if(m_verbosity>0){
      std::cout << "felix_server " << felix_server << " ";
      std::cout << "felix_channel " << felix_channel << " ";
      std::cout << "bco_diff " << bco_diff << std::endl;
    }
  }
  return 0;
}
bool InttBCOMap::IsBad(const int &felix_server, const int &felix_channel, uint64_t const &bco_full, const int &bco)
{
  int bco_diff = (bco_full & 0x7FU) - bco;
  if (bco_diff < 0)
  {
    bco_diff += 128;
  }
  //////////////////////////////////////////////////////////////////////////////
  // Hits belongs to [peak+1,peak-1] (3BCO region) will survive after BCO cut //
  //////////////////////////////////////////////////////////////////////////////
  int bco_peak = m_bco[felix_server][felix_channel];
  int bco_minus = bco_peak - 1;
  if (bco_minus == -1)
  {
    bco_minus = 127;
  }
  int bco_plus = bco_peak + 1;
  if (bco_plus == 128)
  {
    bco_plus = 0;
  }

// -1: m_bco is initial value, not load the parameter. accept all bco
  if( bco_peak == -1 || bco_diff == bco_peak || bco_diff == bco_minus || bco_diff == bco_plus)
  { 
    //std::cout<<"m_bco is initial value, not load the parameter. accept all bco "<<felix_server<<" "<<felix_channel<<std::endl;
    return false;
  }
  else
  {
    return true;
  }
}
bool InttBCOMap::IsBad(InttNameSpace::RawData_s const &raw, uint64_t const &bco_full, const int &bco)
{
  return IsBad(raw.felix_server, raw.felix_channel, bco_full, bco);
}

bool InttBCOMap::IsBad(InttNameSpace::Offline_s const &off, uint64_t const &bco_full, const int &bco)
{
  return IsBad(InttNameSpace::ToRawData(off), bco_full, bco);
}
