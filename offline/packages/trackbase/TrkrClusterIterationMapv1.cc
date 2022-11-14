/**                                                                                                                                                                                                       
 * @file trackbase/TrkrClusterIterationMap.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief TrkrClusterIterationMap implementation
 */

#include "TrkrClusterIterationMapv1.h"

#include <ostream>  // for operator<<, endl, basic_ostream, ostream, basic_o...

void TrkrClusterIterationMapv1::Reset()
{
  m_map.clear();
}

void TrkrClusterIterationMapv1::identify(std::ostream& os) const
{
  os << "-----TrkrClusterIterationMapv1-----" << std::endl;
  os << "Number of associations: " << m_map.size() << std::endl;

  for (auto& entry : m_map)
  {
    // os << "   cluster key: 0x" << std::hex << entry.first << std::dec
    os << "   cluster key: " << entry.first << std::dec
       << " iteration: " << entry.second << std::endl;
  }

  os << "------------------------------" << std::endl;

  return;
}

void TrkrClusterIterationMapv1::addIteration(TrkrDefs::cluskey ckey, short int iter)
{
  m_map.insert(std::make_pair(ckey, iter));
}

short int TrkrClusterIterationMapv1::getIteration(TrkrDefs::cluskey ckey)
{
  Map::iterator iter = m_map.find(ckey);
  if (iter != m_map.end())
  {
    return (*iter).second;
  }
  else
  {
    return 0;
  }
}

unsigned int TrkrClusterIterationMapv1::size() const { return m_map.size(); }
