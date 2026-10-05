/**
 * @file trackbase/LaserAggregatedPadContainerv1.cc
 * @author Ben Kimelman
 * @date February 2024
 * @brief Implementation of LaserAggregatedPadContainerv1
 */
#include "LaserAggregatedPadContainerv1.h"
#include "LaserAggregatedPad.h"

#include <phool/phool.h>

#include <cstdlib>

void LaserAggregatedPadContainerv1::Reset()
{
  for( auto&& [pad, newPad]:m_padMap )
  { delete newPad; }
  
  m_padMap.clear();

  nLaserEvents = 0.0;
}

void LaserAggregatedPadContainerv1::identify(std::ostream& os) const
{
  os << "-----LaserAggregatedPadContainerv1-----" << std::endl;
  ConstIterator iter;
  os << "Number of pads: " << size() << std::endl;
  for (iter = m_padMap.begin(); iter != m_padMap.end(); ++iter)
  {
    os << "pad side: " << padkey::get_side(iter->first) << "   layer: " << padkey::get_layer(iter->first) << "   iphi: " << padkey::get_phibin(iter->first) << std::endl;
    (iter->second)->identify();
  }
  os << "------------------------------" << std::endl;
  return;
}

void LaserAggregatedPadContainerv1::addPad(const padkey::key pad, LaserAggregatedPad* newPad)
{
  auto ret = m_padMap.insert(std::make_pair(pad, newPad));
  if ( !ret.second )
  {
    std::cout << "LaserAggregatedPadContainerv1::AddPad: duplicate pad side: " << padkey::get_side(pad) << "   layer: " << padkey::get_layer(pad) << "   iphi: " << padkey::get_phibin(pad) << ". exiting now" << std::endl;
    exit(1);
  }
}

void LaserAggregatedPadContainerv1::removePad(padkey::key pad)
{ 
  auto *padToRem = findPad(pad);
  delete padToRem;

  m_padMap.erase(pad); 
}

LaserAggregatedPadContainer::ConstRange
LaserAggregatedPadContainerv1::getPads() const
{ return std::make_pair(m_padMap.cbegin(), m_padMap.cend()); }

LaserAggregatedPad*
LaserAggregatedPadContainerv1::findPad(padkey::key pad) const
{
  auto it = m_padMap.find(pad);
  return it == m_padMap.end() ? nullptr:it->second;
}

unsigned int LaserAggregatedPadContainerv1::size() const
{
  return m_padMap.size();
}

void LaserAggregatedPadContainerv1::merge(const LaserAggregatedPadContainer *other)
{
  if (!other || other == this){ return; }

  const auto range = other->getPads();
  auto hint = m_padMap.begin();
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    hint = m_padMap.lower_bound(iter->first);   // or advance manually
    if (hint != m_padMap.end() && hint->first == iter->first)
    {
      hint->second->addAdc(iter->second->getAdc());
      hint->second->addNHits(iter->second->getNHits());
    }
    else
    {
      auto *clone = dynamic_cast<LaserAggregatedPad*>(iter->second->CloneMe());
      if (clone){ hint = m_padMap.emplace_hint(hint, iter->first, clone); }
    }
  }

  setNLaserEvents(getNLaserEvents() + other->getNLaserEvents());
}
