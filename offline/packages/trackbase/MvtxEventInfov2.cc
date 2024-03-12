#include "MvtxEventInfov2.h"

#include <phool/phool.h>

PHObject* MvtxEventInfov2::CloneMe() const
{
  std::cout << PHWHERE << "::" << __func__ << " is not implemented in daughter class" << std::endl;
  return nullptr;
}

void MvtxEventInfov2::Reset()
{
  m_strobe_BCOs.clear();
  m_strobe_BCO_L1_BCO.clear();
  m_StringEventProperties.clear();
  m_IntEventProperties.clear();
  m_Int64EventProperties.clear();
  m_UintEventProperties.clear();
  m_Uint64EventProperties.clear();
  m_FloatEventProperties.clear();

  return;
}

void MvtxEventInfov2::identify(std::ostream& out) const
{
  out << "MvtxEventInfov2 information" << std::endl;

  for (const auto& m_StringEventPropertie : m_StringEventProperties)
  {
    out << m_StringEventPropertie.first << ": " << m_StringEventPropertie.second << std::endl;
  }

  out << "List of strobe BCOs:" << std::endl;
  std::set<uint64_t> strobeList = get_strobe_BCOs();
  for (unsigned long iter : strobeList)
  {
    out << iter << std::endl;
  }

  out << "Number of L1 triggers in this event" << std::endl;
  out << m_number_L1_name << ": " << get_number_L1s() << std::endl;

  out << "Number of heart beats in this event" << std::endl;
  out << m_number_HB_name << ": " << get_number_HB() << std::endl;

  if (get_number_L1s() > 0)
  {
    out << "List of strobe BCOs and L1 BCOs in this event" << std::endl;

    for (unsigned long iterStrobe : strobeList)
    {
      std::set<uint64_t> l1List = get_L1_BCO_from_strobe_BCO(iterStrobe);
      out << "Strobe BCO: " << iterStrobe << std::endl;
      for (unsigned long iterL1 : l1List)
      {
        out << "L1 BCO: " << iterL1 << std::endl;
      }
      out << std::endl;
    }
  }

  for (const auto& m_IntEventPropertie : m_IntEventProperties)
  {
    out << m_IntEventPropertie.first << ": " << m_IntEventPropertie.second << std::endl;
  }

  for (const auto& m_Int64EventPropertie : m_Int64EventProperties)
  {
    out << m_Int64EventPropertie.first << ": " << m_Int64EventPropertie.second << std::endl;
  }

  for (const auto& m_UintEventPropertie : m_UintEventProperties)
  {
    out << m_UintEventPropertie.first << ": " << m_UintEventPropertie.second << std::endl;
  }

  for (const auto& m_Uint64EventPropertie : m_Uint64EventProperties)
  {
    out << m_Uint64EventPropertie.first << ": " << m_Uint64EventPropertie.second << std::endl;
  }

  for (const auto& m_FloatEventPropertie : m_FloatEventProperties)
  {
    out << m_FloatEventPropertie.first << ": " << m_FloatEventPropertie.second << std::endl;
  }
  return;
}

int MvtxEventInfov2::isValid() const
{
  std::cout << PHWHERE << " isValid not implemented by daughter class" << std::endl;
  return 0;
}

void MvtxEventInfov2::set_number_HB(const int ival)
{
  m_number_HB = ival;
}

int MvtxEventInfov2::get_number_HB() const
{
  return m_number_HB;
}

void MvtxEventInfov2::set_strobe_BCO(const uint64_t strobe_BCO)
{
  m_strobe_BCOs.insert(strobe_BCO);
}

void MvtxEventInfov2::set_strobe_BCO_L1_BCO(const uint64_t strobe_BCO, const uint64_t L1_BCO)
{
  strobe_L1_pair strobe_L1_BCO_pair(strobe_BCO, L1_BCO);
  m_strobe_BCO_L1_BCO.insert(strobe_L1_BCO_pair);
}

unsigned int MvtxEventInfov2::get_number_strobes() const
{
  std::set<uint64_t> mySet;
  std::set<strobe_L1_pair>::const_iterator iter;
  for (iter = m_strobe_BCO_L1_BCO.begin(); iter != m_strobe_BCO_L1_BCO.end(); ++iter)
  {
    strobe_L1_pair myPair = *iter;
    mySet.insert(myPair.first);
  }

  return mySet.size();
}

unsigned int MvtxEventInfov2::get_number_L1s() const
{
  std::set<uint64_t> mySet;
  std::set<strobe_L1_pair>::const_iterator iter;
  for (iter = m_strobe_BCO_L1_BCO.begin(); iter != m_strobe_BCO_L1_BCO.end(); ++iter)
  {
    strobe_L1_pair myPair = *iter;
    mySet.insert(myPair.second);
  }

  return mySet.size();
}

std::set<uint64_t> MvtxEventInfov2::get_strobe_BCOs() const
{
  std::set<uint64_t> mySet = m_strobe_BCOs;
  return mySet;
}

std::set<uint64_t> MvtxEventInfov2::get_L1_BCOs() const
{
  std::set<uint64_t> mySet;
  std::set<strobe_L1_pair>::const_iterator iter;
  for (iter = m_strobe_BCO_L1_BCO.begin(); iter != m_strobe_BCO_L1_BCO.end(); ++iter)
  {
    strobe_L1_pair myPair = *iter;
    mySet.insert(myPair.second);
  }

  return mySet;
}

std::set<uint64_t> MvtxEventInfov2::get_strobe_BCO_from_L1_BCO(const uint64_t ival) const
{
  std::set<uint64_t> mySet;
  std::set<strobe_L1_pair>::const_iterator iter;
  for (iter = m_strobe_BCO_L1_BCO.begin(); iter != m_strobe_BCO_L1_BCO.end(); ++iter)
  {
    strobe_L1_pair myPair = *iter;
    if (ival == myPair.second)
    {
      mySet.insert(myPair.first);
    }
  }

  return mySet;
}

std::set<uint64_t> MvtxEventInfov2::get_L1_BCO_from_strobe_BCO(const uint64_t ival) const
{
  std::set<uint64_t> mySet;
  std::set<strobe_L1_pair>::const_iterator iter;
  for (iter = m_strobe_BCO_L1_BCO.begin(); iter != m_strobe_BCO_L1_BCO.end(); ++iter)
  {
    strobe_L1_pair myPair = *iter;
    if (ival == myPair.first)
    {
      mySet.insert(myPair.second);
    }
  }

  return mySet;
}
