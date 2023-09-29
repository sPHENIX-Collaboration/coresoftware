#include "MvtxEventInfov1.h"

//#include <cmath>    // for NAN
//#include <utility>  // for pair
#include <phool/phool.h>

PHObject* MvtxEventInfov1::CloneMe() const
{
  std::cout << PHWHERE << "::" << __func__ << " is not implemented in daughter class" << std::endl;
  return nullptr;
}

void MvtxEventInfov1::Reset()
{
  m_strobe_BCO.clear();
  m_number_LL1 = -1;
  m_StringEventProperties.clear();
  m_IntEventProperties.clear();
  m_Int64EventProperties.clear();
  m_UintEventProperties.clear();
  m_Uint64EventProperties.clear();
  m_FloatEventProperties.clear();

  return;
}


void MvtxEventInfov1::identify(std::ostream &out) const
{
  out << "MvtxEventInfov1 information" << std::endl;

  auto iters = m_StringEventProperties.begin();
  while (iters != m_StringEventProperties.end())
  {
    out << iters->first << ": " << iters->second << std::endl;
    ++iters;
  }

  out << "List of strobe BCOs in this event" << std::endl;
  auto iterBCO = m_strobe_BCO.begin();
  while (iterBCO != m_strobe_BCO.end())
  {
    out << iterBCO->first << ": " << iterBCO->second << std::endl;
    ++iterBCO;
  }
    
  out << "Number of LL1 triggers in this event" << std::endl;
  out << m_number_LL1_name << ": " << m_number_LL1 << std::endl;

  auto iteri = m_IntEventProperties.begin();
  while (iteri != m_IntEventProperties.end())
  {
    out << iteri->first << ": " << iteri->second << std::endl;
    ++iteri;
  }

  auto iteri64 = m_Int64EventProperties.begin();
  while (iteri64 != m_Int64EventProperties.end())
  {
    out << iteri64->first << ": " << iteri64->second << std::endl;
    ++iteri64;
  }

  auto iteru = m_UintEventProperties.begin();
  while (iteru != m_UintEventProperties.end())
  {
    out << iteru->first << ": " << iteru->second << std::endl;
    ++iteru;
  }

  auto iteru64 = m_Uint64EventProperties.begin();
  while (iteru64 != m_Uint64EventProperties.end())
  {
    out << iteru64->first << ": " << iteru64->second << std::endl;
    ++iteru64;
  }

  auto iterf = m_FloatEventProperties.begin();
  while (iterf != m_FloatEventProperties.end())
  {
    out << iterf->first << ": " << iterf->second << std::endl;
    ++iterf;
  }
  return;
}

int MvtxEventInfov1::isValid() const
{
  std::cout << PHWHERE << " isValid not implemented by daughter class" << std::endl;
  return 0;
}

void MvtxEventInfov1::set_number_LL1(const int ival)
{
  m_number_LL1 = ival;
}

int MvtxEventInfov1::get_number_LL1() const
{
  return m_number_LL1;
}

void MvtxEventInfov1::set_strobe_BCO(const uint64_t ival)
{
  m_strobe_BCO[m_strobe_BCO.size()] = ival;
}

uint64_t MvtxEventInfov1::get_strobe_BCO(const uint32_t ival) const
{
  std::map<uint32_t, uint64_t>::const_iterator iter = m_strobe_BCO.find(ival);
  if (iter != m_strobe_BCO.end())
  {
    return iter->second;
  }
  return 0;
}

uint32_t MvtxEventInfov1::get_number_strobe_BCO() const
{
  return m_strobe_BCO.size();
}
