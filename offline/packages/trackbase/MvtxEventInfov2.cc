#include "MvtxEventInfov2.h"

#include <phool/phool.h>

PHObject* MvtxEventInfov2::CloneMe() const
{
  std::cout << PHWHERE << "::" << __func__ << " is not implemented in daughter class" << std::endl;
  return nullptr;
}

void MvtxEventInfov2::Reset()
{
  m_strobe_BCO_L1_BCO.clear();
  m_StringEventProperties.clear();
  m_IntEventProperties.clear();
  m_Int64EventProperties.clear();
  m_UintEventProperties.clear();
  m_Uint64EventProperties.clear();
  m_FloatEventProperties.clear();
  m_strobe_BCO = 0;

  return;
}

void MvtxEventInfov2::identify(std::ostream &out) const
{
  out << "MvtxEventInfov2 information" << std::endl;

  for (auto iters = m_StringEventProperties.begin(); iters != m_StringEventProperties.end(); ++iters)
  {
    out << iters->first << ": " << iters->second << std::endl;
  }
  
  out << "Strobe BCO: " << m_strobe_BCO << std::endl;
 
  out << "Number of L1 triggers in this event" << std::endl;
  out << m_number_L1_name << ": " << get_number_L1s() << std::endl;
    
  out << "Number of heart beats in this event" << std::endl;
  out << m_number_HB_name << ": " << get_number_HB() << std::endl;

  if (get_number_L1s() > 0)
  {
    out << "List of strobe BCOs and L1 BCOs in this event" << std::endl;

    std::set<uint64_t> strobeList = get_strobe_BCOs();
    for (auto iterStrobe = strobeList.begin(); iterStrobe != strobeList.end(); ++iterStrobe)
    { 
      std::set<uint64_t> l1List = get_L1_BCO_from_strobe_BCO(*iterStrobe);
      out << "Strobe BCO: " << *iterStrobe << std::endl; 
      for (auto iterL1 = l1List.begin(); iterL1 != l1List.end(); ++iterL1)
      {
        out << "L1 BCO: " << *iterL1 << std::endl; 
      }
      out << std::endl;
    }
  }

  for (auto iteri = m_IntEventProperties.begin(); iteri != m_IntEventProperties.end(); ++iteri)
  {
    out << iteri->first << ": " << iteri->second << std::endl;
  }

  for (auto iteri64 = m_Int64EventProperties.begin(); iteri64 != m_Int64EventProperties.end(); ++iteri64)
  {
    out << iteri64->first << ": " << iteri64->second << std::endl;
  }

  for (auto iteru = m_UintEventProperties.begin(); iteru != m_UintEventProperties.end(); ++iteru)
  {
    out << iteru->first << ": " << iteru->second << std::endl;
  }

  for (auto iteru64 = m_Uint64EventProperties.begin(); iteru64 != m_Uint64EventProperties.end(); ++iteru64)
  {
    out << iteru64->first << ": " << iteru64->second << std::endl;
  }

  for (auto iterf = m_FloatEventProperties.begin(); iterf != m_FloatEventProperties.end(); ++iterf)
  {
    out << iterf->first << ": " << iterf->second << std::endl;
  }
  return;
}

int MvtxEventInfov2::isValid() const
{
  std::cout << PHWHERE << " isValid not implemented by daughter class" << std::endl;
  return 0;
}
