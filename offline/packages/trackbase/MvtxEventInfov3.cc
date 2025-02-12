#include "MvtxEventInfov3.h"

#include <phool/phool.h>

PHObject* MvtxEventInfov3::CloneMe() const
{
  std::cout << PHWHERE << "::" << __func__ << " is not implemented in daughter class" << std::endl;
  return nullptr;
}

void MvtxEventInfov3::Reset()
{
  m_strobe_BCOs.clear();
  m_L1_BCOs.clear();
  m_Int64EventProperties.clear();
  m_IntEventProperties.clear();
  m_Int64EventProperties.clear();
  m_UintEventProperties.clear();
  m_Uint64EventProperties.clear();
  m_FloatEventProperties.clear();

  return;
}

void MvtxEventInfov3::identify(std::ostream& out) const
{
  out << "MvtxEventInfov1 information" << std::endl;

  for (const auto& m_StringEventPropertie : m_StringEventProperties)
  {
    out << m_StringEventPropertie.first << ": " << m_StringEventPropertie.second << std::endl;
  }

  if (get_number_strobes() > 0)
  {
    out << "List of strobe BCOs in this event" << std::endl;

    for (unsigned long iterStrobe : m_strobe_BCOs)
    {
      out << "Strobe BCO: " << iterStrobe << std::endl;
      out << std::endl;
    }
  }

  if (get_number_L1s() > 0)
  {
    out << "List of L1 BCOs in this event" << std::endl;

    for (unsigned long iterStrobe : m_L1_BCOs)
    {
      out << "L1 BCO: " << iterStrobe << std::endl;
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

int MvtxEventInfov3::isValid() const
{
  std::cout << PHWHERE << " isValid not implemented by daughter class" << std::endl;
  return 0;
}
