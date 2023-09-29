#include "MvtxEventInfo.h"

#include <cmath>    // for NAN
#include <utility>  // for pair
#include <phool/phool.h>

PHObject* MvtxEventInfo::CloneMe() const
{
  std::cout << PHWHERE << "::" << __func__ << " is not implemented in daughter class" << std::endl;
  return nullptr;
}

void MvtxEventInfo::Reset()
{
  std::cout << PHWHERE << "::" << __func__ << " is not implemented in daughter class" << std::endl;
  return;
}


void MvtxEventInfo::identify(std::ostream &out) const
{
  out << "MvtxEventInfo information" << std::endl;

  auto iters = m_StringEventProperties.begin();
  while (iters != m_StringEventProperties.end())
  {
    out << iters->first << ": " << iters->second << std::endl;
    ++iters;
  }

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

int MvtxEventInfo::isValid() const
{
  std::cout << PHWHERE << " isValid not implemented by daughter class" << std::endl;
  return 0;
}

void MvtxEventInfo::set_floatval(const std::string &name, const float fval)
{
  m_FloatEventProperties[name] = fval;
}

float MvtxEventInfo::get_floatval(const std::string &name) const
{
  std::map<std::string, float>::const_iterator iter = m_FloatEventProperties.find(name);
  if (iter != m_FloatEventProperties.end())
  {
    return iter->second;
  }
  return NAN;
}

void MvtxEventInfo::set_intval(const std::string &name, const int32_t ival)
{
  m_IntEventProperties[name] = ival;
}

int MvtxEventInfo::get_intval(const std::string &name) const
{
  std::map<std::string, int32_t>::const_iterator iter = m_IntEventProperties.find(name);
  if (iter != m_IntEventProperties.end())
  {
    return iter->second;
  }
  return -999999;
}

void MvtxEventInfo::set_int64val(const std::string &name, const int64_t ival)
{
  m_Int64EventProperties[name] = ival;
}

int MvtxEventInfo::get_int64val(const std::string &name) const
{
  std::map<std::string, int64_t>::const_iterator iter = m_Int64EventProperties.find(name);
  if (iter != m_Int64EventProperties.end())
  {
    return iter->second;
  }
  return -999999;
}

void MvtxEventInfo::set_uintval(const std::string &name, const uint32_t ival)
{
  m_UintEventProperties[name] = ival;
}

int MvtxEventInfo::get_uintval(const std::string &name) const
{
  std::map<std::string, uint32_t>::const_iterator iter = m_UintEventProperties.find(name);
  if (iter != m_UintEventProperties.end())
  {
    return iter->second;
  }
  return 999999;
}

void MvtxEventInfo::set_uint64val(const std::string &name, const uint64_t ival)
{
  m_Int64EventProperties[name] = ival;
}

int MvtxEventInfo::get_uint64val(const std::string &name) const
{
  std::map<std::string, uint64_t>::const_iterator iter = m_Uint64EventProperties.find(name);
  if (iter != m_Uint64EventProperties.end())
  {
    return iter->second;
  }
  return 999999;
}

void MvtxEventInfo::set_stringval(const std::string &name, const std::string &sval)
{
  m_StringEventProperties[name] = sval;
}

std::string MvtxEventInfo::get_stringval(const std::string &name) const
{
  std::map<std::string, std::string>::const_iterator iter = m_StringEventProperties.find(name);
  if (iter != m_StringEventProperties.end())
  {
    return iter->second;
  }
  return "Empty";
}
