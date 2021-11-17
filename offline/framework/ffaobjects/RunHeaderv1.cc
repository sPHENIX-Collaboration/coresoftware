#include "RunHeaderv1.h"

#include <cmath>    // for NAN
#include <utility>  // for pair

void RunHeaderv1::identify(std::ostream &out) const
{
  out << "identify yourself: I am an RunHeaderv1 Object" << std::endl;
  out << "Run no: " << RunNumber << std::endl;
  auto iter = m_IntRunProperties.begin();
  while (iter != m_IntRunProperties.end())
  {
    out << iter->first << ": " << iter->second << std::endl;
    ++iter;
  }
  auto iterf = m_FloatRunProperties.begin();
  while (iterf != m_FloatRunProperties.end())
  {
    out << iterf->first << ": " << iterf->second << std::endl;
    ++iterf;
  }
  return;
}

int RunHeaderv1::isValid() const
{
  return ((RunNumber) ? 1 : 0);  // return 1 if runnumber not zero
}

void RunHeaderv1::set_floatval(const std::string &name, const float fval)
{
  m_FloatRunProperties[name] = fval;
}

float RunHeaderv1::get_floatval(const std::string &name) const
{
  std::map<std::string, float>::const_iterator iter = m_FloatRunProperties.find(name);
  if (iter != m_FloatRunProperties.end())
  {
    return iter->second;
  }
  return NAN;
}

void RunHeaderv1::set_intval(const std::string &name, const int ival)
{
  m_IntRunProperties[name] = ival;
}

int RunHeaderv1::get_intval(const std::string &name) const
{
  std::map<std::string, int>::const_iterator iter = m_IntRunProperties.find(name);
  if (iter != m_IntRunProperties.end())
  {
    return iter->second;
  }
  return -999999;
}
