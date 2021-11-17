#include "EventHeaderv1.h"

#include <cmath>  // for NAN
#include <iostream>
#include <utility>  // for pair

void EventHeaderv1::Reset()
{
  RunNumber = 0;
  EvtSequence = 0;
  m_IntEventProperties.clear();
  m_FloatEventProperties.clear();
  return;
}

void EventHeaderv1::identify(std::ostream &out) const
{
  out << "identify yourself: I am an EventHeaderv1 Object" << std::endl;
  out << "Run Number: " << RunNumber
      << ", Event no: " << EvtSequence
      << std::endl;
  auto iter = m_IntEventProperties.begin();
  while (iter != m_IntEventProperties.end())
  {
    out << iter->first << ": " << iter->second << std::endl;
    ++iter;
  }
  auto iterf = m_FloatEventProperties.begin();
  while (iterf != m_FloatEventProperties.end())
  {
    out << iterf->first << ": " << iterf->second << std::endl;
    ++iterf;
  }
  return;
}

int EventHeaderv1::isValid() const
{
  return ((RunNumber) ? 1 : 0);  // return 1 if runnumber is not zero
}

void EventHeaderv1::set_floatval(const std::string &name, const float fval)
{
  m_FloatEventProperties[name] = fval;
}

float EventHeaderv1::get_floatval(const std::string &name) const
{
  std::map<std::string, float>::const_iterator iter = m_FloatEventProperties.find(name);
  if (iter != m_FloatEventProperties.end())
  {
    return iter->second;
  }
  return NAN;
}

void EventHeaderv1::set_intval(const std::string &name, const int64_t ival)
{
  m_IntEventProperties[name] = ival;
}

int64_t EventHeaderv1::get_intval(const std::string &name) const
{
  auto iter = m_IntEventProperties.find(name);
  if (iter != m_IntEventProperties.end())
  {
    return iter->second;
  }
  return -999999;
}
