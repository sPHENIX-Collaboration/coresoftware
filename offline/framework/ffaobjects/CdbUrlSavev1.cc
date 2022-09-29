#include "CdbUrlSavev1.h"

#include <phool/phool.h>

#include <iostream>

class PHObject;

PHObject *
CdbUrlSavev1::CloneMe() const
{
  std::cout << "CdbUrlSavev1::CloneMe() is not implemented in daugther class" << std::endl;
  return nullptr;
}

void CdbUrlSavev1::Reset()
{
  std::cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << std::endl;
  return;
}

void CdbUrlSavev1::identify(std::ostream &os) const
{
  os << "identify yourself: CdbUrlSavev1 Object" << std::endl;
  for (auto &iter : m_CdbUrlVector)
  {
    os << "domain: " << std::get<0>(iter)
       << ", url: " << std::get<1>(iter)
       << ", timestamp: " << std::get<2>(iter) << std::endl;
  }
  return;
}

int CdbUrlSavev1::isValid() const
{
  if (!m_CdbUrlVector.empty())
  {
    return 1;
  }
  return 0;
}

std::vector<std::tuple<std::string, std::string, uint64_t>>::const_iterator CdbUrlSavev1::begin() const
{
  return m_CdbUrlVector.begin();
}

std::vector<std::tuple<std::string, std::string, uint64_t>>::const_iterator CdbUrlSavev1::end() const
{
  return m_CdbUrlVector.end();
}

void CdbUrlSavev1::AddUrl(const std::string &domain, const std::string &url, const uint64_t timestamp)
{
  m_CdbUrlVector.emplace_back(domain, url, timestamp);
}

void CdbUrlSavev1::AddUrl(const std::tuple<std::string, std::string, uint64_t> &tup)
{
  m_CdbUrlVector.push_back(tup);
}
