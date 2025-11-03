#include "TruthNeutralMesonContainer.h"
#include "TruthNeutralMeson.h"

#include <iostream>

void TruthNeutralMesonContainer::Reset()
{
  for (auto& pair : _mesons)
  {
    delete pair.second;
  }
  _mesons.clear();
}

void TruthNeutralMesonContainer::identify(std::ostream& os) const
{
  os << "TruthNeutralMesonContainer with " << _mesons.size() << " mesons:" << std::endl;
  for (const auto& pair : _mesons)
  {
    os << "  key = " << pair.first << ": ";
    if (pair.second)
    {
      pair.second->identify(os);
    }
    else
    {
      os << "null pointer!" << std::endl;
    }
  }
}

TruthNeutralMesonContainer::ConstIterator TruthNeutralMesonContainer::AddMeson(TruthNeutralMeson* meson)
{
  unsigned int key = _mesons.size();
  return _mesons.insert(std::make_pair(key, meson)).first;
}

TruthNeutralMeson* TruthNeutralMesonContainer::getMeson(unsigned int id)
{
  auto it = _mesons.find(id);
  if (it != _mesons.end())
  {
    return it->second;
  }
  return nullptr;
}

const TruthNeutralMeson* TruthNeutralMesonContainer::getMeson(unsigned int id) const
{
  auto it = _mesons.find(id);
  if (it != _mesons.end())
  {
    return it->second;
  }
  return nullptr;
}

TruthNeutralMesonContainer::ConstRange TruthNeutralMesonContainer::getMesons() const
{
  return std::make_pair(_mesons.begin(), _mesons.end());
}

TruthNeutralMesonContainer::Range TruthNeutralMesonContainer::getMesons()
{
  return std::make_pair(_mesons.begin(), _mesons.end());
}
