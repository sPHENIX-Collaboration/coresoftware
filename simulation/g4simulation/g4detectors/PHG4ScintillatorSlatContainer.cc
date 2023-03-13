#include "PHG4ScintillatorSlatContainer.h"

#include "PHG4ScintillatorSlat.h"  // for PHG4ScintillatorSlat
#include "PHG4ScintillatorSlatDefs.h"

#include <cstdlib>

void PHG4ScintillatorSlatContainer::Reset()
{
  while (slatmap.begin() != slatmap.end())
  {
    delete slatmap.begin()->second;
    slatmap.erase(slatmap.begin());
  }
  return;
}

void PHG4ScintillatorSlatContainer::identify(std::ostream& os) const
{
  std::map<unsigned int, PHG4ScintillatorSlat*>::const_iterator iter;
  os << "Number of slats: " << size() << std::endl;
  for (iter = slatmap.begin(); iter != slatmap.end(); ++iter)
  {
    os << "slat key 0x" << std::hex << iter->first << std::dec << std::endl;
    (iter->second)->identify();
  }
  return;
}

PHG4ScintillatorSlatContainer::ConstIterator
PHG4ScintillatorSlatContainer::AddScintillatorSlat(const PHG4ScintillatorSlatDefs::keytype key, PHG4ScintillatorSlat* newslat)
{
  if (slatmap.find(key) != slatmap.end())
  {
    std::cout << "PHG4ScintillatorSlatContainer::AddScintillatorSlatSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
  newslat->set_key(key);
  slatmap[key] = newslat;
  return slatmap.find(key);
}

PHG4ScintillatorSlatContainer::ConstRange
PHG4ScintillatorSlatContainer::getScintillatorSlats(const short icolumn) const
{
  if ((icolumn >> PHG4ScintillatorSlatDefs::columnbits) > 0)
  {
    std::cout << " column id too large: " << icolumn << std::endl;
    exit(1);
  }
  //  unsigned int shiftval = detid << slat_idbits;
  PHG4ScintillatorSlatDefs::keytype keylow = icolumn << PHG4ScintillatorSlatDefs::columnbits;
  PHG4ScintillatorSlatDefs::keytype keyup = ((icolumn + 1) << PHG4ScintillatorSlatDefs::columnbits) - 1;
  //   cout << "keylow: 0x" << hex << keylow << dec << std::endl;
  //   cout << "keyup: 0x" << hex << keyup << dec << std::endl;
  ConstRange retpair;
  retpair.first = slatmap.lower_bound(keylow);
  retpair.second = slatmap.upper_bound(keyup);
  return retpair;
}

PHG4ScintillatorSlatContainer::ConstRange
PHG4ScintillatorSlatContainer::getScintillatorSlats() const
{
  return std::make_pair(slatmap.begin(), slatmap.end());
}

PHG4ScintillatorSlat*
PHG4ScintillatorSlatContainer::findScintillatorSlat(PHG4ScintillatorSlatDefs::keytype key)
{
  PHG4ScintillatorSlatContainer::ConstIterator it = slatmap.find(key);

  if (it != slatmap.end())
  {
    return it->second;
  }

  return nullptr;
}

double
PHG4ScintillatorSlatContainer::getTotalEdep() const
{
  ConstIterator iter;
  double totalenergy = 0;
  for (iter = slatmap.begin(); iter != slatmap.end(); ++iter)
  {
    totalenergy += iter->second->get_edep();
  }
  return totalenergy;
}
