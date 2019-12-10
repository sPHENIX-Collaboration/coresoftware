#include "PHG4Showerv1.h"

#include "PHG4HitDefs.h"

#include <algorithm>      // for fill
#include <cmath>
#include <iostream>
#include <iterator>       // for begin, end
#include <utility>

using namespace std;

PHG4Showerv1::PHG4Showerv1()
  : _id(0xFFFFFFFF)
  , _parent_particle_id(0)
  , _parent_shower_id(0)
  , _pos()
  , _covar()
  , _edep()
  , _eion()
  , _light_yield()
  , _eh_ratio()
  , _g4particle_ids()
  , _g4hit_ids()
{
  // with C++11 begin() and end() exist, you do not need the array size anymore
  // to fill an array
  fill(std::begin(_pos), std::end(_pos), NAN);
  fill(std::begin(_covar), std::end(_covar), NAN);
}

void PHG4Showerv1::identify(ostream &os) const
{
  os << "---PHG4Showerv1-------------------------------" << endl;
  os << "id: " << get_id() << endl;
  os << "parent_particle_id: " << get_parent_particle_id() << endl;
  os << "parent_shower_id: " << get_parent_shower_id() << endl;
  os << "x: " << get_x() << endl;
  os << "y: " << get_y() << endl;
  os << "z: " << get_z() << endl;

  os << "         ( ";
  os << get_covar(0, 0) << " , ";
  os << get_covar(0, 1) << " , ";
  os << get_covar(0, 2) << " )" << endl;
  os << " covar = ( ";
  os << get_covar(1, 0) << " , ";
  os << get_covar(1, 1) << " , ";
  os << get_covar(1, 2) << " )" << endl;
  os << "         ( ";
  os << get_covar(2, 0) << " , ";
  os << get_covar(2, 1) << " , ";
  os << get_covar(2, 2) << " )" << endl;

  os << "VOLUME ID : edep eion light_yield" << endl;
  for (std::map<int, float>::const_iterator iter = _edep.begin();
       iter != _edep.end(); ++iter)
  {
    int volid = iter->first;
    os << volid << " : " << get_edep(volid) << " " << get_eion(volid) << " "
       << get_light_yield(volid) << endl;
  }

  os << "G4Particle IDs" << endl;
  for (std::set<int>::const_iterator iter = _g4particle_ids.begin();
       iter != _g4particle_ids.end(); ++iter)
  {
    os << *iter << " ";
  }
  os << endl;

  os << "G4Hit IDs" << endl;
  for (std::map<int, std::set<PHG4HitDefs::keytype> >::const_iterator iter =
           _g4hit_ids.begin();
       iter != _g4hit_ids.end(); ++iter)
  {
    for (std::set<PHG4HitDefs::keytype>::const_iterator jter =
             iter->second.begin();
         jter != iter->second.end(); ++jter)
    {
      os << *jter << " ";
    }
  }
  os << endl;

  os << "-----------------------------------------------" << endl;

  return;
}

int PHG4Showerv1::isValid() const
{
  if (_id == 0)
    return 0;
  for (int i = 0; i < 3; ++i)
  {
    if (isnan(_pos[i]))
      return 0;
  }
  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      if (isnan(get_covar(i, j)))
        return 0;
    }
  }
  return 1;
}

void PHG4Showerv1::set_covar(unsigned int i, unsigned int j, float value)
{
  _covar[covar_index(i, j)] = value;
  return;
}

float PHG4Showerv1::get_covar(unsigned int i, unsigned int j) const
{
  return _covar[covar_index(i, j)];
}

unsigned int PHG4Showerv1::covar_index(unsigned int i, unsigned int j) const
{
  if (i > j)
    std::swap(i, j);
  return i + 1 + (j + 1) * (j) / 2 - 1;
}

unsigned int PHG4Showerv1::get_nhits(int volume) const
{
  std::map<int, unsigned int>::const_iterator citer =
      _nhits.find(volume);
  if (citer == _nhits.end())
    return 0;
  return citer->second;
}

float PHG4Showerv1::get_edep(int volume) const
{
  std::map<int, float>::const_iterator citer =
      _edep.find(volume);
  if (citer == _edep.end())
    return 0.0;
  return citer->second;
}

double PHG4Showerv1::get_edep() const
{
  double sum = 0;
  for (const auto &iter : _edep)
  {
    sum += iter.second;
  }
  return sum;
}

float PHG4Showerv1::get_eion(int volume) const
{
  std::map<int, float>::const_iterator citer =
      _eion.find(volume);
  if (citer == _eion.end())
    return 0.0;
  return citer->second;
}

double PHG4Showerv1::get_eion() const
{
  double sum = 0;
  for (const auto &iter : _eion)
  {
    sum += iter.second;
  }
  return sum;
}

float PHG4Showerv1::get_light_yield(int volume) const
{
  std::map<int, float>::const_iterator citer =
      _light_yield.find(volume);
  if (citer == _light_yield.end())
    return 0.0;
  return citer->second;
}

float PHG4Showerv1::get_eh_ratio(int volume) const
{
  std::map<int, float>::const_iterator citer =
      _eh_ratio.find(volume);
  if (citer == _eh_ratio.end())
    return 0.0;
  return citer->second;
}
