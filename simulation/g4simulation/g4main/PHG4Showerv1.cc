#include "PHG4Showerv1.h"

#include "PHG4HitDefs.h"

#include <algorithm>  // for fill
#include <cmath>
#include <iostream>
#include <iterator>  // for begin, end
#include <utility>

PHG4Showerv1::PHG4Showerv1()
{
  // with C++11 begin() and end() exist, you do not need the array size anymore
  // to fill an array
  std::fill(std::begin(_pos), std::end(_pos), std::numeric_limits<float>::quiet_NaN());
  std::fill(std::begin(_covar), std::end(_covar), std::numeric_limits<float>::quiet_NaN());
}

void PHG4Showerv1::identify(std::ostream &os) const
{
  os << "---PHG4Showerv1-------------------------------" << std::endl;
  os << "id: " << get_id() << std::endl;
  os << "parent_particle_id: " << get_parent_particle_id() << std::endl;
  os << "parent_shower_id: " << get_parent_shower_id() << std::endl;
  os << "x: " << get_x() << std::endl;
  os << "y: " << get_y() << std::endl;
  os << "z: " << get_z() << std::endl;

  os << "         ( ";
  os << get_covar(0, 0) << " , ";
  os << get_covar(0, 1) << " , ";
  os << get_covar(0, 2) << " )" << std::endl;
  os << " covar = ( ";
  os << get_covar(1, 0) << " , ";
  os << get_covar(1, 1) << " , ";
  os << get_covar(1, 2) << " )" << std::endl;
  os << "         ( ";
  os << get_covar(2, 0) << " , ";
  os << get_covar(2, 1) << " , ";
  os << get_covar(2, 2) << " )" << std::endl;

  os << "VOLUME ID : edep eion light_yield" << std::endl;
  for (auto iter : _edep)
  {
    int volid = iter.first;
    os << volid << " : " << get_edep(volid) << " " << get_eion(volid) << " "
       << get_light_yield(volid) << std::endl;
  }

  os << "G4Particle IDs" << std::endl;
  for (int _g4particle_id : _g4particle_ids)
  {
    os << _g4particle_id << " ";
  }
  os << std::endl;

  os << "G4Hit IDs" << std::endl;
  for (const auto &_g4hit_id : _g4hit_ids)
  {
    for (unsigned long long jter : _g4hit_id.second)
    {
      os << jter << " ";
    }
  }
  os << std::endl;

  os << "-----------------------------------------------" << std::endl;

  return;
}

int PHG4Showerv1::isValid() const
{
  if (_id == 0)
  {
    return 0;
  }
  for (float _po : _pos)
  {
    if (std::isnan(_po))
    {
      return 0;
    }
  }
  for (int j = 0; j < 3; ++j)
  {
    for (int i = j; i < 3; ++i)
    {
      if (std::isnan(get_covar(i, j)))
      {
        return 0;
      }
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
  {
    std::swap(i, j);
  }
  return i + 1 + (j + 1) * (j) / 2 - 1;
}

unsigned int PHG4Showerv1::get_nhits(int volume) const
{
  std::map<int, unsigned int>::const_iterator citer =
      _nhits.find(volume);
  if (citer == _nhits.end())
  {
    return 0;
  }
  return citer->second;
}

float PHG4Showerv1::get_edep(int volume) const
{
  std::map<int, float>::const_iterator citer =
      _edep.find(volume);
  if (citer == _edep.end())
  {
    return 0.0;
  }
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
  {
    return 0.0;
  }
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
  {
    return 0.0;
  }
  return citer->second;
}

float PHG4Showerv1::get_eh_ratio(int volume) const
{
  std::map<int, float>::const_iterator citer =
      _eh_ratio.find(volume);
  if (citer == _eh_ratio.end())
  {
    return 0.0;
  }
  return citer->second;
}
