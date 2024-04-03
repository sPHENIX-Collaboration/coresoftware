#include "JetMapv1.h"

#include "Jet.h"

#include <phool/PHObject.h>  // for PHObject

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>  // for reverse_iterator
#include <ostream>   // for operator<<, endl, ostream, basic_ostream::operat...
#include <utility>   // for pair, make_pair
#include <vector>

JetMapv1::JetMapv1(const JetMap& jets)
  : _algo(jets.get_algo())
  , _par(jets.get_par())
{
  for (ConstSrcIter iter = jets.begin_src();
       iter != jets.end_src();
       ++iter)
  {
    _src.insert(*iter);
  }

  for (auto iter : jets)
  {
    Jet* jet = dynamic_cast<Jet*>((iter.second)->CloneMe());
    assert(jet);
    _map.insert(std::make_pair(jet->get_id(), jet));
  }
}

JetMapv1& JetMapv1::operator=(const JetMap& jets)
{
  Reset();

  _algo = jets.get_algo();
  _par = jets.get_par();

  for (ConstSrcIter iter = jets.begin_src();
       iter != jets.end_src();
       ++iter)
  {
    _src.insert(*iter);
  }

  for (auto iter : jets)
  {
    Jet* jet = dynamic_cast<Jet*>((iter.second)->CloneMe());
    assert(jet);
    _map.insert(std::make_pair(jet->get_id(), jet));
  }

  return *this;
}

JetMapv1::~JetMapv1()
{
  JetMapv1::Reset();
}

void JetMapv1::Reset()
{
  _algo = Jet::NONE;
  _par = NAN;
  _src.clear();

  while (_map.begin() != _map.end())
  {
    delete _map.begin()->second;
    _map.erase(_map.begin());
  }
}

void JetMapv1::identify(std::ostream& os) const
{
  os << "JetMapv1: size = " << _map.size() << std::endl;
  os << "          par = " << _par << std::endl;
  os << "          source = ";
  for (ConstSrcIter i = begin_src(); i != end_src(); ++i)
  {
    os << (*i) << ",";
  }
  os << std::endl;

  return;
}

const Jet* JetMapv1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return nullptr;
  return iter->second;
}

Jet* JetMapv1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end()) return nullptr;
  return iter->second;
}

Jet* JetMapv1::insert(Jet* jet)
{
  unsigned int index = 0;
  if (!_map.empty()) index = _map.rbegin()->first + 1;
  _map.insert(std::make_pair(index, jet));
  _map[index]->set_id(index);
  return (_map[index]);
}

std::vector<Jet*> JetMapv1::vec()
{
  std::vector<Jet*> v_data;
  for (auto& _ : _map)
  {
    v_data.push_back(_.second);
  }
  return v_data;
}
