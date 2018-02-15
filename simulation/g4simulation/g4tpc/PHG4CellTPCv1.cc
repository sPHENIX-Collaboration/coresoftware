#include "PHG4CellTPCv1.h"

#include <phool/phool.h>

#include <iostream>

using namespace std;

PHG4CellTPCv1::PHG4CellTPCv1():
  cellid(~0x0)
{}

PHG4CellTPCv1::PHG4CellTPCv1(const PHG4CellDefs::keytype g4cellid):
  cellid(g4cellid)
{}

PHG4CellTPCv1::~PHG4CellTPCv1()
{
  timeseq.clear();
  prop_map.clear();
  return;
}

void
PHG4CellTPCv1::add_edep(const PHG4HitDefs::keytype g4hitid, const int tbin, const float edep)
{
// If the returned boolean from the map.insert is
// true: the element was inserted
// false: the element existed 
// in any case an iterator to the element in the map is returned
 
    EdepMap edepmap;
    pair<map<int,EdepMap>::iterator, bool> ret = timeseq.insert(make_pair(tbin,edepmap));
     map<int,EdepMap>::iterator mapiter = ret.first;
// insert 0 edep so edep can be added no matter if hits was found or not
  pair<EdepIterator, bool> edepret = (mapiter->second).insert(make_pair(g4hitid,0.));
  EdepIterator eiter = edepret.first;
  eiter->second += edep;
  return;
}

bool
PHG4CellTPCv1::has_binning(const PHG4CellDefs::CellBinning binning) const
{
  return PHG4CellDefs::has_binning(cellid, binning);
}

short int
PHG4CellTPCv1::get_detid() const
{
  return PHG4CellDefs::get_detid(cellid);
}

bool
PHG4CellTPCv1::has_property(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return i!=prop_map.end();
}

float
PHG4CellTPCv1::get_property_float(const PROPERTY prop_id) const
{
  if (!check_property(prop_id,type_float))
    {
      pair<const string,PROPERTY_TYPE> property_info =get_property_info(prop_id); 
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second) 
	   << " not " << get_property_type(type_float) << endl; 
      exit(1);
    }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i!=prop_map.end()) return u_property(i->second).fdata;

  return   NAN ;
}

int
PHG4CellTPCv1::get_property_int(const PROPERTY prop_id) const
{
  if (!check_property(prop_id,type_int))
    {
      pair<const string,PROPERTY_TYPE> property_info =get_property_info(prop_id); 
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second) 
	   << " not " << get_property_type(type_int) << endl; 
      exit(1);
    }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i!=prop_map.end()) return u_property(i->second).idata;

  return INT_MIN;
}

unsigned int
PHG4CellTPCv1::get_property_uint(const PROPERTY prop_id) const
{
  if (!check_property(prop_id,type_uint))
    {
      pair<const string,PROPERTY_TYPE> property_info =get_property_info(prop_id); 
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second) 
	   << " not " << get_property_type(type_uint) << endl; 
      exit(1);
    }
  prop_map_t::const_iterator i = prop_map.find(prop_id);

  if (i!=prop_map.end()) return u_property(i->second).uidata;

  return UINT_MAX ;
}

void
PHG4CellTPCv1::add_property(const PROPERTY prop_id, const float value)
{
  if (!check_property(prop_id,type_float))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id); 
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second) 
	   << " not " << get_property_type(type_float) << endl; 
      exit(1);
    }
  float val = value;
  if (prop_map.find(prop_id) != prop_map.end())
    {
      val += get_property_float(prop_id);
    }
  prop_map[prop_id] = u_property(val).get_storage();
}

void
PHG4CellTPCv1::add_property(const PROPERTY prop_id, const int value)
{
  if (!check_property(prop_id,type_int))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id); 
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second) 
	   << " not " << get_property_type(type_int) << endl; 
      exit(1);
    }
  int val = value;
  if (prop_map.find(prop_id) != prop_map.end())
    {
      val += get_property_int(prop_id);
    }
  prop_map[prop_id] += u_property(val).get_storage();
}

void
PHG4CellTPCv1::add_property(const PROPERTY prop_id, const unsigned int value)
{
  if (!check_property(prop_id,type_uint))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id); 
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second) 
	   << " not " << get_property_type(type_uint) << endl; 
      exit(1);
    }
  unsigned int val = value;
  if (prop_map.find(prop_id) != prop_map.end())
    {
      val += get_property_uint(prop_id);
    }
  prop_map[prop_id] += u_property(val).get_storage();
}

void
PHG4CellTPCv1::set_property(const PROPERTY prop_id, const float value)
{
  if (!check_property(prop_id,type_float))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id); 
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second) 
	   << " not " << get_property_type(type_float) << endl; 
      exit(1);
    }
  prop_map[prop_id] = u_property(value).get_storage();
}

void
PHG4CellTPCv1::set_property(const PROPERTY prop_id, const int value)
{
  if (!check_property(prop_id,type_int))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id); 
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second) 
	   << " not " << get_property_type(type_int) << endl; 
      exit(1);
    }
  prop_map[prop_id] += u_property(value).get_storage();
}

void
PHG4CellTPCv1::set_property(const PROPERTY prop_id, const unsigned int value)
{
  if (!check_property(prop_id,type_uint))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id); 
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second) 
	   << " not " << get_property_type(type_uint) << endl; 
      exit(1);
    }
  prop_map[prop_id] += u_property(value).get_storage();
}

unsigned int
PHG4CellTPCv1::get_property_nocheck(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator iter = prop_map.find(prop_id);
  if (iter != prop_map.end())
    {
      return iter->second;
    }
  return UINT_MAX;
}

void
PHG4CellTPCv1::print() const
{
  identify(cout);
}

void
PHG4CellTPCv1::Reset()
{
  timeseq.clear();
  prop_map.clear();
  return;
}

void PHG4CellTPCv1::identify(std::ostream& os) const
{
  os << "New PHG4CellTPCv1  0x" << hex << cellid << dec << endl;

//  os <<"Associated to "<<hitedeps.size()<<" hits"<<endl;
  for (const auto pair :timeseq)
  {
    os <<" TimeBin "<< pair.first << " size: " << (pair.second).size() << endl;
    for (const auto pair2 :pair.second)
    {
      os << "\t PHG4Hit " << pair2.first << " -> " << pair2.second << " GeV" << endl;
    }
  }


//  os <<"Contains to "<<trainOfDigits.size()<<" TPC digitization chain"<<endl;

  for (prop_map_t::const_iterator i = prop_map.begin(); i != prop_map.end(); ++i)
  {
    PROPERTY prop_id = static_cast<PROPERTY>(i->first);
    pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
    os << "\t" << prop_id << ":\t" << property_info.first << " = \t";
    switch (property_info.second)
    {
    case type_int:
      os << get_property_int(prop_id);
      break;
    case type_uint:
      os << get_property_uint(prop_id);
      break;
    case type_float:
      os << get_property_float(prop_id);
      break;
    default:
      os << " unknown type ";
    }
    os << endl;
  }
}
