#include "PHG4Hitv1.h"
#include "PHG4HitDefs.h"

#include <phool/phool.h>

#include <cstdlib>

using namespace std;

ClassImp(PHG4Hitv1)

PHG4Hitv1::PHG4Hitv1():
 hitid(ULONG_LONG_MAX),
 trackid(INT_MIN),
 edep(NAN)
{
  for (int i = 0; i<2;i++)
    {
      set_x(i,NAN);
      set_y(i,NAN);
      set_z(i,NAN);
      set_t(i,NAN);
    }

}

PHG4Hitv1::PHG4Hitv1(PHG4Hit const &g4hit)
{
  Copy(g4hit);
}

void
PHG4Hitv1::print() const {
  std::cout<<"New Hitv1  0x"<< hex << hitid 
	   << dec << "  on track "<<trackid<<" EDep "<<edep<<std::endl;
  std::cout<<"Location: X "<<x[0]<<"/"<<x[1]<<"  Y "<<y[0]<<"/"<<y[1]<<"  Z "<<z[0]<<"/"<<z[1]<<std::endl;
  std::cout<<"Time        "<<t[0]<<"/"<<t[1]<<std::endl;

  for (prop_map_t::const_iterator i = prop_map.begin(); i!= prop_map.end(); ++i)
    {
      PROPERTY prop_id = static_cast<PROPERTY>(i->first);
      pair<const string, PROPERTY_TYPE> property_info = get_property_info(prop_id);
      cout << "\t" << i->first << ":\t" << property_info.first << " = \t";
      switch(property_info.second)
	{
	case type_int:
	  cout << get_property_int(prop_id);
	  break;
	case type_uint:
	  cout << get_property_uint(prop_id);
	  break;
	case type_float:
	  cout << get_property_float(prop_id);
	  break;
	default:
	  cout << " unknown type ";
	}
      cout <<endl;
    }
}

bool
PHG4Hitv1::has_property(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return i!=prop_map.end();
}

float
PHG4Hitv1::get_property_float(const PROPERTY prop_id) const
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
  return (i==prop_map.end())? NAN : u_property(i->second).fdata ;
}

int
PHG4Hitv1::get_property_int(const PROPERTY prop_id) const
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
  return (i==prop_map.end())? NAN : u_property(i->second).idata ;
}

unsigned int
PHG4Hitv1::get_property_uint(const PROPERTY prop_id) const
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
  return (i==prop_map.end())? NAN : u_property(i->second).uidata ;
}

void
PHG4Hitv1::set_property(const PROPERTY prop_id, const float value)
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
PHG4Hitv1::set_property(const PROPERTY prop_id, const int value)
{
  if (!check_property(prop_id,type_int))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id); 
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second) 
	   << " not " << get_property_type(type_int) << endl; 
      exit(1);
    }
  prop_map[prop_id] = u_property(value).get_storage();
}

void
PHG4Hitv1::set_property(const PROPERTY prop_id, const unsigned int value)
{
  if (!check_property(prop_id,type_uint))
    {
      pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id); 
      cout << PHWHERE << " Property " << property_info.first << " with id "
           << prop_id << " is of type " << get_property_type(property_info.second) 
	   << " not " << get_property_type(type_uint) << endl; 
      exit(1);
    }
  prop_map[prop_id] = u_property(value).get_storage();
}

unsigned int
PHG4Hitv1::get_property_nocheck(const PROPERTY prop_id) const
{
  prop_map_t::const_iterator iter = prop_map.find(prop_id);
  if (iter != prop_map.end())
    {
      return iter->second;
    }
  return UINT_MAX;
}

float
PHG4Hitv1::get_px(const int i) const
{
  switch(i)
    {
    case 0:
      return  get_property_float(prop_px_0);
    case 1:
      return  get_property_float(prop_px_1);
    default:
      cout << "Invalid index in get_px: " << i << endl;
      exit(1);
    }
}

float
PHG4Hitv1::get_py(const int i) const
{
  switch(i)
    {
    case 0:
      return  get_property_float(prop_py_0);
    case 1:
      return  get_property_float(prop_py_1);
    default:
      cout << "Invalid index in get_py: " << i << endl;
      exit(1);
    }
}

float
PHG4Hitv1::get_pz(const int i) const
{
  switch(i)
    {
    case 0:
      return  get_property_float(prop_pz_0);
    case 1:
      return  get_property_float(prop_pz_1);
    default:
      cout << "Invalid index in get_pz: " << i << endl;
      exit(1);
    }
}

void
PHG4Hitv1::set_px(const int i, const float f)
{
  switch(i)
    {
    case 0:
      set_property(prop_px_0,f);
      return;
    case 1:
      set_property(prop_px_1,f);
      return;
    default:
      cout << "Invalid index in set_px: " << i << endl;
      exit(1);
    }
}

void
PHG4Hitv1::set_py(const int i, const float f)
{
  switch(i)
    {
    case 0:
      set_property(prop_py_0,f);
      return;
    case 1:
      set_property(prop_py_1,f);
      return;
    default:
      cout << "Invalid index in set_py: " << i << endl;
      exit(1);
    }
}

void
PHG4Hitv1::set_pz(const int i, const float f)
{
  switch(i)
    {
    case 0:
      set_property(prop_pz_0,f);
      return;
    case 1:
      set_property(prop_pz_1,f);
      return;
    default:
      cout << "Invalid index in set_pz: " << i << endl;
      exit(1);
    }
}
