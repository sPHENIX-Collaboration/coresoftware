#include "PHG4Cell.h"

#include <cstdlib>

using namespace std;

void
PHG4Cell::Copy(PHG4Cell const &g4cell)
{
  set_cellid(g4cell.get_cellid());
  for (unsigned char ic = 0; ic < UCHAR_MAX; ic++)
    {
      PROPERTY prop_id = static_cast<PHG4Cell::PROPERTY> (ic);
      if (g4cell.has_property(prop_id))
        {
	  set_property_nocheck(prop_id,g4cell.get_property_nocheck(prop_id));
	}
    }
}


void
PHG4Cell::identify(ostream& os) const
{
  cout << "Class " << this->ClassName() << endl;
  return;
}

ostream& operator<<(ostream& stream, const PHG4Cell * cell){
  stream << "PHG4Cell"  << endl;
  return stream;
}

void
PHG4Cell::Reset()
{
  cout << "Reset not implemented by daughter class" << endl;
  return;
}

std::pair<const std::string,PHG4Cell::PROPERTY_TYPE> 
PHG4Cell::get_property_info(const PROPERTY prop_id) 
{
  switch (prop_id)
  {
  case  prop_edep:
    return make_pair("energy deposition",PHG4Cell::type_float);
  case  prop_eion:
    return make_pair("ionizing energy loss",PHG4Cell::type_float);
  case   prop_light_yield:
    return make_pair("light yield",PHG4Cell::type_float);
  default:
    cout << "PHG4Cell::get_property_info - Fatal Error - unknown index " << prop_id << endl;
    exit(1);
  }
}


bool
PHG4Cell::check_property(const PROPERTY prop_id, const PROPERTY_TYPE prop_type)
{
  pair<const string,PROPERTY_TYPE> property_info = get_property_info(prop_id);
  if (property_info.second != prop_type)
    {
      return false;
    }
  return true;
}

string
PHG4Cell::get_property_type(const PROPERTY_TYPE prop_type)
{
  switch(prop_type)
    {
    case type_int:
      return "int";
    case type_uint:
      return "unsigned int";
    case type_float:
      return "float";
    default:
      return "unkown";
    }
}
