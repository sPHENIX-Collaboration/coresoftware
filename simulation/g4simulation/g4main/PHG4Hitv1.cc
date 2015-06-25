#include "PHG4Hitv1.h"
#include "PHG4HitDefs.h"

using namespace std;

//G4Allocator<PHG4Hitv1> PHG4Hitv1Allocator;

ClassImp(PHG4Hitv1)

PHG4Hitv1::PHG4Hitv1():
 hitid(UINT_MAX),
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
  std::cout<<"New Hitv1   "<<hitid<<"  on track "<<trackid<<" EDep "<<edep<<std::endl;
  std::cout<<"Location: X "<<x[0]<<"/"<<x[1]<<"  Y "<<y[0]<<"/"<<y[1]<<"  Z "<<z[0]<<"/"<<z[1]<<std::endl;
  std::cout<<"Time        "<<t[0]<<"/"<<t[1]<<std::endl;

  for (prop_map_t::const_iterator i = prop_map.begin(); i!= prop_map.end(); ++i)
    {
      u_property p( i->second);

      std::cout <<"\t" << static_cast<const int>(i->first) <<":\t" <<get_property_name(static_cast<PROPERTY>(i->first))<<" = "
          <<"\t"<<p.fdata<<" (float)"
          <<"\t"<<p.idata<<" (int)"
          <<"\t"<<p.uidata<<" (uint)"
          <<"\t"<<sizeof(u_property)<<" (memory size)"
          <<endl;
    }
}

bool
PHG4Hitv1::has_property(PHG4Hitv1::PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return i!=prop_map.end();
}

float
PHG4Hitv1::get_property_float(PHG4Hitv1::PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return (i==prop_map.end())? NAN : u_property(i->second).fdata ;
}

int
PHG4Hitv1::get_property_int(PHG4Hitv1::PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return (i==prop_map.end())? NAN : u_property(i->second).idata ;
}

unsigned int
PHG4Hitv1::get_property_uint(PHG4Hitv1::PROPERTY prop_id) const
{
  prop_map_t::const_iterator i = prop_map.find(prop_id);
  return (i==prop_map.end())? NAN : u_property(i->second).uidata ;
}

void
PHG4Hitv1::set_property(PHG4Hitv1::PROPERTY prop_id, float value)
{
  prop_map[prop_id] = u_property(value).get_storage();
}

void
PHG4Hitv1::set_property(PHG4Hitv1::PROPERTY prop_id, int value)
{
  prop_map[prop_id] = u_property(value).get_storage();
}

void
PHG4Hitv1::set_property(PHG4Hitv1::PROPERTY prop_id, unsigned int value)
{
  prop_map[prop_id] = u_property(value).get_storage();
}
