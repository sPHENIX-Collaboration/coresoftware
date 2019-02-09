#include "RawTower_Temperature.h"
#include <calobase/RawTowerDefs.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <map>
#include <cassert>

#include "PROTOTYPE3_FEM.h"

using namespace std;

RawTower_Temperature::RawTower_Temperature() :
    towerid(~0) // initialize all bits on
{}

// we can copy only from another  RawTower_Temperature, not a generic tower
RawTower_Temperature::RawTower_Temperature(const RawTower_Temperature & tower)
{
  towerid = tower.get_id();

  for ( int i = 0; i < tower.get_nr_entries(); i++)
    {
      add_entry( tower.get_eventnumber_from_entry(i)
		 ,tower.get_time_from_entry(i)
		 ,tower.get_temperature_from_entry(i) );
    }
}

RawTower_Temperature::RawTower_Temperature(RawTowerDefs::keytype id) :
    towerid(id)
{

}

RawTower_Temperature::RawTower_Temperature(const unsigned int icol, const unsigned int irow)
{
  towerid = RawTowerDefs::encode_towerid(RawTowerDefs::NONE, icol, irow);
}


RawTower_Temperature::~RawTower_Temperature()
{
}

void
RawTower_Temperature::Reset()
{
  eventnumbers.clear();
  times.clear();
  temperatures.clear();
}

float  RawTower_Temperature::get_temperature_from_time(const time_t t) const
{
  if ( !isValid() ) return -1;

  if ( t < get_time_from_entry(0) ) // if we ask for a time before the start time, we return the first reading
    {
      return get_temperature_from_entry( 0);
    }


  int lowest_entry=0;
  int above_entry=0;

  for ( int i = 0; i < get_nr_entries() ; i++)
    {
      if ( get_time_from_entry(i) < t )
	{
	  lowest_entry = i;
	}
      else
	{
	  if ( !above_entry) above_entry = i;
	}
    }

  if ( !above_entry  ) // we didn't find a entry later than this
    {
      return get_temperature_from_entry( lowest_entry);
    }

  double m = ( get_temperature_from_entry(above_entry) - get_temperature_from_entry(lowest_entry)) /
    ( get_time_from_entry(above_entry) - get_time_from_entry(lowest_entry));

  return get_temperature_from_entry(lowest_entry) + m * ( t-get_time_from_entry(lowest_entry) );

}
void
RawTower_Temperature::identify(std::ostream& os) const
{
  os << "RawTower_Temperature col=" << get_column() << " row=" << get_row() << ":  " << temperatures.size() << " entries"  << std::endl;
}

void
RawTower_Temperature::print(std::ostream& os) const
{
  identify(os);

  cout << "entry    event     time        T" << endl;
  for ( int i = 0; i < get_nr_entries(); i++)
    {
      os << setw( 4) << i << "  " 
	 << setw(7) <<  get_eventnumber_from_entry(i) << "  "
	 << setw(7) <<  get_time_from_entry(i)        << "  "
	 << setw(5) <<  get_temperature_from_entry(i) 
	 << endl;
    }
}

