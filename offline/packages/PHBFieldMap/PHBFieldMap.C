
#include <iostream>
#include "PHBFieldMap.hh"
#include "PHBFieldMapObj.h"

using namespace std;

PHBFieldMap::PHBFieldMap(const std::string& mapfile)
{

  _scale_factor = 1.0;

#ifndef MAP_VOXEL
  _fieldmap = (void *) new PHBFieldMapObj(mapfile,false);
  const PHBFieldMapObj::BField& bval = ((PHBFieldMapObj *)_fieldmap)->get_field(0.0,0.0,0.0);
#else
  _fieldmap = (void *) new PHBFieldMapObj(mapfile,true);
  const PHBFieldMapObj::BField& bval = ((PHBFieldMapObj *)_fieldmap)->get_field_voxel(0.0,0.0,0.0);
#endif
  
  cout << "Magnetic field at 0,0,0 (Gauss) = " 
       <<  _scale_factor*bval.bx << " " 
       <<  _scale_factor*bval.by << " "
       <<  _scale_factor*bval.bz << endl;

}

PHBFieldMap::~PHBFieldMap()
{
  delete( (PHBFieldMapObj *)_fieldmap );  
}

bool
PHBFieldMap::get_bfield(const double *loc, double *field)
{

#ifdef MAP_VOXEL
  const PHBFieldMapObj::BField& bval = ((PHBFieldMapObj *)_fieldmap)->get_field_voxel(loc[0],loc[1],loc[2]);
#else
  const PHBFieldMapObj::BField& bval = ((PHBFieldMapObj *)_fieldmap)->get_field(loc[0],loc[1],loc[2]);
#endif

  field[0] = _scale_factor*bval.bx;
  field[1] = _scale_factor*bval.by;
  field[2] = _scale_factor*bval.bz;

  return True;

}

void 
PHBFieldMap::set_scale_factor(double new_factor){ 

  _scale_factor = new_factor; 

}
