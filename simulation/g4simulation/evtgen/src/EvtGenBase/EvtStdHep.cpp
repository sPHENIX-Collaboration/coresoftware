
/***********************************************************************
* Copyright 1998-2020 CERN for the benefit of the EvtGen authors       *
*                                                                      *
* This file is part of EvtGen.                                         *
*                                                                      *
* EvtGen is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* EvtGen is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     *
***********************************************************************/

#include "EvtGenBase/EvtStdHep.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <iomanip>
#include <iostream>
using namespace std;

void EvtStdHep::init()
{
    _npart = 0;
}

int EvtStdHep::getNPart()
{
    return _npart;
}

void EvtStdHep::createParticle( EvtVector4R p4, EvtVector4R x, int prntfirst,
                                int prntlast, int id )
{
    _p4[_npart] = p4;
    _x[_npart] = x;
    _prntfirst[_npart] = prntfirst;
    _prntlast[_npart] = prntlast;
    _daugfirst[_npart] = -1;
    _dauglast[_npart] = -1;
    _id[_npart] = id;
    _istat[_npart] = 1;

    //we also need to fix up the parents pointer to the daughter!

    if ( prntfirst >= 0 ) {
        int i;
        for ( i = prntfirst; i <= prntlast; i++ ) {
            _istat[i] = 2;
            if ( _daugfirst[i] == -1 )
                _daugfirst[i] = _npart;
            if ( _dauglast[i] < _npart )
                _dauglast[i] = _npart;
        }
    }

    _npart++;
}

void EvtStdHep::translate( EvtVector4R d )
{
    int i;
    for ( i = 0; i < _npart; i++ ) {
        _x[i] += d;
    }
}

/*
ostream& operator<<(ostream& s, const EvtStdHep& stdhep){

  int w=s.width();
  int p=s.precision();
  std::ios::fmtflags f=s.flags();


  s <<endl;
  s << "  N      Id Ist   M1   M2   DF   DL      px      py      pz       E       t       x       y       z"<<endl;
  int i;
  for(i=0;i<stdhep._npart;i++){
    
    s.width(3);
    s<<i<<" ";
    s.width(7);
    s<<stdhep._id[i]<<" ";
    s.width(3);
    s<<stdhep._istat[i]<<" ";
    s.width(4);
    s<<stdhep._prntfirst[i]<<" ";
    s.width(4);
    s<<stdhep._prntlast[i]<<" ";
    s.width(4);
    s<<stdhep._daugfirst[i]<<" ";
    s.width(4);
    s<<stdhep._dauglast[i]<<" ";
    s.width(7);
    s.precision(4);
    s<<setiosflags( ios::right|ios::fixed );
    s<<stdhep._p4[i].get(1)<<" ";
    s.width(7);
    s.precision(4);
    s<<setiosflags( ios::right|ios::fixed );
    s<<stdhep._p4[i].get(2)<<" ";
    s.width(7);
    s.precision(4);
    s<<setiosflags( ios::right|ios::fixed );
    s<<stdhep._p4[i].get(3)<<" ";
    s.width(7);
    s.precision(4);
    s<<setiosflags( ios::right|ios::fixed );
    s<<stdhep._p4[i].get(0)<<" ";
    s.width(7);
    s.precision(4);
    s<<setiosflags( ios::right|ios::fixed );
    s<<stdhep._x[i].get(0)<<" ";
    s.width(7);
    s.precision(4);
    s<<setiosflags( ios::right|ios::fixed );
    s<<stdhep._x[i].get(1)<<" ";
    s.width(7);
    s.precision(4);
    s<<setiosflags( ios::right|ios::fixed );
    s<<stdhep._x[i].get(2)<<" ";
    s.width(7);
    s.precision(4);
    s<<setiosflags( ios::right|ios::fixed );
    s<<stdhep._x[i].get(3)<<endl;
    s.width(0);
  }
  
  s<<endl;

  s.width(w);
  s.precision(p);
  s.flags((std::ios::fmtflags)f);
  
  return s;

}  

*/
