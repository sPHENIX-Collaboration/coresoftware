
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

#include "EvtGenBase/EvtAmpSubIndex.hh"

#include "EvtGenBase/EvtAmpIndex.hh"
#include "EvtGenBase/EvtPatches.hh"

#include <vector>
using std::vector;

EvtAmpSubIndex::EvtAmpSubIndex( EvtAmpIndex* ind, std::vector<int> sub ) :
    _ind( ind ), _sub( sub ), _size( sub.size() ), _nstate( sub.size() )
{
    int i;

    for ( i = 0; i < _size; i++ ) {
        if ( i == 0 ) {
            _nstate[i] = 1;
        } else {
            _nstate[i] = _nstate[i - 1] * _ind->_ind[sub[i - 1]];
        }
    }
}

int EvtAmpSubIndex::index()
{
    int i;
    int ind = 0;

    for ( i = 0; i < _size; i++ ) {
        ind += _ind->_state[_ind->_ind[i]] * _nstate[i];
    }

    return ind;
}
