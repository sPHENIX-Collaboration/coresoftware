
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

#include "EvtGenBase/EvtAmpIndex.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <vector>
using std::vector;

EvtAmpIndex::EvtAmpIndex( std::vector<int> ind ) :
    _ind( ind ), _size( ind.size() ), _state( ind.size() ), _nstate( ind.size() )
{
    int i;

    for ( i = 0; i < _size; i++ ) {
        _state[i] = 0;
        if ( i == 0 ) {
            _nstate[i] = 1;
        } else {
            _nstate[i] = _nstate[i - 1] * _ind[i];
        }
    }
}

void EvtAmpIndex::reset()
{
    int i;
    for ( i = 0; i < _size; i++ ) {
        _state[i] = 0;
    }
}

bool EvtAmpIndex::next()
{
    int i;
    for ( i = 0; i < _size; i++ ) {
        _state[i]++;
        if ( _state[i] < _ind[i] ) {
            return true;
        } else {
            _state[i] = 0;
        }
    }
    return false;
}

int EvtAmpIndex::index()
{
    int i;
    int ind = 0;

    for ( i = 0; i < _size; i++ ) {
        ind += _state[i] * _nstate[i];
    }

    return ind;
}
