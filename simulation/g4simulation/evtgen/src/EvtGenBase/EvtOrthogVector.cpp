
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

#include "EvtGenBase/EvtOrthogVector.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>
using std::fstream;

EvtOrthogVector::EvtOrthogVector( int n, std::vector<double>* vectors )
{
    _dimen = n;
    _holder.resize( n );

    std::vector<int> temp;

    int i;
    for ( i = 0; i < n; i++ ) {
        _orthogVector.push_back( 0. );
        temp.push_back( i );
    }

    findOrthog( _dimen, temp, vectors );
}

void EvtOrthogVector::findOrthog( int dim, std::vector<int> invect,
                                  std::vector<double>* vectors )
{
    if ( dim == 2 ) {
        _holder[0] = invect[0];
        _holder[1] = invect[1];
        int sign = findEvenOddSwaps();
        {
            double addition = 1;
            int i;
            for ( i = 1; i < _dimen; i++ ) {
                addition *= vectors[i - 1][_holder[i]];
            }
            addition *= sign;
            _orthogVector[_holder[0]] += addition;
        }

        _holder[0] = invect[1];
        _holder[1] = invect[0];

        {
            double addition = 1;
            int i;
            for ( i = 1; i < _dimen; i++ ) {
                addition *= vectors[i - 1][_holder[i]];
            }
            addition *= sign;
            _orthogVector[_holder[0]] -= addition;
        }

        return;
    } else {
        std::vector<int> temp( ( 2 * dim ) );

        int i;
        for ( i = 0; i < dim; i++ )
            temp[i] = invect[i];
        for ( i = 0; i < dim; i++ )
            temp[i + dim] = invect[i];

        for ( i = 0; i < dim; i++ ) {
            _holder[dim - 1] = temp[dim - 1 + i];
            std::vector<int> tempDim( ( dim - 1 ) );

            int j;
            for ( j = 0; j < ( dim - 1 ); j++ )
                tempDim[j] = temp[j + i];
            findOrthog( dim - 1, tempDim, vectors );
        }
    }

    return;
}

int EvtOrthogVector::findEvenOddSwaps()
{
    std::vector<int> temp( _dimen );

    int i, j, nSwap;
    for ( i = 0; i < _dimen; i++ )
        temp[i] = _holder[i];

    nSwap = 0;
    for ( i = 0; i < ( _dimen - 1 ); i++ ) {
        for ( j = i + 1; j < _dimen; j++ ) {
            if ( temp[i] > temp[j] ) {
                int duh = temp[j];
                temp[j] = temp[i];
                temp[i] = duh;
                nSwap += 1;
            }
        }
    }
    nSwap -= ( nSwap / 2 ) * 2;

    if ( nSwap )
        return -1;

    return 1;
}
