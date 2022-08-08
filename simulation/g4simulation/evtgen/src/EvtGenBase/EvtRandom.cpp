
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

#include "EvtGenBase/EvtRandom.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtReport.hh"

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using std::endl;

EvtRandomEngine* EvtRandom::_randomEngine = 0;

void EvtRandom::setRandomEngine( EvtRandomEngine* randomEngine )
{
    _randomEngine = randomEngine;
}

double EvtRandom::random()
{
    if ( _randomEngine == 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "No random engine available in "
            << "EvtRandom::random()." << endl;
        ::abort();
    }

    return _randomEngine->random();
}

// Random number routine to generate numbers between
// min and max.  By djl on July 27, 1995.
double EvtRandom::Flat( double min, double max )
{
    if ( min > max ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "min>max in EvtRandom::Flat(" << min << "," << max << ")" << endl;
        ::abort();
    }

    return EvtRandom::random() * ( max - min ) + min;
}

double EvtRandom::Flat( double max )
{
    return max * EvtRandom::random();
}

double EvtRandom::Flat()
{
    return EvtRandom::random();
}

double EvtRandom::Gaussian()
{
    double x = EvtRandom::random();
    double y = EvtRandom::random();

    return cos( x * EvtConst::twoPi ) * sqrt( -2.0 * log( 1 - y ) );
}
