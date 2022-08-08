
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

#include "EvtGenBase/EvtMTRandomEngine.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <iostream>

EvtMTRandomEngine::EvtMTRandomEngine( unsigned int seed ) :
    engine_( seed ), distribution_( URDist( 0.0, 1.0 ) )
{
    EvtGenReport( EVTGEN_INFO, "EvtMTRandomEngine" )
        << "Mersenne-Twister random number generator with seed = " << seed
        << std::endl;
}

double EvtMTRandomEngine::random()
{
    return distribution_( engine_ );
}
