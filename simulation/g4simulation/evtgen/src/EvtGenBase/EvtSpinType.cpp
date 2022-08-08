
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

#include "EvtGenBase/EvtSpinType.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

int EvtSpinType::getSpin2( spintype stype )
{
    switch ( stype ) {
        case SCALAR:
        case STRING:
            return 0;
        case DIRAC:
        case NEUTRINO:
            return 1;
        case VECTOR:
        case PHOTON:
            return 2;
        case RARITASCHWINGER:
            return 3;
        case TENSOR:
            return 4;
        case SPIN5HALF:
            return 5;
        case SPIN3:
            return 6;
        case SPIN7HALF:
            return 7;
        case SPIN4:
            return 8;
        default:
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Unknown spintype in EvtSpinType!" << std::endl;
            return 0;
    }
}

int EvtSpinType::getSpinStates( spintype stype )
{
    switch ( stype ) {
        case SCALAR:
        case STRING:
        case NEUTRINO:
            return 1;
        case DIRAC:
        case PHOTON:
            return 2;
        case VECTOR:
            return 3;
        case RARITASCHWINGER:
            return 4;
        case TENSOR:
            return 5;
        case SPIN5HALF:
            return 6;
        case SPIN3:
            return 7;
        case SPIN7HALF:
            return 8;
        case SPIN4:
            return 9;
        default:
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Unknown spintype in EvtSpinType!" << std::endl;
            return 0;
    }
}
