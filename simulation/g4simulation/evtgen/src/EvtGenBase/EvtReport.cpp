
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

#include "EvtGenBase/EvtReport.hh"

#include "EvtGenBase/EvtPatches.hh"
using std::cerr;
using std::cout;
using std::endl;
using std::ostream;

//
// constants, enums and typedefs
//

ostream& EvtGenReport( EvtGenSeverity severity, const char* facility )
{
    int printNoFacility = 1;

    if ( ( facility == 0 ) && ( printNoFacility == 1 ) ) {
        cout << "There is no `facility' implemented in `report'" << endl;
        printNoFacility = 0;
    }
    if ( severity < EVTGEN_WARNING ) {
        if ( facility[0] != 0 ) {
            cerr << facility << ":";
        }
        return ( cerr );
    }
    if ( facility[0] != 0 ) {
        cout << facility << ":";
    }
    return cout;
}
