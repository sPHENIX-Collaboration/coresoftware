
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

#if !defined( TOOLBOX_FUNCTIONS_HH )
#define TOOLBOX_FUNCTIONS_HH

#if !defined( FILENAME_ONLY ) /* relative path includes */

// system include files
#include <iostream>

// user include files

#else /* filename-only includes */
#include <iostream>
#include <types.h>
#endif /* filename-only includes */
// system include files

// user include files

// forward declarations

//
// constants, enums and typedefs
//
enum EvtGenSeverity
{
    EVTGEN_EMERGENCY,    // fatal
    EVTGEN_ALERT,        // requires immediate action
    EVTGEN_CRITICAL,     // serious
    EVTGEN_ERROR,
    EVTGEN_WARNING,
    EVTGEN_NOTICE,    // "normal but significant"
    EVTGEN_INFO,      // informational
    EVTGEN_DEBUG      // debug
};

// function declaration
std::ostream& EvtGenReport( EvtGenSeverity severity, const char* facility = 0 );

// inline function definitions

#endif /* TOOLBOX_FUNCTIONS_HH */
