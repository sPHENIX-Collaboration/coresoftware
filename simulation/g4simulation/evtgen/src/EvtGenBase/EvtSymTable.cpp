
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

#include "EvtGenBase/EvtSymTable.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
using std::endl;
using std::fstream;

std::map<std::string, std::string> EvtSymTable::_symMap;

EvtSymTable::EvtSymTable()
{
}

void EvtSymTable::define( const std::string& symname, std::string d )
{
    if ( _symMap.find( symname ) != _symMap.end() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Symbol:" << symname.c_str()
            << " redefined, old value:" << _symMap[symname].c_str()
            << " new value:" << d.c_str() << endl;
        _symMap[symname] = d;
        return;
    }

    _symMap[symname] = d;
    return;
}

std::string EvtSymTable::get( const std::string& symname, int& ierr )
{
    ierr = 0;

    if ( _symMap.find( symname ) != _symMap.end() )
        return _symMap[symname];

    // If no matching symbol found just return the string

    return symname;
}
