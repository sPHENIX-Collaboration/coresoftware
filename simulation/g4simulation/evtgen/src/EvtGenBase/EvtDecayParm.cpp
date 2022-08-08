
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

#include "EvtGenBase/EvtDecayParm.hh"

#include "EvtGenBase/EvtPatches.hh"

#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
using std::fstream;

void EvtDecayParm::init( fcnPtr pfcn, int ndaug, int* daugs, int narg,
                         double* args, std::string name )
{
    int i;

    itsfcn = pfcn;
    itsndaug = ndaug;
    itsnarg = narg;

    itsdaugs = new int[itsndaug];
    for ( i = 0; i < itsndaug; i++ ) {
        itsdaugs[i] = daugs[i];
    }
    itsargs = new double[itsnarg];
    for ( i = 0; i < itsnarg; i++ ) {
        itsargs[i] = args[i];
    }
    modelname = name;
}

EvtDecayParm::EvtDecayParm()
{
    itsfcn = 0;
    itsndaug = 0;
    itsnarg = 0;
    itsdaugs = 0;
    itsargs = 0;

    modelname = "**********";
}

EvtDecayParm::~EvtDecayParm()
{
    if ( itsdaugs != 0 ) {
        delete[] itsdaugs;
    }

    if ( itsargs != 0 ) {
        delete[] itsargs;
    }
}
