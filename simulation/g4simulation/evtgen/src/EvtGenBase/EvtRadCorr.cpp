
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

#include "EvtGenBase/EvtRadCorr.hh"

#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <iostream>
#include <stdlib.h>
using std::endl;

EvtAbsRadCorr* EvtRadCorr::_fsrEngine = 0;
bool EvtRadCorr::_alwaysRadCorr = false;
bool EvtRadCorr::_neverRadCorr = false;

EvtRadCorr::EvtRadCorr()
{
    _fsrEngine = 0;
    _alwaysRadCorr = false;
    _neverRadCorr = false;
}

EvtRadCorr::~EvtRadCorr()
{
    if ( _fsrEngine )
        delete _fsrEngine;
    _fsrEngine = 0;
}

void EvtRadCorr::setRadCorrEngine( EvtAbsRadCorr* fsrEngine )
{
    _fsrEngine = fsrEngine;
}

void EvtRadCorr::doRadCorr( EvtParticle* p )
{
    if ( _fsrEngine == 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "No RadCorr model available in "
            << "EvtRadCorr::doRadCorr()." << endl;
        ::abort();
    }

    if ( !_neverRadCorr )
        _fsrEngine->doRadCorr( p );
    return;
}

bool EvtRadCorr::alwaysRadCorr()
{
    return _alwaysRadCorr;
}
bool EvtRadCorr::neverRadCorr()
{
    return _neverRadCorr;
}

void EvtRadCorr::setAlwaysRadCorr()
{
    _alwaysRadCorr = true;
    _neverRadCorr = false;
}
void EvtRadCorr::setNeverRadCorr()
{
    _alwaysRadCorr = false;
    _neverRadCorr = true;
}
void EvtRadCorr::setNormalRadCorr()
{
    _alwaysRadCorr = false;
    _neverRadCorr = false;
}
