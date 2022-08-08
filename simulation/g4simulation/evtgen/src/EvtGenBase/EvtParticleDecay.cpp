
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

#include "EvtGenBase/EvtParticleDecay.hh"

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"

#include <ctype.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
using std::fstream;
void EvtParticleDecay::printSummary()
{
    if ( _decay != 0 )
        _decay->printSummary();
}

void EvtParticleDecay::chargeConj( EvtParticleDecay* decay )
{
    _brfrsum = decay->_brfrsum;
    _massmin = decay->_massmin;

    _decay = decay->_decay->clone();

    int ndaug = decay->_decay->getNDaug();
    int narg = decay->_decay->getNArg();
    double brfr = decay->_decay->getBranchingFraction();
    std::string name = decay->_decay->getName();
    EvtId ipar = EvtPDL::chargeConj( decay->_decay->getParentId() );
    int i;
    EvtId* daug = new EvtId[ndaug];
    for ( i = 0; i < ndaug; i++ ) {
        daug[i] = EvtPDL::chargeConj( decay->_decay->getDaug( i ) );
    }
    //Had to add 1 to make sure the vector is not empty!
    std::vector<std::string> args;
    for ( i = 0; i < narg; i++ ) {
        args.push_back( decay->_decay->getArgStr( i ) );
    }

    _decay->saveDecayInfo( ipar, ndaug, daug, narg, args, name, brfr );

    if ( decay->_decay->getPHOTOS() )
        _decay->setPHOTOS();

    delete[] daug;
}
