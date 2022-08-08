
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

#include "EvtGenModels/EvtDMix.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"

#include <stdlib.h>
#include <string>

std::string EvtDMix::getName()
{
    return "DMIX";
}

EvtDecayBase* EvtDMix::clone()
{
    return new EvtDMix;
}

void EvtDMix::init()
{
    // check that there are 0 arguments
    checkNArg( 3 );
    _rd = getArg( 0 );
    _xpr = getArg( 1 );
    _ypr = getArg( 2 );
}

void EvtDMix::initProbMax()
{
    noProbMax();
}

void EvtDMix::decay( EvtParticle* p )
{
    //unneeded - lange - may13-02
    //if ( p->getNDaug() != 0 ) {
    //Will end up here because maxrate multiplies by 1.2
    //  EvtGenReport(EVTGEN_DEBUG,"EvtGen") << "In EvtDMix: has "
    //			   <<" daugthers should not be here!"<<endl;
    //  return;
    //}

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    double ctau = EvtPDL::getctau( p->getId() );
    if ( ctau == 0. )
        return;

    double pdf, random, gt, weight;

    double maxPdf = _rd + sqrt( _rd ) * _ypr * 50. +
                    2500.0 * ( _xpr * _xpr + _ypr * _ypr ) / 4.0;
    bool keepGoing = true;
    while ( keepGoing ) {
        random = EvtRandom::Flat();
        gt = -log( random );
        weight = random;
        pdf = _rd + sqrt( _rd ) * _ypr * gt +
              gt * gt * ( _xpr * _xpr + _ypr * _ypr ) / 4.0;
        pdf *= exp( -1.0 * gt );
        pdf /= weight;
        if ( pdf > maxPdf )
            std::cout << pdf << " " << weight << " " << maxPdf << " " << gt
                      << std::endl;
        if ( pdf > maxPdf * EvtRandom::Flat() )
            keepGoing = false;
    }

    p->setLifetime( gt * ctau );

    return;
}
