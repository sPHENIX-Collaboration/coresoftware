
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

#include "EvtGenModels/EvtBtoXsgamma.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include "EvtGenModels/EvtBtoXsgammaAliGreub.hh"
#include "EvtGenModels/EvtBtoXsgammaFixedMass.hh"
#include "EvtGenModels/EvtBtoXsgammaFlatEnergy.hh"
#include "EvtGenModels/EvtBtoXsgammaKagan.hh"

#include <stdlib.h>
#include <string>
using std::endl;

std::string EvtBtoXsgamma::getName()
{
    return "BTOXSGAMMA";
}

EvtDecayBase* EvtBtoXsgamma::clone()
{
    return new EvtBtoXsgamma;
}

void EvtBtoXsgamma::init()
{
    //Arguments:
    // 0: Ali-Greub model = 1, Kagan model = 2
    //No more arguments for Ali-Greub model
    // 1:
    // 2:
    // 3:

    // check that at least one b->sg model has been selected
    if ( getNArg() == 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBtoXsgamma generator expected "
            << " at least 1 argument but found: " << getNArg() << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}

void EvtBtoXsgamma::initProbMax()
{
    noProbMax();
}

void EvtBtoXsgamma::decay( EvtParticle* p )
{
    //initialize here. -- its too damn slow otherwise.

    if ( _model == 0 ) {
        if ( getArg( 0 ) == 1 )
            _model = std::make_unique<EvtBtoXsgammaAliGreub>();
        else if ( getArg( 0 ) == 2 )
            _model = std::make_unique<EvtBtoXsgammaKagan>();
        else if ( getArg( 0 ) == 3 )
            _model = std::make_unique<EvtBtoXsgammaFixedMass>();
        else if ( getArg( 0 ) == 4 )
            _model = std::make_unique<EvtBtoXsgammaFlatEnergy>();
        else {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "No valid EvtBtoXsgamma generator model selected "
                << "Set arg(0) to 1 for Ali-Greub model or 2 for "
                << " Kagan model or 3 for a fixed mass" << endl;
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Will terminate execution!" << endl;
            ::abort();
        }
        _model->init( getNArg(), getArgs() );
    }

    //  if ( p->getNDaug() != 0 ) {
    //Will end up here because maxrate multiplies by 1.2
    //  EvtGenReport(EVTGEN_DEBUG,"EvtGen") << "In EvtBtoXsgamma: X_s daughters should not be here!"<<endl;
    //  return;
    //}

    double m_b;
    int i;
    p->makeDaughters( getNDaug(), getDaugs() );
    EvtParticle* pdaug[MAX_DAUG];

    for ( i = 0; i < getNDaug(); i++ ) {
        pdaug[i] = p->getDaug( i );
    }

    static EvtVector4R p4[MAX_DAUG];
    static double mass[MAX_DAUG];

    m_b = p->mass();

    mass[1] = EvtPDL::getMass( getDaug( 1 ) );

    int Xscode = EvtPDL::getStdHep( getDaug( 0 ) );

    mass[0] = _model->GetMass( Xscode );

    EvtGenKine::PhaseSpace( getNDaug(), mass, p4, m_b );

    for ( i = 0; i < getNDaug(); i++ ) {
        pdaug[i]->init( getDaugs()[i], p4[i] );
    }
}
