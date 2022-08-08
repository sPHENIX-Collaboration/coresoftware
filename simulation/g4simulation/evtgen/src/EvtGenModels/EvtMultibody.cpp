
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

#include "EvtGenModels/EvtMultibody.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtResonance.hh"
#include "EvtGenBase/EvtResonance2.hh"
#include "EvtGenBase/EvtdFunction.hh"

EvtMultibody::~EvtMultibody()
{
    if ( _decayTree != NULL )
        delete _decayTree;
    _decayTree = NULL;
    if ( _ilist != NULL )
        delete[] _ilist;
    _ilist = NULL;
}

std::string EvtMultibody::getName()
{
    return "D_MULTIBODY";
}

EvtDecayBase* EvtMultibody::clone()
{
    return new EvtMultibody;
}

void EvtMultibody::init()
{
    int N = getNArg();

    _decayTree = new EvtMTree( getDaugs(), getNDaug() );
    _ilist = new int[getNDaug() + 1];

    for ( int i = 0; i < N - 1; ++i ) {
        if ( getArgStr( i ) == "RESONANCE" ) {
            _decayTree->addtree( getArgStr( ++i ) );
        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Syntax error at " << getArgStr( i ) << std::endl;
            ::abort();
        }
    }
}

// Set the maximum probability amplitude - if function is left blank then the
// program will search for it.  This however is not deterministic and therefore
// in the release cannot be in place.
void EvtMultibody::initProbMax()
{
    // setProbMax(1.0);
}

void EvtMultibody::decay( EvtParticle* p )
{
    // Initialize the phase space before doing anything else!
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtSpinAmp amp = _decayTree->amplitude( p );

    vector<int> index = amp.iterallowedinit();
    vector<unsigned int> spins = amp.dims();

    do {
        for ( size_t i = 0; i < index.size(); ++i ) {
            _ilist[i] = index[i] + spins[i];
        }

        vertex( _ilist, amp( index ) );
    } while ( amp.iterateallowed( index ) );
}
