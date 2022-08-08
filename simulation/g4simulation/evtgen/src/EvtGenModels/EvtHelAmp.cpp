
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

#include "EvtGenModels/EvtHelAmp.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtEvalHelAmp.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <string>
#include <vector>
using std::endl;

std::string EvtHelAmp::getName()
{
    return "HELAMP";
}

EvtDecayBase* EvtHelAmp::clone()
{
    return new EvtHelAmp;
}

void EvtHelAmp::init()
{
    checkNDaug( 2 );

    //find out how many states each particle have
    int _nA = EvtSpinType::getSpinStates( EvtPDL::getSpinType( getParentId() ) );
    int _nB = EvtSpinType::getSpinStates( EvtPDL::getSpinType( getDaug( 0 ) ) );
    int _nC = EvtSpinType::getSpinStates( EvtPDL::getSpinType( getDaug( 1 ) ) );

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "_nA,_nB,_nC:" << _nA << "," << _nB << "," << _nC << endl;
    }

    //find out what 2 times the spin is
    int _JA2 = EvtSpinType::getSpin2( EvtPDL::getSpinType( getParentId() ) );
    int _JB2 = EvtSpinType::getSpin2( EvtPDL::getSpinType( getDaug( 0 ) ) );
    int _JC2 = EvtSpinType::getSpin2( EvtPDL::getSpinType( getDaug( 1 ) ) );

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "_JA2,_JB2,_JC2:" << _JA2 << "," << _JB2 << "," << _JC2 << endl;
    }

    //allocate memory
    std::vector<int> _lambdaA2( _nA );
    std::vector<int> _lambdaB2( _nB );
    std::vector<int> _lambdaC2( _nC );

    EvtComplexPtr* _HBC = new EvtComplexPtr[_nB];
    for ( int ib = 0; ib < _nB; ib++ ) {
        _HBC[ib] = new EvtComplex[_nC];
    }

    int i;
    //find the allowed helicities (actually 2*times the helicity!)

    fillHelicity( _lambdaA2.data(), _nA, _JA2, getParentId() );
    fillHelicity( _lambdaB2.data(), _nB, _JB2, getDaug( 0 ) );
    fillHelicity( _lambdaC2.data(), _nC, _JC2, getDaug( 1 ) );

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Helicity states of particle A:" << endl;
        for ( i = 0; i < _nA; i++ ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << _lambdaA2[i] << endl;
        }

        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Helicity states of particle B:" << endl;
        for ( i = 0; i < _nB; i++ ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << _lambdaB2[i] << endl;
        }

        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Helicity states of particle C:" << endl;
        for ( i = 0; i < _nC; i++ ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << _lambdaC2[i] << endl;
        }
    }

    //now read in the helicity amplitudes

    int argcounter = 0;

    for ( int ib = 0; ib < _nB; ib++ ) {
        for ( int ic = 0; ic < _nC; ic++ ) {
            _HBC[ib][ic] = 0.0;
            if ( abs( _lambdaB2[ib] - _lambdaC2[ic] ) <= _JA2 )
                argcounter += 2;
        }
    }

    checkNArg( argcounter );

    argcounter = 0;

    for ( int ib = 0; ib < _nB; ib++ ) {
        for ( int ic = 0; ic < _nC; ic++ ) {
            if ( abs( _lambdaB2[ib] - _lambdaC2[ic] ) <= _JA2 ) {
                _HBC[ib][ic] = getArg( argcounter ) *
                               exp( EvtComplex( 0.0, getArg( argcounter + 1 ) ) );
                ;
                argcounter += 2;
                if ( verbose() ) {
                    EvtGenReport( EVTGEN_INFO, "EvtGen" )
                        << "_HBC[" << ib << "][" << ic << "]=" << _HBC[ib][ic]
                        << endl;
                }
            }
        }
    }

    _evalHelAmp = std::make_unique<EvtEvalHelAmp>( getParentId(), getDaug( 0 ),
                                                   getDaug( 1 ), _HBC );

    // Note: these are not class data members but local variables.
    for ( int ib = 0; ib < _nB; ib++ ) {
        delete[] _HBC[ib];
    }
    delete[] _HBC;    // _HBC is copied in ctor of EvtEvalHelAmp above.
}

void EvtHelAmp::initProbMax()
{
    double maxprob = _evalHelAmp->probMax();

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Calculated probmax" << maxprob << endl;
    }

    setProbMax( maxprob );
}

void EvtHelAmp::decay( EvtParticle* p )
{
    //first generate simple phase space
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    _evalHelAmp->evalAmp( p, _amp2 );
}

void EvtHelAmp::fillHelicity( int* lambda2, int n, int J2, EvtId id )
{
    int i;

    //photon is special case!
    if ( n == 2 && J2 == 2 ) {
        lambda2[0] = 2;
        lambda2[1] = -2;
        return;
    }

    //and so is the neutrino!
    if ( n == 1 && J2 == 1 ) {
        if ( EvtPDL::getStdHep( id ) > 0 ) {
            //particle i.e. lefthanded
            lambda2[0] = -1;
        } else {
            //anti particle i.e. righthanded
            lambda2[0] = 1;
        }
        return;
    }

    assert( n == J2 + 1 );

    for ( i = 0; i < n; i++ ) {
        lambda2[i] = n - i * 2 - 1;
    }

    return;
}
