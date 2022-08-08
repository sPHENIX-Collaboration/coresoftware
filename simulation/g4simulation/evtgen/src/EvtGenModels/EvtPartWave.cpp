
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

#include "EvtGenModels/EvtPartWave.hh"

#include "EvtGenBase/EvtCGCoefSingle.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <algorithm>
#include <stdlib.h>
#include <string>
using std::endl;

std::string EvtPartWave::getName()
{
    return "PARTWAVE";
}

EvtDecayBase* EvtPartWave::clone()
{
    return new EvtPartWave;
}

void EvtPartWave::init()
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
    int* _lambdaA2 = new int[_nA];
    int* _lambdaB2 = new int[_nB];
    int* _lambdaC2 = new int[_nC];

    EvtComplexPtr* _HBC = new EvtComplexPtr[_nB];
    int ib, ic;
    for ( ib = 0; ib < _nB; ib++ ) {
        _HBC[ib] = new EvtComplex[_nC];
    }

    int i;
    //find the allowed helicities (actually 2*times the helicity!)

    fillHelicity( _lambdaA2, _nA, _JA2 );
    fillHelicity( _lambdaB2, _nB, _JB2 );
    fillHelicity( _lambdaC2, _nC, _JC2 );

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

        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Will now figure out the valid (M_LS) states:" << endl;
    }

    int Lmin = std::max( _JA2 - _JB2 - _JC2,
                         std::max( _JB2 - _JA2 - _JC2, _JC2 - _JA2 - _JB2 ) );
    if ( Lmin < 0 )
        Lmin = 0;
    //int Lmin=_JA2-_JB2-_JC2;
    int Lmax = _JA2 + _JB2 + _JC2;

    int L;

    int _nPartialWaveAmp = 0;

    int _nL[50];
    int _nS[50];

    for ( L = Lmin; L <= Lmax; L += 2 ) {
        int Smin = abs( L - _JA2 );
        if ( Smin < abs( _JB2 - _JC2 ) )
            Smin = abs( _JB2 - _JC2 );
        int Smax = L + _JA2;
        if ( Smax > abs( _JB2 + _JC2 ) )
            Smax = abs( _JB2 + _JC2 );
        int S;
        for ( S = Smin; S <= Smax; S += 2 ) {
            _nL[_nPartialWaveAmp] = L;
            _nS[_nPartialWaveAmp] = S;

            _nPartialWaveAmp++;
            if ( verbose() ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" )
                    << "M[" << L << "][" << S << "]" << endl;
            }
        }
    }

    checkNArg( _nPartialWaveAmp * 2 );

    int argcounter = 0;

    EvtComplex _M[50];

    double partampsqtot = 0.0;

    for ( i = 0; i < _nPartialWaveAmp; i++ ) {
        _M[i] = getArg( argcounter ) *
                exp( EvtComplex( 0.0, getArg( argcounter + 1 ) ) );
        ;
        argcounter += 2;
        partampsqtot += abs2( _M[i] );
        if ( verbose() ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" )
                << "M[" << _nL[i] << "][" << _nS[i] << "]=" << _M[i] << endl;
        }
    }

    //Now calculate the helicity amplitudes

    double helampsqtot = 0.0;

    for ( ib = 0; ib < _nB; ib++ ) {
        for ( ic = 0; ic < _nC; ic++ ) {
            _HBC[ib][ic] = 0.0;
            if ( abs( _lambdaB2[ib] - _lambdaC2[ic] ) <= _JA2 ) {
                for ( i = 0; i < _nPartialWaveAmp; i++ ) {
                    int L = _nL[i];
                    int S = _nS[i];
                    int lambda2 = _lambdaB2[ib];
                    int lambda3 = _lambdaC2[ic];
                    int s1 = _JA2;
                    int s2 = _JB2;
                    int s3 = _JC2;
                    int m1 = lambda2 - lambda3;
                    EvtCGCoefSingle c1( s2, s3 );
                    EvtCGCoefSingle c2( L, S );

                    if ( verbose() ) {
                        EvtGenReport( EVTGEN_INFO, "EvtGen" )
                            << "s2,lambda2:" << s2 << " " << lambda2 << endl;
                    }
                    //fkw changes to satisfy KCC
                    double fkwTmp = ( L + 1.0 ) / ( s1 + 1.0 );

                    if ( S >= abs( m1 ) ) {
                        EvtComplex tmp = sqrt( fkwTmp ) *
                                         c1.coef( S, m1, s2, s3, lambda2,
                                                  -lambda3 ) *
                                         c2.coef( s1, m1, L, S, 0, m1 ) * _M[i];
                        _HBC[ib][ic] += tmp;
                    }
                }
                if ( verbose() ) {
                    EvtGenReport( EVTGEN_INFO, "EvtGen" )
                        << "_HBC[" << ib << "][" << ic << "]=" << _HBC[ib][ic]
                        << endl;
                }
            }
            helampsqtot += abs2( _HBC[ib][ic] );
        }
    }

    if ( fabs( helampsqtot - partampsqtot ) / ( helampsqtot + partampsqtot ) >
         1e-6 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "In EvtPartWave for decay " << EvtPDL::name( getParentId() )
            << " -> " << EvtPDL::name( getDaug( 0 ) ) << " "
            << EvtPDL::name( getDaug( 1 ) ) << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "With arguments: " << std::endl;
        for ( i = 0; i * 2 < getNArg(); i++ ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "M(" << _nL[i] << "," << _nS[i]
                << ")="
                //				 <<getArg(2*i)<<" "<<getArg(2*i+1)<<std::endl;
                << _M[i] << std::endl;
        }
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "The total probability in the partwave basis is: " << partampsqtot
            << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "The total probability in the helamp basis is: " << helampsqtot
            << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Most likely this is because the specified partwave amplitudes "
            << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "project onto unphysical helicities of photons or neutrinos. "
            << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Seriously consider if your specified amplitudes are correct. "
            << std::endl;
    }

    _evalHelAmp = std::make_unique<EvtEvalHelAmp>( getParentId(), getDaug( 0 ),
                                                   getDaug( 1 ), _HBC );
}

void EvtPartWave::initProbMax()
{
    double maxprob = _evalHelAmp->probMax();

    if ( verbose() ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "Calculated probmax" << maxprob << endl;
    }

    setProbMax( maxprob );
}

void EvtPartWave::decay( EvtParticle* p )
{
    //first generate simple phase space
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    _evalHelAmp->evalAmp( p, _amp2 );

    return;
}

void EvtPartWave::fillHelicity( int* lambda2, int n, int J2 )
{
    int i;

    //photon is special case!
    if ( n == 2 && J2 == 2 ) {
        lambda2[0] = 2;
        lambda2[1] = -2;
        return;
    }

    assert( n == J2 + 1 );

    for ( i = 0; i < n; i++ ) {
        lambda2[i] = n - i * 2 - 1;
    }

    return;
}
