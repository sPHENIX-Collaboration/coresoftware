
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

#include "EvtGenBase/EvtAmp.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSpinDensity.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
using std::endl;

EvtAmp::EvtAmp()
{
    _ndaug = 0;
    _pstates = 0;
    _nontrivial = 0;
}

EvtAmp::EvtAmp( const EvtAmp& amp )
{
    int i;

    _ndaug = amp._ndaug;
    _pstates = amp._pstates;
    for ( i = 0; i < _ndaug; i++ ) {
        dstates[i] = amp.dstates[i];
        _dnontrivial[i] = amp._dnontrivial[i];
    }
    _nontrivial = amp._nontrivial;

    int namp = 1;

    for ( i = 0; i < _nontrivial; i++ ) {
        _nstate[i] = amp._nstate[i];
        namp *= _nstate[i];
    }

    for ( i = 0; i < namp; i++ ) {
        assert( i < 125 );
        _amp[i] = amp._amp[i];
    }
}

void EvtAmp::init( EvtId p, int ndaugs, EvtId* daug )
{
    setNDaug( ndaugs );
    int ichild;
    int daug_states[100], parstates;
    for ( ichild = 0; ichild < ndaugs; ichild++ ) {
        daug_states[ichild] = EvtSpinType::getSpinStates(
            EvtPDL::getSpinType( daug[ichild] ) );
    }

    parstates = EvtSpinType::getSpinStates( EvtPDL::getSpinType( p ) );

    setNState( parstates, daug_states );
}

void EvtAmp::setNDaug( int n )
{
    _ndaug = n;
}

void EvtAmp::setNState( int parent_states, int* daug_states )
{
    _nontrivial = 0;
    _pstates = parent_states;

    if ( _pstates > 1 ) {
        _nstate[_nontrivial] = _pstates;
        _nontrivial++;
    }

    int i;

    for ( i = 0; i < _ndaug; i++ ) {
        dstates[i] = daug_states[i];
        _dnontrivial[i] = -1;
        if ( daug_states[i] > 1 ) {
            _nstate[_nontrivial] = daug_states[i];
            _dnontrivial[i] = _nontrivial;
            _nontrivial++;
        }
    }

    if ( _nontrivial > 5 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Too many nontrivial states in EvtAmp!" << endl;
    }
}

void EvtAmp::setAmp( int* ind, const EvtComplex& a )
{
    int nstatepad = 1;
    int position = ind[0];

    for ( int i = 1; i < _nontrivial; i++ ) {
        nstatepad *= _nstate[i - 1];
        position += nstatepad * ind[i];
    }
    assert( position < 125 );
    _amp[position] = a;
}

const EvtComplex& EvtAmp::getAmp( int* ind ) const
{
    int nstatepad = 1;
    int position = ind[0];

    for ( int i = 1; i < _nontrivial; i++ ) {
        nstatepad *= _nstate[i - 1];
        position += nstatepad * ind[i];
    }

    return _amp[position];
}

EvtSpinDensity EvtAmp::getSpinDensity()
{
    EvtSpinDensity rho;
    rho.setDim( _pstates );

    EvtComplex temp;

    int i, j, n;

    if ( _pstates == 1 ) {
        if ( _nontrivial == 0 ) {
            rho.set( 0, 0, _amp[0] * conj( _amp[0] ) );
            return rho;
        }

        n = 1;

        temp = EvtComplex( 0.0 );

        for ( i = 0; i < _nontrivial; i++ ) {
            n *= _nstate[i];
        }

        for ( i = 0; i < n; i++ ) {
            temp += _amp[i] * conj( _amp[i] );
        }

        rho.set( 0, 0, temp );
        ;

        return rho;

    }

    else {
        for ( i = 0; i < _pstates; i++ ) {
            for ( j = 0; j < _pstates; j++ ) {
                temp = EvtComplex( 0.0 );

                int kk;

                int allloop = 1;
                for ( kk = 0; kk < _ndaug; kk++ ) {
                    allloop *= dstates[kk];
                }

                for ( kk = 0; kk < allloop; kk++ ) {
                    temp += _amp[_pstates * kk + i] *
                            conj( _amp[_pstates * kk + j] );
                }

                //        if (_nontrivial>3){
                //EvtGenReport(EVTGEN_ERROR,"EvtGen") << "Can't handle so many states in EvtAmp!"<<endl;
                //}

                rho.set( i, j, temp );
            }
        }
        return rho;
    }
}

EvtSpinDensity EvtAmp::getBackwardSpinDensity( EvtSpinDensity* rho_list )
{
    EvtSpinDensity rho;

    rho.setDim( _pstates );

    if ( _pstates == 1 ) {
        rho.set( 0, 0, EvtComplex( 1.0, 0.0 ) );
        return rho;
    }

    int k;

    EvtAmp ampprime;

    ampprime = ( *this );

    for ( k = 0; k < _ndaug; k++ ) {
        if ( dstates[k] != 1 ) {
            ampprime = ampprime.contract( _dnontrivial[k], rho_list[k + 1] );
        }
    }

    return ampprime.contract( 0, ( *this ) );
}

EvtSpinDensity EvtAmp::getForwardSpinDensity( EvtSpinDensity* rho_list, int i )
{
    EvtSpinDensity rho;

    rho.setDim( dstates[i] );

    int k;

    if ( dstates[i] == 1 ) {
        rho.set( 0, 0, EvtComplex( 1.0, 0.0 ) );

        return rho;
    }

    EvtAmp ampprime;

    ampprime = ( *this );

    if ( _pstates != 1 ) {
        ampprime = ampprime.contract( 0, rho_list[0] );
    }

    for ( k = 0; k < i; k++ ) {
        if ( dstates[k] != 1 ) {
            ampprime = ampprime.contract( _dnontrivial[k], rho_list[k + 1] );
        }
    }

    return ampprime.contract( _dnontrivial[i], ( *this ) );
}

EvtAmp EvtAmp::contract( int k, const EvtSpinDensity& rho )
{
    EvtAmp temp;

    int i, j;
    temp._ndaug = _ndaug;
    temp._pstates = _pstates;
    temp._nontrivial = _nontrivial;

    for ( i = 0; i < _ndaug; i++ ) {
        temp.dstates[i] = dstates[i];
        temp._dnontrivial[i] = _dnontrivial[i];
    }

    if ( _nontrivial == 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Should not be here EvtAmp!" << endl;
    }

    for ( i = 0; i < _nontrivial; i++ ) {
        temp._nstate[i] = _nstate[i];
    }

    EvtComplex c;

    int index[10];
    for ( i = 0; i < 10; i++ ) {
        index[i] = 0;
    }

    int allloop = 1;
    int indflag, ii;
    for ( i = 0; i < _nontrivial; i++ ) {
        allloop *= _nstate[i];
    }

    for ( i = 0; i < allloop; i++ ) {
        c = EvtComplex( 0.0 );
        int tempint = index[k];
        for ( j = 0; j < _nstate[k]; j++ ) {
            index[k] = j;
            c += rho.get( j, tempint ) * getAmp( index );
        }
        index[k] = tempint;

        temp.setAmp( index, c );

        indflag = 0;
        for ( ii = 0; ii < _nontrivial; ii++ ) {
            if ( indflag == 0 ) {
                if ( index[ii] == ( _nstate[ii] - 1 ) ) {
                    index[ii] = 0;
                } else {
                    indflag = 1;
                    index[ii] += 1;
                }
            }
        }
    }
    return temp;
}

EvtSpinDensity EvtAmp::contract( int k, const EvtAmp& amp2 )
{
    int i, j, l;

    EvtComplex temp;
    EvtSpinDensity rho;

    rho.setDim( _nstate[k] );

    int allloop = 1;
    int indflag, ii;
    for ( i = 0; i < _nontrivial; i++ ) {
        allloop *= _nstate[i];
    }

    int index[10];
    int index1[10];
    //  int l;
    for ( i = 0; i < _nstate[k]; i++ ) {
        for ( j = 0; j < _nstate[k]; j++ ) {
            if ( _nontrivial == 0 ) {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "Should not be here1 EvtAmp!" << endl;
                rho.set( 0, 0, EvtComplex( 1.0, 0.0 ) );
                return rho;
            }

            for ( ii = 0; ii < 10; ii++ ) {
                index[ii] = 0;
                index1[ii] = 0;
            }
            index[k] = i;
            index1[k] = j;

            temp = EvtComplex( 0.0 );

            for ( l = 0; l < int( allloop / _nstate[k] ); l++ ) {
                temp += getAmp( index ) * conj( amp2.getAmp( index1 ) );
                indflag = 0;
                for ( ii = 0; ii < _nontrivial; ii++ ) {
                    if ( ii != k ) {
                        if ( indflag == 0 ) {
                            if ( index[ii] == ( _nstate[ii] - 1 ) ) {
                                index[ii] = 0;
                                index1[ii] = 0;
                            } else {
                                indflag = 1;
                                index[ii] += 1;
                                index1[ii] += 1;
                            }
                        }
                    }
                }
            }
            rho.set( i, j, temp );
        }
    }

    return rho;
}

EvtAmp EvtAmp::contract( int, const EvtAmp&, const EvtAmp& )
{
    //Do we need this method?
    EvtAmp tmp;
    EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
        << "EvtAmp::contract not written yet" << endl;
    return tmp;
}

void EvtAmp::dump()
{
    int i, list[10];
    for ( i = 0; i < 10; i++ ) {
        list[i] = 0;
    }

    EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
        << "Number of daugthers:" << _ndaug << endl;
    EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
        << "Number of states of the parent:" << _pstates << endl;
    EvtGenReport( EVTGEN_DEBUG, "EvtGen" ) << "Number of states on daughters:";
    for ( i = 0; i < _ndaug; i++ ) {
        EvtGenReport( EVTGEN_DEBUG, "EvtGen" ) << dstates[i] << " ";
    }
    EvtGenReport( EVTGEN_DEBUG, "EvtGen" ) << endl;
    EvtGenReport( EVTGEN_DEBUG, "EvtGen" ) << "Nontrivial index of  daughters:";
    for ( i = 0; i < _ndaug; i++ ) {
        EvtGenReport( EVTGEN_DEBUG, "EvtGen" ) << _dnontrivial[i] << " ";
    }
    EvtGenReport( EVTGEN_DEBUG, "EvtGen" ) << endl;
    EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
        << "number of nontrivial states:" << _nontrivial << endl;
    EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
        << "Nontrivial particles number of states:";
    for ( i = 0; i < _nontrivial; i++ ) {
        EvtGenReport( EVTGEN_DEBUG, "EvtGen" ) << _nstate[i] << " ";
    }
    EvtGenReport( EVTGEN_DEBUG, "EvtGen" ) << endl;
    EvtGenReport( EVTGEN_DEBUG, "EvtGen" ) << "Amplitudes:" << endl;
    if ( _nontrivial == 0 ) {
        list[0] = 0;
        EvtGenReport( EVTGEN_DEBUG, "EvtGen" ) << getAmp( list ) << endl;
    }

    int allloop[10];
    for ( i = 0; i < 10; i++ ) {
        allloop[i] = 0;
    }

    allloop[0] = 1;
    for ( i = 0; i < _nontrivial; i++ ) {
        if ( i == 0 ) {
            allloop[i] *= _nstate[i];
        } else {
            allloop[i] = allloop[i - 1] * _nstate[i];
        }
    }
    int index = 0;
    for ( i = 0; i < allloop[_nontrivial - 1]; i++ ) {
        EvtGenReport( EVTGEN_DEBUG, "EvtGen" ) << getAmp( list ) << " ";
        if ( i == allloop[index] - 1 ) {
            index++;
            EvtGenReport( EVTGEN_DEBUG, "EvtGen" ) << endl;
        }
    }

    EvtGenReport( EVTGEN_DEBUG, "EvtGen" )
        << "-----------------------------------" << endl;
}

void EvtAmp::vertex( const EvtComplex& c )
{
    int list[1];
    list[0] = 0;
    setAmp( list, c );
}

void EvtAmp::vertex( int i, const EvtComplex& c )
{
    int list[1];
    list[0] = i;
    setAmp( list, c );
}

void EvtAmp::vertex( int i, int j, const EvtComplex& c )
{
    int list[2];
    list[0] = i;
    list[1] = j;
    setAmp( list, c );
}

void EvtAmp::vertex( int i, int j, int k, const EvtComplex& c )
{
    int list[3];
    list[0] = i;
    list[1] = j;
    list[2] = k;
    setAmp( list, c );
}

void EvtAmp::vertex( int* i1, const EvtComplex& c )
{
    setAmp( i1, c );
}

EvtAmp& EvtAmp::operator=( const EvtAmp& amp )
{
    int i;

    _ndaug = amp._ndaug;
    _pstates = amp._pstates;
    for ( i = 0; i < _ndaug; i++ ) {
        dstates[i] = amp.dstates[i];
        _dnontrivial[i] = amp._dnontrivial[i];
    }
    _nontrivial = amp._nontrivial;

    int namp = 1;

    for ( i = 0; i < _nontrivial; i++ ) {
        _nstate[i] = amp._nstate[i];
        namp *= _nstate[i];
    }

    for ( i = 0; i < namp; i++ ) {
        assert( i < 125 );
        _amp[i] = amp._amp[i];
    }

    return *this;
}
