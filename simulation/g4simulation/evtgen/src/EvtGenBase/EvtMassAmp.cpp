
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

#include "EvtGenBase/EvtMassAmp.hh"

#include "EvtGenBase/EvtPatches.hh"

EvtMassAmp::EvtMassAmp( const EvtPropBreitWignerRel& prop,
                        const EvtTwoBodyVertex& vd ) :
    EvtAmplitude<EvtPoint1D>(),
    _prop( prop ),
    _vd( vd ),
    _useBirthFact( false ),
    _useDeathFact( false ),
    _useBirthFactFF( false ),
    _useDeathFactFF( false )
{
}

EvtMassAmp::EvtMassAmp( const EvtMassAmp& other ) :
    EvtAmplitude<EvtPoint1D>( other ),
    _prop( other._prop ),
    _vd( other._vd ),
    _vb( other._vb ? new EvtTwoBodyVertex( *other._vb ) : nullptr ),
    _useBirthFact( other._useBirthFact ),
    _useDeathFact( other._useDeathFact ),
    _useBirthFactFF( other._useBirthFactFF ),
    _useDeathFactFF( other._useDeathFactFF )
{
}

EvtMassAmp& EvtMassAmp::operator=( const EvtMassAmp& other )
{
    EvtAmplitude<EvtPoint1D>::operator=( other );
    _prop = other._prop;
    _vd = other._vd;
    _vb.reset( other._vb ? new EvtTwoBodyVertex( *other._vb ) : nullptr );
    _useBirthFact = other._useBirthFact;
    _useDeathFact = other._useDeathFact;
    _useBirthFactFF = other._useBirthFactFF;
    _useDeathFactFF = other._useDeathFactFF;
    return *this;
}

EvtComplex EvtMassAmp::amplitude( const EvtPoint1D& p ) const
{
    // Modified vertex

    double m = p.value();
    // keep things from crashing..

    if ( m < ( _vd.mA() + _vd.mB() ) )
        return EvtComplex( 0., 0. );

    EvtTwoBodyKine vd( _vd.mA(), _vd.mB(), m );

    // Compute mass-dependent width for relativistic propagator

    EvtPropBreitWignerRel bw( _prop.m0(), _prop.g0() * _vd.widthFactor( vd ) );
    EvtComplex amp = bw.evaluate( m );

    // Birth vertex factors

    if ( _useBirthFact ) {
        assert( _vb );
        if ( ( m + _vb->mB() ) < _vb->mAB() ) {
            EvtTwoBodyKine vb( m, _vb->mB(), _vb->mAB() );
            amp *= _vb->phaseSpaceFactor( vb, EvtTwoBodyKine::AB );
            amp *= sqrt( ( vb.p() / _vb->pD() ) );

            if ( _useBirthFactFF ) {
                assert( _vb );
                amp *= _vb->formFactor( vb );
            }
        } else {
            if ( _vb->L() != 0 )
                amp = 0.;
        }
    }

    // Decay vertex factors

    if ( _useDeathFact ) {
        amp *= _vd.phaseSpaceFactor( vd, EvtTwoBodyKine::AB );
        amp *= sqrt( ( vd.p() / _vd.pD() ) );
    }
    if ( _useDeathFactFF )
        amp *= _vd.formFactor( vd );

    return amp;
}
