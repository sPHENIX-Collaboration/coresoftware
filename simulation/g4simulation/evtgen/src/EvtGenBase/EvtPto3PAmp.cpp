
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

#include "EvtGenBase/EvtPto3PAmp.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzCoord.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtdFunction.hh"

#include <assert.h>
#include <iostream>
#include <math.h>

using EvtCyclic3::Index;
using EvtCyclic3::Pair;
using std::endl;

EvtPto3PAmp::EvtPto3PAmp( EvtDalitzPlot dp, Pair pairAng, Pair pairRes,
                          EvtSpinType::spintype spin, const EvtPropagator& prop,
                          NumType typeN ) :
    EvtAmplitude<EvtDalitzPoint>(),
    _pairAng( pairAng ),
    _pairRes( pairRes ),
    _spin( spin ),
    _typeN( typeN ),
    _prop( (EvtPropagator*)prop.clone() ),
    _g0( prop.g0() ),
    _min( 0 ),
    _max( 0 ),
    _vb( prop.m0(), dp.m( EvtCyclic3::other( pairRes ) ), dp.bigM(), spin ),
    _vd( dp.m( EvtCyclic3::first( pairRes ) ),
         dp.m( EvtCyclic3::second( pairRes ) ), prop.m0(), spin )
{
}

EvtPto3PAmp::EvtPto3PAmp( const EvtPto3PAmp& other ) :
    EvtAmplitude<EvtDalitzPoint>( other ),
    _pairAng( other._pairAng ),
    _pairRes( other._pairRes ),
    _spin( other._spin ),
    _typeN( other._typeN ),
    _prop( ( other._prop ) ? (EvtPropagator*)other._prop->clone() : 0 ),
    _g0( other._g0 ),
    _min( other._min ),
    _max( other._max ),
    _vb( other._vb ),
    _vd( other._vd )
{
}

EvtPto3PAmp::~EvtPto3PAmp()
{
    if ( _prop )
        delete _prop;
}

void EvtPto3PAmp::set_fd( double R )
{
    _vd.set_f( R );
}

void EvtPto3PAmp::set_fb( double R )
{
    _vb.set_f( R );
}

EvtComplex EvtPto3PAmp::amplitude( const EvtDalitzPoint& x ) const
{
    EvtComplex amp( 1.0, 0.0 );

    double m = sqrt( x.q( _pairRes ) );

    if ( ( _max > 0 && m > _max ) || ( _min > 0 && m < _min ) )
        return EvtComplex( 0.0, 0.0 );

    EvtTwoBodyKine vd( x.m( EvtCyclic3::first( _pairRes ) ),
                       x.m( EvtCyclic3::second( _pairRes ) ), m );
    EvtTwoBodyKine vb( m, x.m( EvtCyclic3::other( _pairRes ) ), x.bigM() );

    // Compute mass-dependent width for relativistic propagators

    if ( _typeN != NBW && _typeN != FLATTE ) {
        _prop->set_g0( _g0 * _vd.widthFactor( vd ) );
    }

    // Compute propagator

    amp *= evalPropagator( m );

    // Compute form-factors

    amp *= _vd.formFactor( vd );
    amp *= _vb.formFactor( vb );

    amp *= numerator( x );

    return amp;
}

EvtComplex EvtPto3PAmp::numerator( const EvtDalitzPoint& x ) const
{
    EvtComplex ret( 0., 0. );
    double m = sqrt( x.q( _pairRes ) );
    EvtTwoBodyKine vd( x.m( EvtCyclic3::first( _pairRes ) ),
                       x.m( EvtCyclic3::second( _pairRes ) ), m );
    EvtTwoBodyKine vb( m, x.m( EvtCyclic3::other( _pairRes ) ), x.bigM() );

    // Non-relativistic Breit-Wigner

    if ( NBW == _typeN ) {
        ret = angDep( x );
    }

    // Standard relativistic Zemach propagator

    else if ( RBW_ZEMACH == _typeN ) {
        ret = _vd.phaseSpaceFactor( vd, EvtTwoBodyKine::AB ) * angDep( x );
    }

    // Kuehn-Santamaria normalization:

    else if ( RBW_KUEHN == _typeN ) {
        ret = _prop->m0() * _prop->m0() * angDep( x );
    }

    // CLEO amplitude is not factorizable
    //
    // The CLEO amplitude numerator is proportional to:
    //
    // m2_AC - m2_BC + (m2_D - m2_C)(m2_B - m2_A)/m2_0
    //
    // m2_AC = (eA + eC)^2 + (P - P_C cosTh(BC))^2
    // m2_BC = (eB + eC)^2 + (P + P_C cosTh(BC))^2
    //
    // The first term m2_AB-m2_BC is therefore a p-wave term
    // - 4PP_C cosTh(BC)
    // The second term is an s-wave, the amplitude
    // does not factorize!
    //
    // The first term is just Zemach. However, the sign is flipped!
    // Let's consistently use the convention in which the amplitude
    // is proportional to +cosTh(BC). In the CLEO expressions, I will
    // therefore exchange AB to get rid of the sign flip.

    if ( RBW_CLEO == _typeN || FLATTE == _typeN || GS == _typeN ) {
        Index iA = EvtCyclic3::other( _pairAng );               // A = other(BC)
        Index iB = EvtCyclic3::common( _pairRes, _pairAng );    // B = common(AB,BC)
        Index iC = EvtCyclic3::other( _pairRes );               // C = other(AB)

        double M = x.bigM();
        double mA = x.m( iA );
        double mB = x.m( iB );
        double mC = x.m( iC );
        double qAB = x.q( EvtCyclic3::combine( iA, iB ) );
        double qBC = x.q( EvtCyclic3::combine( iB, iC ) );
        double qCA = x.q( EvtCyclic3::combine( iC, iA ) );

        //double m0 = _prop->m0();

        if ( _spin == EvtSpinType::SCALAR )
            ret = EvtComplex( 1., 0. );
        else if ( _spin == EvtSpinType::VECTOR ) {
            //ret = qCA - qBC - (M*M - mC*mC)*(mA*mA - mB*mB)/m0/m0;
            ret = qCA - qBC - ( M * M - mC * mC ) * ( mA * mA - mB * mB ) / qAB;
        } else if ( _spin == EvtSpinType::TENSOR ) {
            //double x1 = qBC - qCA + (M*M - mC*mC)*(mA*mA - mB*mB)/m0/m0;
            double x1 = qBC - qCA +
                        ( M * M - mC * mC ) * ( mA * mA - mB * mB ) / qAB;
            double x2 = M * M - mC * mC;
            //double x3 = qAB - 2*M*M - 2*mC*mC + x2*x2/m0/m0;
            double x3 = qAB - 2 * M * M - 2 * mC * mC + x2 * x2 / qAB;
            double x4 = mB * mB - mA * mA;
            //double x5 = qAB - 2*mB*mB - 2*mA*mA + x4*x4/m0/m0;
            double x5 = qAB - 2 * mB * mB - 2 * mA * mA + x4 * x4 / qAB;
            ret = ( x1 * x1 - 1. / 3. * x3 * x5 );
        } else
            assert( 0 );
    }

    return ret;
}

double EvtPto3PAmp::angDep( const EvtDalitzPoint& x ) const
{
    // Angular dependece for factorizable amplitudes
    // unphysical cosines indicate we are in big trouble

    double cosTh = x.cosTh( _pairAng, _pairRes );
    if ( fabs( cosTh ) > 1. ) {
        EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "cosTh " << cosTh << endl;
        assert( 0 );
    }

    // in units of half-spin

    return EvtdFunction::d( EvtSpinType::getSpin2( _spin ), 2 * 0, 2 * 0,
                            acos( cosTh ) );
}
