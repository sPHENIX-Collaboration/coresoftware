
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

#include "EvtGenBase/EvtSpinDensity.hh"

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
using std::endl;
using std::ostream;

EvtSpinDensity::EvtSpinDensity( const EvtSpinDensity& density )
{
    dim = 0;
    rho = 0;

    int i, j;
    setDim( density.dim );

    for ( i = 0; i < dim; i++ ) {
        for ( j = 0; j < dim; j++ ) {
            rho[i][j] = density.rho[i][j];
        }
    }
}

EvtSpinDensity& EvtSpinDensity::operator=( const EvtSpinDensity& density )
{
    int i, j;
    setDim( density.dim );

    for ( i = 0; i < dim; i++ ) {
        for ( j = 0; j < dim; j++ ) {
            rho[i][j] = density.rho[i][j];
        }
    }

    return *this;
}

EvtSpinDensity::~EvtSpinDensity()
{
    if ( dim != 0 ) {
        int i;
        for ( i = 0; i < dim; i++ )
            delete[] rho[i];
    }

    delete[] rho;
}

EvtSpinDensity::EvtSpinDensity()
{
    dim = 0;
    rho = 0;
}

void EvtSpinDensity::setDim( int n )
{
    if ( dim == n )
        return;
    if ( dim != 0 ) {
        int i;
        for ( i = 0; i < dim; i++ )
            delete[] rho[i];
        delete[] rho;
        rho = 0;
        dim = 0;
    }
    if ( n == 0 )
        return;
    dim = n;
    rho = new EvtComplexPtr[n];
    int i;
    for ( i = 0; i < n; i++ ) {
        rho[i] = new EvtComplex[n];
    }
}

int EvtSpinDensity::getDim() const
{
    return dim;
}

void EvtSpinDensity::set( int i, int j, const EvtComplex& rhoij )
{
    assert( i < dim && j < dim );
    rho[i][j] = rhoij;
}

const EvtComplex& EvtSpinDensity::get( int i, int j ) const
{
    assert( i < dim && j < dim );
    return rho[i][j];
}

void EvtSpinDensity::setDiag( int n )
{
    setDim( n );
    int i, j;

    for ( i = 0; i < n; i++ ) {
        for ( j = 0; j < n; j++ ) {
            rho[i][j] = EvtComplex( 0.0 );
        }
        rho[i][i] = EvtComplex( 1.0 );
    }
}

double EvtSpinDensity::normalizedProb( const EvtSpinDensity& d )
{
    int i, j;
    EvtComplex prob( 0.0, 0.0 );
    double norm = 0.0;

    if ( dim != d.dim ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Not matching dimensions in NormalizedProb" << endl;
        ::abort();
    }

    for ( i = 0; i < dim; i++ ) {
        norm += real( rho[i][i] );
        for ( j = 0; j < dim; j++ ) {
            prob += rho[i][j] * d.rho[i][j];
        }
    }

    if ( imag( prob ) > 0.00000001 * real( prob ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Imaginary probability:" << prob << " " << norm << endl;
    }
    if ( real( prob ) < 0.0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Negative probability:" << prob << " " << norm << endl;
    }

    return real( prob ) / norm;
}

int EvtSpinDensity::check()
{
    if ( dim < 1 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "dim=" << dim << "in SpinDensity::Check" << endl;
    }

    int i, j;

    double trace( 0.0 );

    for ( i = 0; i < dim; i++ ) {
        trace += abs( rho[i][i] );
    }

    for ( i = 0; i < dim; i++ ) {
        if ( real( rho[i][i] ) < 0.0 )
            return 0;
        if ( imag( rho[i][i] ) * 1000000.0 > trace ) {
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << *this << endl;
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << trace << endl;
            EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Failing 1" << endl;
            return 0;
        }
    }

    for ( i = 0; i < dim; i++ ) {
        for ( j = i + 1; j < dim; j++ ) {
            if ( fabs( real( rho[i][j] - rho[j][i] ) ) >
                 0.00000001 * ( abs( rho[i][i] ) + abs( rho[j][j] ) ) ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Failing 2" << endl;
                return 0;
            }
            if ( fabs( imag( rho[i][j] + rho[j][i] ) ) >
                 0.00000001 * ( abs( rho[i][i] ) + abs( rho[j][j] ) ) ) {
                EvtGenReport( EVTGEN_INFO, "EvtGen" ) << "Failing 3" << endl;
                return 0;
            }
        }
    }

    return 1;
}

ostream& operator<<( ostream& s, const EvtSpinDensity& d )
{
    int i, j;

    s << endl;
    s << "Dimension:" << d.dim << endl;

    for ( i = 0; i < d.dim; i++ ) {
        for ( j = 0; j < d.dim; j++ ) {
            s << d.rho[i][j] << " ";
        }
        s << endl;
    }

    return s;
}
