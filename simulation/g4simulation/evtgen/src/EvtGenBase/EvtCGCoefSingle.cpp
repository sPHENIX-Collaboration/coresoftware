
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

#include "EvtGenBase/EvtCGCoefSingle.hh"

#include "EvtGenBase/EvtOrthogVector.hh"
#include "EvtGenBase/EvtPatches.hh"

#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>

void EvtCGCoefSingle::init( int j1, int j2 )
{
    _j1 = j1;
    _j2 = j2;

    _Jmax = abs( j1 + j2 );
    _Jmin = abs( j1 - j2 );

    _table.resize( ( _Jmax - _Jmin ) / 2 + 1 );

    int J, M;

    int lenmax = j1 + 1;
    if ( j2 < j1 )
        lenmax = j2 + 1;

    //set vector sizes
    for ( J = _Jmax; J >= _Jmin; J -= 2 ) {
        _table[( J - _Jmin ) / 2].resize( J + 1 );
        for ( M = J; J >= -M; M -= 2 ) {
            int len = ( ( _j1 + _j2 ) - abs( M ) ) / 2 + 1;
            if ( len > lenmax )
                len = lenmax;
            _table[( J - _Jmin ) / 2][( M + J ) / 2].resize( len );
        }
    }

    //now fill the vectors
    for ( J = _Jmax; J >= _Jmin; J -= 2 ) {
        //bootstrap with highest M(=J) as a special case
        if ( J == _Jmax ) {
            cg( J, J, _j1, _j2 ) = 1.0;
        } else {
            int n = ( _Jmax - J ) / 2 + 1;
            std::vector<double>* vectors = new std::vector<double>[n - 1];
            int i, k;
            for ( i = 0; i < n - 1; i++ ) {
                // i corresponds to J=Jmax-2*i
                vectors[i].resize( n );
                for ( k = 0; k < n; k++ ) {
                    double tmp = _table[( _Jmax - _Jmin ) / 2 - i]
                                       [( J + _Jmax - 2 * i ) / 2][k];
                    vectors[i][k] = tmp;
                }
            }
            EvtOrthogVector getOrth( n, vectors );
            std::vector<double> orth = getOrth.getOrthogVector();
            int sign = 1;
            if ( orth[n - 1] < 0.0 )
                sign = -1;
            for ( k = 0; k < n; k++ ) {
                _table[( J - _Jmin ) / 2][J][k] = sign * orth[k];
            }
            delete[] vectors;
        }
        for ( M = J - 2; M >= -J; M -= 2 ) {
            int len = ( ( _j1 + _j2 ) - abs( M ) ) / 2 + 1;
            if ( len > lenmax )
                len = lenmax;
            int mmin = M - j2;
            if ( mmin < -j1 )
                mmin = -j1;
            int m1;
            for ( m1 = mmin; m1 < mmin + len * 2; m1 += 2 ) {
                int m2 = M - m1;
                double sum = 0.0;
                float fkwTmp = _j1 * ( _j1 + 2 ) - ( m1 + 2 ) * m1;
                //fkw 2/2/2001: changes needed to satisfy KCC
                //fkw if (m1+2<=_j1) sum+=0.5*sqrt(_j1*(_j1+2)-(m1+2)*m1)*cg(J,M+2,m1+2,m2);
                //fkw if (m2+2<=_j2) sum+=0.5*sqrt(_j2*(_j2+2)-(m2+2)*m2)*cg(J,M+2,m1,m2+2);
                //fkw sum/=(0.5*sqrt(J*(J+2)-(M+2)*M));
                if ( m1 + 2 <= _j1 )
                    sum += 0.5 * sqrt( fkwTmp ) * cg( J, M + 2, m1 + 2, m2 );
                fkwTmp = _j2 * ( _j2 + 2 ) - ( m2 + 2 ) * m2;
                if ( m2 + 2 <= _j2 )
                    sum += 0.5 * sqrt( fkwTmp ) * cg( J, M + 2, m1, m2 + 2 );
                fkwTmp = J * ( J + 2 ) - ( M + 2 ) * M;
                sum /= ( 0.5 * sqrt( fkwTmp ) );
                cg( J, M, m1, m2 ) = sum;
            }
        }
    }
}

double EvtCGCoefSingle::coef( int J, int M, int j1, int j2, int m1, int m2 )
{
    assert( j1 == _j1 );
    _unused( j1 );
    assert( j2 == _j2 );
    _unused( j2 );

    return cg( J, M, m1, m2 );
}

double& EvtCGCoefSingle::cg( int J, int M, int m1, int m2 )
{
    assert( M == m1 + m2 );
    _unused( m2 );
    assert( abs( M ) <= J );
    assert( J <= _Jmax );
    assert( J >= _Jmin );
    assert( abs( m1 ) <= _j1 );
    assert( abs( m2 ) <= _j2 );

    //find lowest m1 allowed for the given M

    int mmin = M - _j2;

    if ( mmin < -_j1 )
        mmin = -_j1;

    int n = m1 - mmin;

    return _table[( J - _Jmin ) / 2][( M + J ) / 2][n / 2];
}
