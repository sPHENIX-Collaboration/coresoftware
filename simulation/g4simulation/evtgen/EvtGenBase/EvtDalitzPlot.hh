
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

#ifndef EVT_DALITZ_PLOT_HH
#define EVT_DALITZ_PLOT_HH

#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDecayMode.hh"
#include "EvtGenBase/EvtTwoBodyVertex.hh"

#include <assert.h>

class EvtDalitzPlot {
  public:
    EvtDalitzPlot();
    EvtDalitzPlot( double mA, double mB, double mC, double bigM,
                   double ldel = 0., double rdel = 0. );
    EvtDalitzPlot( const EvtDecayMode& mode, double ldel = 0., double rdel = 0. );
    bool operator==( const EvtDalitzPlot& other ) const;
    const EvtDalitzPlot* clone() const;

    // Absolute limits for masses squared in the Dalitz plot
    // e.g. qAbsMin(0) is the lowest possible value
    // for m2 of particles {12}

    double qAbsMin( EvtCyclic3::Pair i ) const;
    double qAbsMax( EvtCyclic3::Pair i ) const;
    double mAbsMin( EvtCyclic3::Pair i ) const;
    double mAbsMax( EvtCyclic3::Pair i ) const;

    // Absolute limits for Zemach coordinate qres and qhel (approximate)
    // qHelAbsMin(BC,CA) means absolute minimum for (qCA-qAB)/2.

    double qResAbsMin( EvtCyclic3::Pair i ) const;
    double qResAbsMax( EvtCyclic3::Pair i ) const;
    double qHelAbsMin( EvtCyclic3::Pair i ) const;
    double qHelAbsMax( EvtCyclic3::Pair i ) const;
    inline double qSumMin() const { return sum() + _ldel; }
    inline double qSumMax() const { return sum() + _rdel; }
    inline bool fuzzy() const { return ( _rdel - _ldel != 0. ); }

    // Find the area of the Dalitz plot by numeric integration. (N bins for variable q(i) are used).
    // Very large numbers of N can result in a very long calculation. It should not
    // matter which two pairs f variables are used. The integral should eventually
    // converge to the same number

    double getArea( int N = 1000, EvtCyclic3::Pair i = EvtCyclic3::AB,
                    EvtCyclic3::Pair j = EvtCyclic3::BC ) const;

    // Limits for masses squared when one mass squared is known

    double qMin( EvtCyclic3::Pair i, EvtCyclic3::Pair j, double q ) const;
    double qMax( EvtCyclic3::Pair i, EvtCyclic3::Pair j, double q ) const;

    // Coordinate transformations

    double cosTh( EvtCyclic3::Pair i1, double q1, EvtCyclic3::Pair i2,
                  double q2 ) const;
    double e( EvtCyclic3::Index i, EvtCyclic3::Pair j, double q ) const;
    double p( EvtCyclic3::Index i, EvtCyclic3::Pair j, double q ) const;

    double q( EvtCyclic3::Pair i1, double cosTh, EvtCyclic3::Pair i2,
              double q2 ) const;

    // |J| of transformation of qi to cosTh in the rest-frame of j

    double jacobian( EvtCyclic3::Pair i, double q ) const;

    // Given resonance index and mass returns decay
    // and birth vertices

    EvtTwoBodyVertex vD( EvtCyclic3::Pair iRes, double m0, int L ) const;
    EvtTwoBodyVertex vB( EvtCyclic3::Pair iRes, double m0, int L ) const;

    // Accessors

    double sum() const;
    inline double bigM() const { return _bigM; }
    inline double mA() const { return _mA; }
    inline double mB() const { return _mB; }
    inline double mC() const { return _mC; }
    double m( EvtCyclic3::Index i ) const;

    void print() const;

    void sanityCheck() const;

  protected:
    // Defines two dimensional dalitz plot

    double _mA;
    double _mB;
    double _mC;
    double _bigM;

    // Defines third dimension, or fuzziness. M^2 + ldel < M^2 < M^2 + rdel

    double _ldel;
    double _rdel;
};

#endif
