
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

#ifndef EVT_DALITZ_POINT_HH
#define EVT_DALITZ_POINT_HH

#include "EvtGenBase/EvtCyclic3.hh"
#include "EvtGenBase/EvtDalitzCoord.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"

// This class describes the complete kinematics of the Dalitz decay.
// It holds all the six invariant momentum products, three daughter
// particle masses and three invariant masses of pairs of particles.
// This description is completely symmetric with respect to particle
// permutations.
//
// Another way to slice the six coordinate is to make a transformation
// to the mass of the decaying particle. The four masses make up a
// Dalitz plot. The other two are coordinates of a point in the plot.

class EvtDalitzPoint final {
  public:
    EvtDalitzPoint();
    EvtDalitzPoint( double mA, double mB, double mC, double qAB, double qBC,
                    double qCA );
    EvtDalitzPoint( double mA, double mB, double mC, EvtCyclic3::Pair i,
                    double qres, double qhel, double qsum );
    EvtDalitzPoint( const EvtDalitzPlot&, const EvtDalitzCoord& );

    EvtDalitzCoord getDalitzPoint( EvtCyclic3::Pair i, EvtCyclic3::Pair j ) const;
    EvtDalitzPlot getDalitzPlot() const;

    double q( EvtCyclic3::Pair ) const;
    double bigM() const;
    double m( EvtCyclic3::Index ) const;

    // Zemach variables

    double qres( EvtCyclic3::Pair i ) const;
    double qhel( EvtCyclic3::Pair i ) const;
    double qsum() const;

    // Kinematic quantities
    //
    // pp  - 4 momentum product
    // e,p,cosTh - energy/moementum in rest-frame of j

    double qMin( EvtCyclic3::Pair i, EvtCyclic3::Pair j ) const;
    double qMax( EvtCyclic3::Pair i, EvtCyclic3::Pair j ) const;
    double pp( EvtCyclic3::Index i, EvtCyclic3::Index j ) const;
    double e( EvtCyclic3::Index i, EvtCyclic3::Pair j ) const;
    double p( EvtCyclic3::Index i, EvtCyclic3::Pair j ) const;
    double cosTh( EvtCyclic3::Pair pairAng, EvtCyclic3::Pair pairRes ) const;

    bool isValid() const;

    void print() const;

  private:
    double _mA, _mB, _mC;       // masses
    double _qAB, _qBC, _qCA;    // masses squared
};

#endif
