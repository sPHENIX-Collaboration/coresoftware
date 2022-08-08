
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

#ifndef EVTGENKINE_HH
#define EVTGENKINE_HH

class EvtVector4R;
class EvtParticle;

class EvtGenKine {
  public:
    static double PhaseSpace( int ndaug, double mass[30], EvtVector4R p4[30],
                              double mp );

    static double PhaseSpacePole( double M, double m1, double m2, double m3,
                                  double a, EvtVector4R p4[10] );

    /*
     * Function which takes two invariant masses squared in 3-body decay and
     * parent after makeDaughters() and generateMassTree() and
     * calculates/generates momenta of daughters and sets those.
     */
    static void ThreeBodyKine( const double m12Sq, const double m23Sq,
                               EvtParticle* p );
};

#endif
