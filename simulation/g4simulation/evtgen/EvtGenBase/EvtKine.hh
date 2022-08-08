
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

#ifndef EVTKINE_HH
#define EVTKINE_HH

class EvtVector4R;
class EvtComplex;

double EvtDecayAngle( const EvtVector4R&, const EvtVector4R&, const EvtVector4R& );

double EvtDecayAngleChi( const EvtVector4R&, const EvtVector4R&,
                         const EvtVector4R&, const EvtVector4R&,
                         const EvtVector4R& );

//
// This routine calculates the cosine of the angle between
// the normal of the decay plane and the flight direction of particle q
// in the parent frame.
//
double EvtDecayPlaneNormalAngle( const EvtVector4R& p, const EvtVector4R& q,
                                 const EvtVector4R& d1, const EvtVector4R& d2 );

// Added by AJB
//
// Calculate phi (between 0 and 2 pi) of the daughter given the 4-momentum of
// the grandparent, parent, resonance and the daughter.  g, p, q and d need to
// be in the same rest frame.  Note that for the first level of the tree there
// is no grandparent and thus <0,0,0,1> should be passed in for g.  When there
// is no parent the angles need to be calculated by simply by calculating polar
// and azymuthal angles in the rest frame of the resonance (since this will
// generally be the root particle and is generally at rest the polar and
// azymuthal angels can simply be calculated.
//
double EvtDecayAnglePhi( const EvtVector4R& g, const EvtVector4R& p,
                         const EvtVector4R& q, const EvtVector4R& d );

// Wigner big-D function in Jackson convention
//
// XXX NOTE XXX
//  - EvtDecayAngle returns the cos \theta and EvtdFunction requires theta
//  - In EvtdFunction j m1 and m2 are really 2 * j, 2 * m1, 2*m2 to deal with
//    spin 1/2 particles
//
EvtComplex wignerD( int j, int m1, int m2, double phi, double theta,
                    double gamma );

//
// Function to calculate momentum of daughters in two body decay in mothers
// rest frame.
double twoBodyMomentum( const double M, const double m1, const double m2 );

#endif
