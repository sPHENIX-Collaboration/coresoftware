
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

#ifndef EvtWNPI_HH
#define EvtWNPI_HH

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

// Description: Routine to calculate W -> (n pi) current
//			according to [Kuhn, Was, Acta.Phys.Polon B39 (2008) 147]

class EvtWnPi {
  public:
    EvtVector4C WCurrent( EvtVector4R q1 );
    EvtVector4C WCurrent( EvtVector4R q1, EvtVector4R q2 );
    EvtVector4C WCurrent( EvtVector4R q1, EvtVector4R q2, EvtVector4R q3 );
    EvtVector4C WCurrent( EvtVector4R q1, EvtVector4R q2, EvtVector4R q3,
                          EvtVector4R q4, EvtVector4R q5 );

  protected:
    EvtVector4C JB( EvtVector4R q1, EvtVector4R q2, EvtVector4R q3,
                    EvtVector4R q4, EvtVector4R q5 );
    EvtComplex BWa( EvtVector4R q );
    EvtComplex BWf( EvtVector4R q );
    EvtComplex BWr( EvtVector4R q );
    double pi3G( double Q2 );
};

#endif
