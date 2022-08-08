
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

#ifndef EvtWHad_HH
#define EvtWHad_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <vector>

// Description: Routine to calculate W -> (n pi) + (m K) current
//			according to [Kuhn, Was, Acta.Phys.Polon B39 (2008) 147]

class EvtWHad {
  public:
    EvtWHad();

    EvtVector4C WCurrent( const EvtVector4R& q1 ) const;

    EvtVector4C WCurrent( const EvtVector4R& q1, const EvtVector4R& q2 ) const;

    EvtVector4C WCurrent( const EvtVector4R& q1, const EvtVector4R& q2,
                          const EvtVector4R& q3 ) const;

    EvtVector4C WCurrent( const EvtVector4R& q1, const EvtVector4R& q2,
                          const EvtVector4R& q3, const EvtVector4R& q4,
                          const EvtVector4R& q5 ) const;

    EvtVector4C WCurrent_KKP( const EvtVector4R& pKplus,
                              const EvtVector4R& pKminus,
                              const EvtVector4R& pPiPlus ) const;

    EvtVector4C WCurrent_KPP( const EvtVector4R& pKplus,
                              const EvtVector4R& pPiPlus,
                              const EvtVector4R& pPiMinus ) const;

    EvtVector4C WCurrent_KSK( const EvtVector4R& pKS,
                              const EvtVector4R& pKplus ) const;

  protected:
    EvtVector4C JB( const EvtVector4R& q1, const EvtVector4R& q2,
                    const EvtVector4R& q3, const EvtVector4R& q4,
                    const EvtVector4R& q5 ) const;

    EvtComplex BWa( const EvtVector4R& q ) const;

    EvtComplex BWf( const EvtVector4R& q ) const;

    EvtComplex BWr( const EvtVector4R& q ) const;

    EvtComplex BWKK( double s, int i ) const;

    double pi3G( double Q2 ) const;

    EvtComplex pcm( double s ) const;

  private:
    std::vector<double> mRho_, gamma0_, cK_;
};

#endif
