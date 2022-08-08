
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

#ifndef EVTVECTORPARTICLE_HH
#define EVTVECTORPARTICLE_HH

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"

class EvtId;

class EvtVectorParticle : public EvtParticle {
  public:
    EvtVectorParticle() = default;

    void init( EvtId part_n, double e, double px, double py, double pz );
    void init( EvtId part_n, const EvtVector4R& p ) override;
    void init( EvtId part_n, const EvtVector4R& p, const EvtVector4C&,
               const EvtVector4C&, const EvtVector4C& );
    EvtVector4C epsParent( int i ) const override
    {
        return boostTo( _eps[i], this->getP4() );
    }
    EvtVector4C eps( int i ) const override { return _eps[i]; }
    EvtSpinDensity rotateToHelicityBasis() const override;
    EvtSpinDensity rotateToHelicityBasis( double alpha, double beta,
                                          double gamma ) const override;

  private:
    EvtVector4C _eps[3];

    EvtVectorParticle( const EvtVectorParticle& vector );
    EvtVectorParticle& operator=( const EvtVectorParticle& vector );
};

#endif
