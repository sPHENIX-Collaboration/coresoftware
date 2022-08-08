
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

#ifndef EVTRARITASCHWINGERPARTICLE_HH
#define EVTRARITASCHWINGERPARTICLE_HH

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRaritaSchwinger.hh"

class EvtVector4R;

class EvtRaritaSchwingerParticle : public EvtParticle {
  public:
    EvtRaritaSchwingerParticle() = default;
    void init( EvtId id, const EvtVector4R& p4 ) override;
    void init( EvtId id, const EvtVector4R& p4, const EvtRaritaSchwinger&,
               const EvtRaritaSchwinger&, const EvtRaritaSchwinger&,
               const EvtRaritaSchwinger&, const EvtRaritaSchwinger&,
               const EvtRaritaSchwinger&, const EvtRaritaSchwinger&,
               const EvtRaritaSchwinger& );
    EvtRaritaSchwinger spRSParent( int ) const override;
    EvtRaritaSchwinger spRS( int ) const override;
    EvtSpinDensity rotateToHelicityBasis() const override;
    EvtSpinDensity rotateToHelicityBasis( double alpha, double beta,
                                          double gamma ) const override;

  private:
    EvtRaritaSchwinger _spinorRest[4];
    EvtRaritaSchwinger _spinor[4];
    EvtRaritaSchwingerParticle( const EvtRaritaSchwingerParticle& d );
    EvtRaritaSchwingerParticle& operator=( const EvtRaritaSchwingerParticle& d );
};
#endif
