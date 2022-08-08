
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

#ifndef EVTDIRACPARTICLE_HH
#define EVTDIRACPARTICLE_HH

#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtParticle.hh"

class EvtId;
class EvtVector4R;

class EvtDiracParticle : public EvtParticle {
  public:
    EvtDiracParticle() = default;
    void init( EvtId part_n, const EvtVector4R& p4 ) override;
    void init( EvtId part_n, const EvtVector4R& p4, const EvtDiracSpinor&,
               const EvtDiracSpinor&, const EvtDiracSpinor&,
               const EvtDiracSpinor& );
    EvtDiracSpinor spParent( int i ) const override { return _spinorParent[i]; }
    EvtDiracSpinor sp( int i ) const override { return _spinorRest[i]; }
    EvtSpinDensity rotateToHelicityBasis() const override;
    EvtSpinDensity rotateToHelicityBasis( double alpha, double beta,
                                          double gamma ) const override;

  private:
    EvtDiracSpinor _spinorRest[2];
    EvtDiracSpinor _spinorParent[2];
    EvtDiracParticle( const EvtDiracParticle& d );
    EvtDiracParticle& operator=( const EvtDiracParticle& d );
};
#endif
