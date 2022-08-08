
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

#ifndef EVTPHOTONPARTICLE_HH
#define EVTPHOTONPARTICLE_HH

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtVector4C.hh"
class EvtId;

//Class to handle massless spin 1 particles.

class EvtPhotonParticle : public EvtParticle {
  public:
    EvtPhotonParticle() = default;

    void init( EvtId part_n, double e, double px, double py, double pz );
    void init( EvtId part_n, const EvtVector4R& p4 ) override;

    //Return polarization vectors
    EvtVector4C epsParentPhoton( int i ) override;
    EvtVector4C epsPhoton( int i ) override;

    EvtSpinDensity rotateToHelicityBasis() const override;
    EvtSpinDensity rotateToHelicityBasis( double alpha, double beta,
                                          double gamma ) const override;

  private:
    EvtVector4C eps1, eps2;
    int _evalBasis;

    EvtPhotonParticle( const EvtPhotonParticle& photon );
    EvtPhotonParticle& operator=( const EvtPhotonParticle& photon );
};

#endif
