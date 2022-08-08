
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

#ifndef EVTTENSORPARTICLE_HH
#define EVTTENSORPARTICLE_HH

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtTensor4C.hh"

class EvtTensorParticle : public EvtParticle {
  public:
    EvtTensorParticle() = default;

    void init( EvtId part_n, double e, double px, double py, double pz );
    void init( EvtId part_n, const EvtVector4R& p4 ) override;
    void init( EvtId part_n, const EvtVector4R& p4, const EvtTensor4C&,
               const EvtTensor4C&, const EvtTensor4C&, const EvtTensor4C&,
               const EvtTensor4C& );
    //Returns polarization tensors.
    EvtTensor4C epsTensorParent( int i ) const override;
    EvtTensor4C epsTensor( int i ) const override;

    EvtSpinDensity rotateToHelicityBasis() const override;
    EvtSpinDensity rotateToHelicityBasis( double alpha, double beta,
                                          double gamma ) const override;

  private:
    EvtTensor4C eps[5];    //eps1,eps2,eps3,eps4,eps5;

    EvtTensorParticle( const EvtTensorParticle& tensor );
    EvtTensorParticle& operator=( const EvtTensorParticle& tensor );
};

#endif
