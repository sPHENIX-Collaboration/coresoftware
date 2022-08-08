
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

#ifndef EVT_MASSAMP_HH
#define EVT_MASSAMP_HH

#include "EvtGenBase/EvtAmplitude.hh"
#include "EvtGenBase/EvtPoint1D.hh"
#include "EvtGenBase/EvtPropBreitWignerRel.hh"
#include "EvtGenBase/EvtTwoBodyVertex.hh"

// Relativistic lineshape for a two-body decay of a resonance to two
// pseudoscalars. The mass dependence of the width and the vertex factors
// are included in the calculation.

class EvtMassAmp : public EvtAmplitude<EvtPoint1D> {
  public:
    EvtMassAmp( const EvtPropBreitWignerRel& prop, const EvtTwoBodyVertex& vd );
    EvtMassAmp( const EvtMassAmp& other );
    EvtMassAmp& operator=( const EvtMassAmp& other );

    EvtComplex amplitude( const EvtPoint1D& p ) const override;

    EvtAmplitude<EvtPoint1D>* clone() const override
    {
        return new EvtMassAmp( *this );
    }

    void setBirthVtx( const EvtTwoBodyVertex& vb )
    {
        _vb = std::make_unique<EvtTwoBodyVertex>( vb );
    }

    void addBirthFact() { _useBirthFact = true; }
    void addDeathFact() { _useDeathFact = true; }
    void addBirthFactFF() { _useBirthFactFF = true; }
    void addDeathFactFF() { _useDeathFactFF = true; }

  private:
    EvtPropBreitWignerRel _prop;
    EvtTwoBodyVertex _vd;
    std::unique_ptr<EvtTwoBodyVertex> _vb;

    bool _useBirthFact;
    bool _useDeathFact;
    bool _useBirthFactFF;
    bool _useDeathFactFF;
};

#endif
