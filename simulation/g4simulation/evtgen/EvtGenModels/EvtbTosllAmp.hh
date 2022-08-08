
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

#ifndef EVTBTOSLLAMP_HH
#define EVTBTOSLLAMP_HH

class EvtAmp;
class EvtId;
class EvtbTosllFF;
class EvtParticle;
class EvtComplex;

class EvtbTosllAmp {
  public:
    virtual ~EvtbTosllAmp(){};

    //Daughters are initialized and have been added to the parent.
    //No need to carry around the daughters seperately!

    virtual void CalcAmp( EvtParticle* parent, EvtAmp& amp,
                          EvtbTosllFF* formFactors ) = 0;

    double CalcMaxProb( EvtId parent, EvtId meson, EvtId lepton, EvtId nudaug,
                        EvtbTosllFF* formFactors, double& poleSize );

    EvtComplex GetC7Eff( double q2, bool nnlo = true );
    EvtComplex GetC9Eff( double q2, bool nnlo = true, bool btod = false );
    EvtComplex GetC10Eff( double q2, bool nnlo = true );

    double dGdsProb( double mb, double ms, double ml, double s );

    double dGdsdupProb( double mb, double ms, double ml, double s, double u );
};

#endif
