
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

#ifndef EVTITGSIMPSONINTEGRATOR_HH
#define EVTITGSIMPSONINTEGRATOR_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

#include "EvtGenModels/EvtItgAbsIntegrator.hh"

// Description:
//      Simpson integrator (Stolen and modified from
//      the BaBar IntegrationUtils package - author: Phil Strother).

class EvtItgSimpsonIntegrator : public EvtItgAbsIntegrator {
  public:
    EvtItgSimpsonIntegrator( const EvtItgAbsFunction&,
                             double precision = 1.0e-5, int maxLoop = 20 );

  protected:
    double evaluateIt( double, double ) const override;

  private:
    double _precision;
    double _maxLoop;

    //EvtItgSimpsonIntegrator( const EvtItgSimpsonIntegrator& );                //// Copy Constructor
    //EvtItgSimpsonIntegrator& operator= ( const EvtItgSimpsonIntegrator& );    // Assignment op
};

#endif    // ITGSIMPSONINTEGRATOR_HH
