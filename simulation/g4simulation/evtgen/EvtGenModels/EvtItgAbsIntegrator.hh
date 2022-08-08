
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

#ifndef EVTITGABSINTEGRATOR_HH
#define EVTITGABSINTEGRATOR_HH

#include "EvtGenModels/EvtItgAbsFunction.hh"

// Description:
//      Abstraction of a generic integrator (Stolen and modified from
//      the BaBar IntegrationUtils package - author: Phil Strother).

class EvtItgAbsIntegrator {
  public:
    EvtItgAbsIntegrator( const EvtItgAbsFunction& );

    virtual ~EvtItgAbsIntegrator() = default;

    double evaluate( double lower, double upper ) const;

    double normalisation() const;

  protected:
    double trapezoid( double lower, double higher, int n, double& result ) const;

    virtual double evaluateIt( double lower, double higher ) const = 0;

    double myFunction( double x ) const { return _myFunction( x ); }

  private:
    const EvtItgAbsFunction& _myFunction;

    void boundsCheck( double&, double& ) const;

    // Note: if your class needs a copy constructor or an assignment operator,
    //  make one of the following public and implement it.
    EvtItgAbsIntegrator();
    //EvtItgAbsIntegrator( const EvtItgAbsIntegrator& );                // Copy Constructor
    //EvtItgAbsIntegrator& operator= ( const EvtItgAbsIntegrator& );    // Assignment op
};

#endif    // EVTITGABSINTEGRATOR_HH
