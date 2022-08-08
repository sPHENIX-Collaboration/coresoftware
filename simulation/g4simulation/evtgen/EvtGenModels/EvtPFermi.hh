
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

#ifndef EVTPFERMI_HH
#define EVTPFERMI_HH

// Description:
//   DFN model:
//      F(k+) = N (1-x)^a exp((1+a)x) ,x=k+/(mB-mb)
//      the fermi motion distribution according to
//      hep-ph/9905351 v2
//   BLNP model:
//      F(what,Lambda,b) = pow(_b,_b)/(tgamma(_b)*_Lambda)*pow(what/_Lambda,_b-1)*
//                           exp(-_b*what/Lambda);
//      the leading order shape function (exp) (hep-ph/0504071)

class EvtPFermi final {
  public:
    // Constructors

    EvtPFermi( const double& a, const double& mB, const double& mb );
    EvtPFermi( const double& Lambda, const double& b );

    // Operators

    // Selectors

    // Modifiers

    // Methods

    double getFPFermi( const double& kplus );
    double getSFBLNP( const double& what );

  protected:
    // Helper functions

  private:
    // Friends

    // Data members

    double _a;
    double _mb;
    double _mB;
    double _Lambda;
    double _b;
};

#endif    // EVTPFERMI_HH
