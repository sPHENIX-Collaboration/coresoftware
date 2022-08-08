
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

#include "EvtGenModels/EvtBToDiBaryonlnupQCDFF.hh"

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"

EvtBToDiBaryonlnupQCDFF::EvtBToDiBaryonlnupQCDFF() : DPars(), nDPars( 0 )
{
}

EvtBToDiBaryonlnupQCDFF::EvtBToDiBaryonlnupQCDFF( std::vector<double>& DParameters ) :
    DPars( DParameters ), nDPars( DParameters.size() )
{
}

void EvtBToDiBaryonlnupQCDFF::getFF( EvtParticle*, double dibaryonMass,
                                     EvtBToDiBaryonlnupQCDFF::FormFactors& FF ) const
{
    if ( nDPars == 6 && dibaryonMass > 0.0 ) {
        // 5/3*[1/M^2]^3
        double t = 5.0 / ( 3.0 * pow( dibaryonMass, 6.0 ) );

        double Dp = DPars[0];
        double Dpb = DPars[1];
        double D2 = DPars[2];
        double D3 = DPars[3];
        double D4 = DPars[4];
        double D5 = DPars[5];

        FF.F1 = ( Dp + 0.2 * Dpb ) * t;
        FF.F2 = -D2 * t;
        FF.F3 = -D3 * t;
        FF.F4 = -D4 * t;
        FF.F5 = -D5 * t;

        FF.G1 = ( Dp - 0.2 * Dpb ) * t;
        FF.G2 = -FF.F2;
        FF.G3 = -FF.F3;
        FF.G4 = -FF.F4;
        FF.G5 = -FF.F5;
    }
}

void EvtBToDiBaryonlnupQCDFF::getDiracFF(
    EvtParticle* parent, double dibaryonMass,
    EvtBToDiBaryonlnupQCDFF::FormFactors& FF ) const
{
    this->getFF( parent, dibaryonMass, FF );
}

void EvtBToDiBaryonlnupQCDFF::getRaritaFF(
    EvtParticle* parent, double dibaryonMass,
    EvtBToDiBaryonlnupQCDFF::FormFactors& FF ) const
{
    this->getFF( parent, dibaryonMass, FF );
}
