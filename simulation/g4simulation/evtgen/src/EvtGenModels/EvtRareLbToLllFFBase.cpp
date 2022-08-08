
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

#include "EvtGenModels/EvtRareLbToLllFFBase.hh"

EvtRareLbToLllFFBase::FormFactors::FormFactors()
{
    areZero();
}

void EvtRareLbToLllFFBase::FormFactors::areZero()
{
    for ( unsigned int i = 0; i < 4; ++i ) {
        F_[i] = 0;
        G_[i] = 0;
        FT_[i] = 0;
        GT_[i] = 0;
    }
}

EvtRareLbToLllFFBase::EvtRareLbToLllFFBase() :
    natural_( "Lambda0", "anti-Lambda0", "Lambda(1520)0", "anti-Lambda(1520)0",
              "Lambda(1600)0", "anti-Lambda(1600)0" )
{
}

bool EvtRareLbToLllFFBase::isNatural( EvtParticle* lambda )
{
    return natural_.contains( lambda->getId() );
}

double EvtRareLbToLllFFBase::calculateVdotV( EvtParticle* parent,
                                             EvtParticle* lambda ) const
{
    EvtVector4R p4parent;
    p4parent.set( parent->mass(), 0, 0, 0 );

    EvtVector4R p4lambda = lambda->getP4();

    double M = lambda->mass();
    double MB = parent->mass();

    return p4parent.cont( p4lambda ) / ( MB * M );

    // return E_Lambda/M_Lambda
}

double EvtRareLbToLllFFBase::calculateVdotV( EvtParticle* parent,
                                             EvtParticle* lambda,
                                             const double qsq ) const
{
    double M = lambda->mass();
    double MB = parent->mass();

    double E = ( MB * MB - M * M - qsq ) / ( 2. * MB );

    return E / M;
}
