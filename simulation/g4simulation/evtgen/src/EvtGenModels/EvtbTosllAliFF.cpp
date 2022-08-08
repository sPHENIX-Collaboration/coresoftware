
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

#include "EvtGenModels/EvtbTosllAliFF.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtPatches.hh"

#include <math.h>

EvtbTosllAliFF::EvtbTosllAliFF()
{
}

void EvtbTosllAliFF::getScalarFF( EvtId parent, EvtId /*daught*/, double t,
                                  double /*mass*/, double& fp, double& f0,
                                  double& ft )
{
    double m = EvtPDL::getMeanMass( parent );
    //double md=EvtPDL::getMeanMass(daught);

    double shat = t / ( m * m );
    double shat2 = shat * shat;
    double shat3 = shat2 * shat;

    fp = 0.278 * exp( 1.568 * shat + 0.470 * shat2 + 0.885 * shat3 );
    f0 = 0.278 * exp( 0.740 * shat + 0.080 * shat2 + 0.425 * shat3 );
    ft = 0.300 * exp( 1.600 * shat + 0.501 * shat2 + 0.796 * shat3 );
}

void EvtbTosllAliFF::getVectorFF( EvtId parent, EvtId /*daught*/, double t,
                                  double /*mass*/, double& a1, double& a2,
                                  double& a0, double& v, double& t1, double& t2,
                                  double& t3 )
{
    double m = EvtPDL::getMeanMass( parent );

    double shat = t / ( m * m );
    double shat2 = shat * shat;

    //this is Ali 'minimum allowed form factors'
    a1 = 0.294 * exp( 0.656 * shat + 0.456 * shat2 );
    a2 = 0.246 * exp( 1.237 * shat + 0.822 * shat2 );
    a0 = 0.412 * exp( 1.543 * shat + 0.954 * shat2 );
    v = 0.399 * exp( 1.537 * shat + 1.123 * shat2 );

    t1 = 0.334 * exp( 1.575 * shat + 1.140 * shat2 );
    t2 = 0.334 * exp( 0.562 * shat + 0.481 * shat2 );
    t3 = 0.234 * exp( 1.230 * shat + 1.089 * shat2 );
}
