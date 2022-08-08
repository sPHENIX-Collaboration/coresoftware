
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

#ifndef EVTRARELBTOLLLFFLQCD_HH
#define EVTRARELBTOLLLFFLQCD_HH 1

// Include files

/** @class EvtRareLbToLllFF EvtRareLbToLllFF.hh EvtGenModels/EvtRareLbToLllFF.hh
 *
 *
 *  @author Michal Kreps
 *  @date   2016-04-19
 *  @date   2014-10-23
 */

#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtParticle.hh"

#include "EvtGenModels/EvtRareLbToLllFFBase.hh"

#include <map>
#include <string>

class EvtRareLbToLllFFlQCD : public EvtRareLbToLllFFBase {
  public:
    /// Standard constructor
    EvtRareLbToLllFFlQCD() = default;

    void init() override;

    void getFF( EvtParticle* parent, EvtParticle* lambda,
                EvtRareLbToLllFFBase::FormFactors& FF ) override;

  private:
    double formFactorParametrization( double q2, double a0, double a1,
                                      double pole );
    double zvar( double q2 );

    double fconsts[3][3];
    double gconsts[3][3];
    double hconsts[3][3];
    double htildaconsts[3][3];

    double t0;
    double tplus;
};

#endif    // EVTRARELBTOLLLFF_HH
