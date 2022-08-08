
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

#ifndef EVTRARELBTOLLLFFGUTSCHE_HH
#define EVTRARELBTOLLLFFGUTSCHE_HH 1

// Include files

/** @class EvtRareLbToLllFF EvtRareLbToLllFF.hh EvtGenModels/EvtRareLbToLllFF.hh
 *
 *
 *  @author Michal Kreps
 *  @date   2014-10-21
 */

#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtParticle.hh"

#include "EvtGenModels/EvtRareLbToLllFFBase.hh"

#include <map>
#include <string>

class EvtRareLbToLllFFGutsche : public EvtRareLbToLllFFBase {
  public:
    void init() override;

    void getFF( EvtParticle* parent, EvtParticle* lambda,
                EvtRareLbToLllFFBase::FormFactors& FF ) override;

  private:
    double formFactorParametrization( double s, double f0, double a, double b );

    double fVconsts[3][3];
    double fAconsts[3][3];
    double fTVconsts[3][3];
    double fTAconsts[3][3];

    static EvtIdSet fParents;
    static EvtIdSet fDaughters;
};

#endif    // EVTRARELBTOLLLFF_HH
