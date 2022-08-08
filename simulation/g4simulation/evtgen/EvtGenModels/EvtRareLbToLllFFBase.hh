
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

#ifndef EVTRARELBTOLLLFFBASE_HH
#define EVTRARELBTOLLLFFBASE_HH 1

// Include files

/** @class 
 *  
 *
 *  @author Michal Kreps
 *  @date   2014-10-20
 */

#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtParticle.hh"

#include <map>
#include <string>

class EvtRareLbToLllFFBase {
  public:
    class FormFactors {
      public:
        FormFactors();

        virtual ~FormFactors(){};

        void areZero();

        double F_[4];
        double G_[4];
        double FT_[4];
        double GT_[4];
    };

    virtual void init() = 0;

    virtual void getFF( EvtParticle* parent, EvtParticle* lambda,
                        EvtRareLbToLllFFBase::FormFactors& FF ) = 0;

    bool isNatural( EvtParticle* lambda );

    EvtRareLbToLllFFBase();
    virtual ~EvtRareLbToLllFFBase(){};

  protected:
    double calculateVdotV( EvtParticle* parent, EvtParticle* lambda ) const;
    double calculateVdotV( EvtParticle*, EvtParticle*, double qsq ) const;

    EvtIdSet natural_;
};

#endif
