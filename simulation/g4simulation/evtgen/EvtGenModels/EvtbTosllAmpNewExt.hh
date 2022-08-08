
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

#ifndef EVTBTOSLLAMPNEWEXT_HH
#define EVTBTOSLLAMPNEWEXT_HH

#include "EvtGenBase/EvtId.hh"

class EvtAmp;
class EvtbTosllFFNew;
class EvtParticle;
class EvtbTosllWilsCoeffNLO;

// all definicions !
class EvtbTosllAmpNewExt {
  public:
    virtual ~EvtbTosllAmpNewExt(){};

    virtual void CalcAmp( EvtParticle* parent, EvtAmp& amp,
                          EvtbTosllFFNew* formFactors,
                          EvtbTosllWilsCoeffNLO* WilsCoeff, double /* mu */,
                          int /* Nf */, int /* res_swch */, int /* ias */,
                          double /* CKM_A */, double /* CKM_lambda */,
                          double /* CKM_barrho */, double /* CKM_bareta */,
                          double /* ReA7 */, double /* ImA7 */,
                          double /* ReA10 */, double /*ImA10 */ ) = 0;

    virtual double CalcMaxProb( EvtId /* parnum */, EvtId /* mesnum */,
                                EvtId /* l1num */, EvtId /* l2num */,
                                EvtbTosllFFNew* /* formFactors */,
                                EvtbTosllWilsCoeffNLO* /*WilsCoeff */,
                                double /* mu */, int /* Nf */,
                                int /* res_swch */, int /* ias */,
                                double /* CKM_A */, double /* CKM_lambda */,
                                double /* CKM_barrho */, double /* CKM_bareta */,
                                double /* ReA7 */, double /* ImA7 */,
                                double /* ReA10 */, double /*ImA10 */ )
    {
        return 0.0;
    };

    virtual double lambda( double /* a */, double /* b */, double /* c */ )
    {
        return 0.0;
    };
};

#endif
