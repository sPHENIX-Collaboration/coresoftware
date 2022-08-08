
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

#include "EvtGenModels/EvtbsToLLLL.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"

#include "EvtGenModels/EvtbTosllWilsCoeffNLO.hh"
#include "EvtGenModels/Evtbs2llGammaFFMNT.hh"
#include "EvtGenModels/EvtbsToLLLLAmp.hh"

#include <stdlib.h>
#include <string.h>

EvtbsToLLLL::~EvtbsToLLLL()
{
    delete _mntffmodel;
    if ( _calcamp )
        delete _calcamp;
}

// The module name specification
std::string EvtbsToLLLL::getName()
{
    return "BQTOLLLL";
}

// The implementation of the clone() method
EvtDecayBase* EvtbsToLLLL::clone()
{
    return new EvtbsToLLLL;
}

// The inicialization of the decay model
//
// Tn the our model we have are following 4 arguments:
//
//           mu          - the scale parameter, GeV;
//           Nf          - number of "effective" flavors (for b-quark Nf=5);
//           res_swch    - resonant switching parametr:
//                         = 0 the resonant contribution switched OFF,
//                         = 1 the resonant contribution switched ON;
//           ias         - switching parametr for \alpha_s(M_Z) value:
//                         = 0 PDG 1sigma minimal alpha_s(M_Z),
//                         = 1 PDG average value  alpha_s(M_Z),
//                         = 2 PDG 1sigma maximal alpha_s(M_Z).
//           Wolfenstein parameterization for CKM matrix
//                         CKM_A, CKM_lambda, CKM_barrho, CKM_bareta
//
void EvtbsToLLLL::init()
{
    // check that there are 8 arguments
    checkNArg( 8 );
    // check that there are 4 daughteres
    checkNDaug( 4 );

    // We expect that the parent to be a scalar (B-meson)
    // and the daughters to be l^+, l^-, l^+ and l^-
    checkSpinParent( EvtSpinType::SCALAR );

    // We expect that the all daughters are the ell+ or ell- == DIRAC
    checkSpinDaughter( 0, EvtSpinType::DIRAC );
    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::DIRAC );
    checkSpinDaughter( 3, EvtSpinType::DIRAC );

    _mntffmodel = new Evtbs2llGammaFFMNT();
    _wilscoeff = new EvtbTosllWilsCoeffNLO();
    _calcamp = new EvtbsToLLLLAmp();
}

// Set the maximum probability of the decay
void EvtbsToLLLL::initProbMax()
{
    double mymaxprob = -10.0;    // maximum of the probability

    //  EvtId parnum, l1num, l2num, l3num, l4num;

    //  parnum = getParentId();
    //  l1num  = getDaug(0);
    //  l2num  = getDaug(1);
    //  l3num  = getDaug(2);
    //  l4num  = getDaug(3);

    // EvtSpinType::spintype mesontype=EvtPDL::getSpinType(getDaug(0));

    //  double         mu = getArg(0);        // the scale parameter
    //  int            Nf = (int) getArg(1);  // number of "effective" flavors
    //  int      res_swch = (int) getArg(2);  // resonant switching parametr
    //  int           ias = (int) getArg(3);  // switching parametr for \alpha_s(M_Z)
    //  double      CKM_A = getArg(4);
    //  double CKM_lambda = getArg(5);
    //  double CKM_barrho = getArg(6);
    //  double CKM_bareta = getArg(7);

    mymaxprob = _calcamp->CalcMaxProb(
        //                                   parnum, l1num,l2num, l3num,l4num,
        //		                    _mntffmodel, _wilscoeff,
        //                                   mu, Nf, res_swch, ias,
        //                                   CKM_A,CKM_lambda,CKM_barrho,CKM_bareta
    );

    if ( mymaxprob <= 0.0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "The function void EvtbsToLLLL::initProbMax()"
            << "\n Unexpected value of the probability maximum!"
            << "\n mymaxprob = " << mymaxprob << std::endl;
        ::abort();
    }

    setProbMax( mymaxprob );
}

void EvtbsToLLLL::decay( EvtParticle* p )
{
    double mu = getArg( 0 );            // the scale parameter
    int Nf = (int)getArg( 1 );          // number of "effective" flavors
    int res_swch = (int)getArg( 2 );    // resonant switching parametr
    int ias = (int)getArg( 3 );         // switching parametr for \alpha_s(M_Z)
    double CKM_A = getArg( 4 );
    double CKM_lambda = getArg( 5 );
    double CKM_barrho = getArg( 6 );
    double CKM_bareta = getArg( 7 );

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    _calcamp->CalcAmp( p, _amp2, _mntffmodel, _wilscoeff, mu, Nf, res_swch, ias,
                       CKM_A, CKM_lambda, CKM_barrho, CKM_bareta );

    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen") << "\n The function EvtbTosllMSExt::decay(...) passed with arguments:"
    //                        << "\n mu = " << mu << " Nf =" << Nf
    //                        << " res_swch = " << res_swch
    //                        << " ias = " << ias
    //                        << " CKM_A = " << CKM_A
    //                        << " CKM_lambda = " << CKM_lambda
    //                        << " CKM_barrho = " << CKM_barrho
    //                        << " CKM_bareta = " << CKM_bareta
    //                        << " ReA7 = " << ReA7
    //                        << " ImA7 = " << ImA7
    //                        << " ReA10 = " << ReA10
    //                        << " ImA10 = " << ImA10 << std::endl;
}
