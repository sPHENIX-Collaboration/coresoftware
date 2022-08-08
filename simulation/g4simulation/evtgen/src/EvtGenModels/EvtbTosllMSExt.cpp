
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

#include "EvtGenModels/EvtbTosllMSExt.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"

#include "EvtGenModels/EvtbTosllAmpNewExt.hh"
#include "EvtGenModels/EvtbTosllMSFF.hh"
#include "EvtGenModels/EvtbTosllScalarAmpNewExt.hh"
#include "EvtGenModels/EvtbTosllVectorAmpNewExt.hh"
#include "EvtGenModels/EvtbTosllWilsCoeffNLO.hh"

#include <stdlib.h>
#include <string.h>

EvtbTosllMSExt::~EvtbTosllMSExt()
{
    delete _msffmodel;
    if ( _calcamp )
        delete _calcamp;
}

// The module name specification
std::string EvtbTosllMSExt::getName()
{
    return "BTOSLLMSEXT";
}

// The implementation of the clone() method
EvtDecayBase* EvtbTosllMSExt::clone()
{
    return new EvtbTosllMSExt;
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
//           Complex multiplication coefficients
//                         A7  = ReA7  + i*ImA7
//                         A10 = ReA10 + i*ImA10
//
void EvtbTosllMSExt::init()
{
    // check that there are 12 arguments
    checkNArg( 12 );
    // check that there are 3 daughteres
    checkNDaug( 3 );

    // We expect that the parent to be a scalar (B-meson)
    // and the daughters to be K^*, l^+ and l^-
    checkSpinParent( EvtSpinType::SCALAR );

    // We expect that the first daughter is the K* == VECTOR
    EvtSpinType::spintype mesontype = EvtPDL::getSpinType( getDaug( 0 ) );

    if ( !( mesontype == EvtSpinType::VECTOR ||
            mesontype == EvtSpinType::SCALAR ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtbTosllMSExt generator expected "
            << " a SCALAR or VECTOR 1st daughter, found:"
            << EvtPDL::name( getDaug( 0 ) ).c_str() << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << std::endl;
        ::abort();
    }

    // We expect that the second and third daughters
    // are the ell+ and ell- == DIRAC
    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::DIRAC );

    _msffmodel = new EvtbTosllMSFF();
    _wilscoeff = new EvtbTosllWilsCoeffNLO();
    if ( mesontype == EvtSpinType::VECTOR ) {
        _calcamp = new EvtbTosllVectorAmpNewExt();
    }
    if ( mesontype == EvtSpinType::SCALAR ) {
        _calcamp = new EvtbTosllScalarAmpNewExt();
    }
}

// Set the maximum probability of the decay
// differencial distribution d^2\Gamma/d\hat s d\cos\theta
void EvtbTosllMSExt::initProbMax()
{
    double mymaxprob = -10.0;    // maximum of the probability

    EvtId parnum, mesnum, l1num, l2num;

    parnum = getParentId();
    mesnum = getDaug( 0 );
    l1num = getDaug( 1 );
    l2num = getDaug( 2 );

    // EvtSpinType::spintype mesontype=EvtPDL::getSpinType(getDaug(0));

    double mu = getArg( 0 );            // the scale parameter
    int Nf = (int)getArg( 1 );          // number of "effective" flavors
    int res_swch = (int)getArg( 2 );    // resonant switching parametr
    int ias = (int)getArg( 3 );         // switching parametr for \alpha_s(M_Z)
    double CKM_A = getArg( 4 );
    double CKM_lambda = getArg( 5 );
    double CKM_barrho = getArg( 6 );
    double CKM_bareta = getArg( 7 );
    double ReA7 = getArg( 8 );
    double ImA7 = getArg( 9 );
    double ReA10 = getArg( 10 );
    double ImA10 = getArg( 11 );

    mymaxprob = _calcamp->CalcMaxProb( parnum, mesnum, l1num, l2num, _msffmodel,
                                       _wilscoeff, mu, Nf, res_swch, ias, CKM_A,
                                       CKM_lambda, CKM_barrho, CKM_bareta, ReA7,
                                       ImA7, ReA10, ImA10 );

    if ( mymaxprob <= 0.0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "The function void EvtbTosllMSExt::initProbMax()"
            << "\n Unexpected value of the probability maximum!"
            << "\n mymaxprob = " << mymaxprob << std::endl;
        ::abort();
    }

    setProbMax( mymaxprob );
}

void EvtbTosllMSExt::decay( EvtParticle* p )
{
    double mu = getArg( 0 );            // the scale parameter
    int Nf = (int)getArg( 1 );          // number of "effective" flavors
    int res_swch = (int)getArg( 2 );    // resonant switching parametr
    int ias = (int)getArg( 3 );         // switching parametr for \alpha_s(M_Z)
    double CKM_A = getArg( 4 );
    double CKM_lambda = getArg( 5 );
    double CKM_barrho = getArg( 6 );
    double CKM_bareta = getArg( 7 );
    double ReA7 = getArg( 8 );
    double ImA7 = getArg( 9 );
    double ReA10 = getArg( 10 );
    double ImA10 = getArg( 11 );

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    // The class "EvtbTosllVectorAmpNewExt" is the derived class of the
    // class  "EvtbTosllAmpNewExt" (see the file "EvtbTosllVectorAmpNewExt.hh")
    // and
    // the class "EvtbTosllMSFF" is the derived class of the
    // class  "EvtbTosllFFNew" (see the file "EvtbTosllMSFF.hh")
    _calcamp->CalcAmp( p, _amp2, _msffmodel, _wilscoeff, mu, Nf, res_swch, ias,
                       CKM_A, CKM_lambda, CKM_barrho, CKM_bareta, ReA7, ImA7,
                       ReA10, ImA10 );

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
