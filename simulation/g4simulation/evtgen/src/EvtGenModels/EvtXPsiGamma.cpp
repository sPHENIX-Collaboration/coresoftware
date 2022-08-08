
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

#include "EvtGenModels/EvtXPsiGamma.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtTensorParticle.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <iostream>
#include <stdlib.h>
#include <string>

using namespace std;

std::string EvtXPsiGamma::getName()
{
    return "X38722-+_PSI_GAMMA";
}

EvtDecayBase* EvtXPsiGamma::clone()
{
    //  cout<<" (* AVL: === EvtXPsiGamma::clone() ============ *)"<<endl;
    return new EvtXPsiGamma;
}

EvtComplex EvtXPsiGamma::fT2( EvtVector4R p, EvtVector4R q, EvtTensor4C epsPI,
                              EvtVector4C epsEps, EvtVector4C epsEta )
{
    // T2 term from [Bazi](10)
    EvtTensor4C epsPQ =
        EvtGenFunctions::directProd( q, p );    // e_{mu nu a b} p^a q^b;
    epsPQ = dual( epsPQ );

    EvtVector4C tmp1 = epsPI.cont1( epsEps );
    EvtVector4C tmp2 = epsPQ.cont1( tmp1 );
    EvtComplex T2 =
        tmp2 * epsEta;    // epa^a pi_{a mu} e_{mu nu rho si} p_nu q_rho eta_si

    tmp1 = epsPI.cont1( epsEta );
    tmp2 = epsPQ.cont1( tmp1 );
    T2 += tmp2 *
          epsEps;    // T2 - eta^a pi_{a mu} e_{mu nu rho si} q_nu p_rhi eps_si

    return T2;
}

EvtComplex EvtXPsiGamma::fT3( EvtVector4R p, EvtVector4R q, EvtTensor4C epsPI,
                              EvtVector4C epsEps, EvtVector4C epsEta )
{
    // T3 term from [Bazi](11)
    EvtVector4R Q = p - q, P = p + q;
    EvtVector4C tmp1 = epsPI.cont1( Q );    // Q_a pi_{a mu}
    EvtTensor4C tmp3 = dual( EvtGenFunctions::directProd(
        P, epsEps ) );    // e_{mu nu rho si} P^rho eps^si
    EvtVector4C tmp4 = tmp3.cont1( tmp1 );
    EvtComplex T3 =
        tmp4 * epsEta;    // Q_a pi_{a mu} e_{mu nu rho si} P^rho eps_si eta_nu
    return T3;
}

void EvtXPsiGamma::decay( EvtParticle* root )
{
    ncall++;
    root->initializePhaseSpace( getNDaug(), getDaugs() );

    double gOmega = 1.58,
           gPOmega = -0.74;    // X -> omega psi couplings from table II
    double gRho = 1.58, gPRho = -0.74;    // X -> omega psi couplings from table II
    double fRho = 0.121, mRho2 = 0.770 * 0.770, fOmega = 0.036,
           mOmega2 = 0.782 * 0.782;

    EvtComplex amp;

    if ( _ID0 == EvtPDL::getId( "gamma" ) ) {
        for ( int iPsi = 0; iPsi < 4; iPsi++ ) {
            for ( int iGamma = 0; iGamma < 1; iGamma++ ) {
                for ( int iChi = 0; iChi < 4; iChi++ ) {
                    EvtComplex T2 = fT2(
                        root->getDaug( 1 )->getP4(),
                        root->getDaug( 0 )->getP4(), root->epsTensor( iChi ),
                        root->getDaug( 1 )->epsParent( iPsi ).conj(),
                        root->getDaug( 0 )->epsParentPhoton( iGamma ).conj() );
                    EvtComplex T3 = fT3(
                        root->getDaug( 1 )->getP4(),
                        root->getDaug( 0 )->getP4(), root->epsTensor( iChi ),
                        root->getDaug( 1 )->epsParent( iPsi ).conj(),
                        root->getDaug( 0 )->epsParentPhoton( iGamma ).conj() );
                    amp = ( fOmega / mOmega2 * gOmega + fRho / mRho2 * gRho ) *
                              T2 +
                          ( fOmega / mOmega2 * gPOmega + fRho / mRho2 * gPRho ) *
                              T3;
                    vertex( iChi, iGamma, iPsi, amp );
                };
            };
        };
    } else if ( _ID0 == EvtPDL::getId( "omega" ) ) {
        for ( int iPsi = 0; iPsi < 4; iPsi++ ) {
            for ( int iGamma = 0; iGamma < 4; iGamma++ ) {
                for ( int iChi = 0; iChi < 4; iChi++ ) {
                    EvtComplex T2 = fT2(
                        root->getDaug( 1 )->getP4(),
                        root->getDaug( 0 )->getP4(), root->epsTensor( iChi ),
                        root->getDaug( 1 )->epsParent( iPsi ).conj(),
                        root->getDaug( 0 )->epsParent( iGamma ).conj() );
                    EvtComplex T3 = fT3(
                        root->getDaug( 1 )->getP4(),
                        root->getDaug( 0 )->getP4(), root->epsTensor( iChi ),
                        root->getDaug( 1 )->epsParent( iPsi ).conj(),
                        root->getDaug( 0 )->epsParent( iGamma ).conj() );
                    //	  cout << "AVL:: omega"<<endl;
                    amp = gOmega * T2 + gPOmega * T3;
                    vertex( iChi, iGamma, iPsi, amp );
                };
            };
        };
    } else if ( _ID0 == EvtPDL::getId( "rho0" ) ) {
        for ( int iPsi = 0; iPsi < 4; iPsi++ ) {
            for ( int iGamma = 0; iGamma < 4; iGamma++ ) {
                for ( int iChi = 0; iChi < 4; iChi++ ) {
                    EvtComplex T2 = fT2(
                        root->getDaug( 1 )->getP4(),
                        root->getDaug( 0 )->getP4(), root->epsTensor( iChi ),
                        root->getDaug( 1 )->epsParent( iPsi ).conj(),
                        root->getDaug( 0 )->epsParent( iGamma ).conj() );
                    EvtComplex T3 = fT3(
                        root->getDaug( 1 )->getP4(),
                        root->getDaug( 0 )->getP4(), root->epsTensor( iChi ),
                        root->getDaug( 1 )->epsParent( iPsi ).conj(),
                        root->getDaug( 0 )->epsParent( iGamma ).conj() );
                    //	  cout << "AVL:: rho"<<endl;
                    amp = gRho * T2 + gPRho * T3;
                    vertex( iChi, iGamma, iPsi, amp );
                };
            };
        };
    } else {
        cout << "AVL:: Not realized yet" << endl;
    };
}

void EvtXPsiGamma::init()
{
    //  cout<<" (* AVL: ==== EvtXPsiGamma::init() ============ *)"<<endl;

    ncall = 0;

    checkNArg( 0 );
    checkNDaug( 2 );

    checkSpinParent( EvtSpinType::TENSOR );

    //  checkSpinDaughter(0,EvtSpinType::PHOTON);
    checkSpinDaughter( 1, EvtSpinType::VECTOR );

    _ID0 = getDaug( 0 );
    /*  if(_ID0 == EvtPDL::getId("gamma") ) {
    cout << "AVL:: gamma"<<endl;
  }
  else if(_ID0 == EvtPDL::getId("omega") ) {
    cout << "AVL:: omega"<<endl;
  }
  else if(_ID0 == EvtPDL::getId("rho0") ) {
    cout << "AVL:: rho"<<endl;
  };
*/
}

void EvtXPsiGamma::initProbMax()
{
    if ( _ID0 == EvtPDL::getId( "gamma" ) )
        setProbMax( 2.400 );
    else if ( _ID0 == EvtPDL::getId( "omega" ) )
        setProbMax( 16. );
    else if ( _ID0 == EvtPDL::getId( "rho0" ) )
        setProbMax( 70. );
}
