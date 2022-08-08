
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

#include "EvtGenModels/Evtbs2llGammaISRFSR.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"

#include "EvtGenModels/EvtbTosllWilsCoeffNLO.hh"
#include "EvtGenModels/Evtbs2llGammaFFMNT.hh"
#include "EvtGenModels/Evtbs2llGammaISRFSRAmp.hh"

#include <stdlib.h>
#include <string.h>

Evtbs2llGammaISRFSR::~Evtbs2llGammaISRFSR()
{
    delete _mntffmodel;
    if ( _calcamp )
        delete _calcamp;
}

// The module name specification
std::string Evtbs2llGammaISRFSR::getName()
{
    return "BSTOGLLISRFSR";
}

// The implementation of the clone() method
EvtDecayBase* Evtbs2llGammaISRFSR::clone()
{
    return new Evtbs2llGammaISRFSR;
}

// The inicialization of the decay model
//
// Tn the our model we have are following 4 arguments:
//
//           mu          - the scale parameter, GeV;
//           Nf          - number of "effective" flavors (for b-quark Nf=5);
//           sr          - state radiation type:
//                         = 0 ISR only,
//                         = 1 FSR only,
//                         = 2 ISR + FSR (== BSTOGLLMNT model);
//           res_swch    - resonant switching parametr:
//                         = 0 the resonant contribution switched OFF,
//                         = 1 the resonant contribution switched ON;
//           ias         - switching parametr for \alpha_s(M_Z) value:
//                         = 0 PDG 1sigma minimal alpha_s(M_Z),
//                         = 1 PDG average value  alpha_s(M_Z),
//                         = 2 PDG 1sigma maximal alpha_s(M_Z).
//	     Egamma_min  - photon energy cut, GeV;
//           Wolfenstein parameterization for CKM matrix
//                         CKM_A, CKM_lambda, CKM_barrho, CKM_bareta
//
void Evtbs2llGammaISRFSR::init()
{
    // check that there are 10 arguments
    checkNArg( 10, 11 );
    // check that there are 3 daughteres
    checkNDaug( 3 );

    // We expect that the parent to be a scalar (B-meson)
    // and the daughters to be Gamma, l^+ and l^-
    checkSpinParent( EvtSpinType::SCALAR );

    // We expect that the first daughter is the photon
    EvtSpinType::spintype photontype = EvtPDL::getSpinType( getDaug( 0 ) );

    if ( !( photontype == EvtSpinType::PHOTON ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Evtbs2llGammaISRFSR generator expected "
            << " a PHOTON 1st daughter, found:"
            << EvtPDL::name( getDaug( 0 ) ).c_str() << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << std::endl;
        ::abort();
    }

    // We expect that the second and third daughters
    // are the ell+ and ell- == DIRAC
    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::DIRAC );

    _mntffmodel = new Evtbs2llGammaFFMNT();
    _wilscoeff = new EvtbTosllWilsCoeffNLO();
    if ( photontype == EvtSpinType::PHOTON ) {
        _calcamp = new Evtbs2llGammaISRFSRAmp();
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "The init()-function in the Evtbs2llGammaISRFSR generator:"
            << "The absence in the radiative decay!" << std::endl;
        ::abort();
    }
}

// Set the maximum probability of the decay
// differencial distribution d^2\Gamma/d\hat s d\cos\theta
void Evtbs2llGammaISRFSR::initProbMax()
{
    double mymaxprob = -10.0;    // maximum of the probability

    EvtId parnum, photnum, l1num, l2num;

    parnum = getParentId();
    photnum = getDaug( 0 );
    l1num = getDaug( 1 );
    l2num = getDaug( 2 );

    double mu = getArg( 0 );            // the scale parameter
    int Nf = (int)getArg( 1 );          // number of "effective" flavors
    int sr = (int)getArg( 2 );          // state radiation type
    int res_swch = (int)getArg( 3 );    // resonant switching parametr
    int ias = (int)getArg( 4 );         // switching parametr for \alpha_s(M_Z)
    double Egamma_min = getArg( 5 );    // photon energy cut
    double CKM_A = getArg( 6 );
    double CKM_lambda = getArg( 7 );
    double CKM_barrho = getArg( 8 );
    double CKM_bareta = getArg( 9 );
    double mumumass_min = 0.;
    if ( getNArg() == 11 )
        mumumass_min = getArg( 10 );

    mymaxprob = _calcamp->CalcMaxProb( parnum, photnum, l1num, l2num, _mntffmodel,
                                       _wilscoeff, mu, Nf, sr, res_swch, ias,
                                       Egamma_min, CKM_A, CKM_lambda,
                                       CKM_barrho, CKM_bareta, mumumass_min );

    if ( mymaxprob <= 0.0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "The function void Evtbs2llGammaISRFSR::initProbMax()"
            << "\n Unexpected value of the probability maximum!"
            << "\n mymaxprob = " << mymaxprob << std::endl;
        ::abort();
    }

    setProbMax( mymaxprob );
}

void Evtbs2llGammaISRFSR::decay( EvtParticle* p )
{
    double mu = getArg( 0 );            // the scale parameter
    int Nf = (int)getArg( 1 );          // number of "effective" flavors
    int sr = (int)getArg( 2 );          // state radiation type
    int res_swch = (int)getArg( 3 );    // resonant switching parametr
    int ias = (int)getArg( 4 );         // switching parametr for \alpha_s(M_Z)
    double Egamma_min = getArg( 5 );    // photon energy cut
    double CKM_A = getArg( 6 );
    double CKM_lambda = getArg( 7 );
    double CKM_barrho = getArg( 8 );
    double CKM_bareta = getArg( 9 );
    double mumumass_min = 0.;
    if ( getNArg() == 11 )
        mumumass_min = getArg( 10 );

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    // The class "Evtbs2llGammaFFMNT" is the derived class of the
    // class  "Evtbs2llGammaFF" (see the file "Evtbs2llGammaFF.hh")
    _calcamp->CalcAmp( p, _amp2, _mntffmodel, _wilscoeff, mu, Nf, sr, res_swch,
                       ias, Egamma_min, CKM_A, CKM_lambda, CKM_barrho,
                       CKM_bareta, mumumass_min );

    //  EvtGenReport(EVTGEN_NOTICE,"EvtGen") << "\n "
    //<< "\n The function Evtbs2llGammaISRFSR::decay(...) passed with arguments:"
    //                          << "\n mu = " << mu << " Nf =" << Nf
    //                          << "\n res_swch = " << res_swch
    //                          << "\n sr = " << sr << "\n ias = " << ias
    //		 	              << "\n Egamma_min =" << Egamma_min
    //                          << "\n CKM_A = " << CKM_A
    //                          << "\n CKM_lambda = " << CKM_lambda
    //                          << "\n CKM_barrho = " << CKM_barrho
    //                          << "\n CKM_bareta = " << CKM_bareta
    //                          << "\n "
    //                          << std::endl;
}
