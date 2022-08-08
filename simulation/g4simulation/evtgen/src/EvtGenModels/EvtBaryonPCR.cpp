
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

#include "EvtGenModels/EvtBaryonPCR.hh"

#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include "EvtGenModels/EvtBaryonPCRFF.hh"

#include <stdlib.h>
#include <string>

using namespace std;
#ifdef D0
#undef D0
#endif

std::string EvtBaryonPCR::getName()
{
    return "BaryonPCR";
}

EvtBaryonPCR* EvtBaryonPCR::clone()
{
    return new EvtBaryonPCR;
}

void EvtBaryonPCR::decay( EvtParticle* p )
{
    //This is a kludge to avoid warnings because the K_2* mass becomes to large.
    static EvtIdSet regenerateMasses( "K_2*+", "K_2*-", "K_2*0", "anti-K_2*0",
                                      "K_1+", "K_1-", "K_10", "anti-K_10",
                                      "D'_1+", "D'_1-", "D'_10", "anti-D'_10" );

    if ( regenerateMasses.contains( getDaug( 0 ) ) ) {
        p->resetFirstOrNot();
    }

    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtComplex r00( getArg( 0 ), 0.0 );
    EvtComplex r01( getArg( 1 ), 0.0 );
    EvtComplex r10( getArg( 2 ), 0.0 );
    EvtComplex r11( getArg( 3 ), 0.0 );

    calcamp->CalcAmp( p, _amp2, baryonpcrffmodel.get(), r00, r01, r10, r11 );
}

void EvtBaryonPCR::initProbMax()
{
    // Baryons (partial list 5/28/04)

    static EvtId SIGC0 = EvtPDL::getId( "Sigma_c0" );
    static EvtId SIGC0B = EvtPDL::getId( "anti-Sigma_c0" );
    static EvtId SIGCP = EvtPDL::getId( "Sigma_c+" );
    static EvtId SIGCM = EvtPDL::getId( "anti-Sigma_c-" );
    static EvtId SIGCPP = EvtPDL::getId( "Sigma_c++" );
    static EvtId SIGCMM = EvtPDL::getId( "anti-Sigma_c--" );
    static EvtId LAMCP = EvtPDL::getId( "Lambda_c+" );
    static EvtId LAMCM = EvtPDL::getId( "anti-Lambda_c-" );
    static EvtId LAMC1P = EvtPDL::getId( "Lambda_c(2593)+" );
    static EvtId LAMC1M = EvtPDL::getId( "anti-Lambda_c(2593)-" );
    static EvtId LAMC2P = EvtPDL::getId( "Lambda_c(2625)+" );
    static EvtId LAMC2M = EvtPDL::getId( "anti-Lambda_c(2625)-" );
    static EvtId LAMB = EvtPDL::getId( "Lambda_b0" );
    static EvtId LAMBB = EvtPDL::getId( "anti-Lambda_b0" );

    EvtId parnum, barnum, lnum;

    parnum = getParentId();
    barnum = getDaug( 0 );
    lnum = getDaug( 1 );

    if ( parnum == LAMB || parnum == LAMBB ) {
        if ( barnum == LAMCP || barnum == LAMCM || barnum == LAMC1P ||
             barnum == LAMC1M || barnum == LAMC2P || barnum == LAMC2M ||
             barnum == SIGC0 || barnum == SIGC0B || barnum == SIGCP ||
             barnum == SIGCM || barnum == SIGCPP || barnum == SIGCMM ) {
            setProbMax( 22000.0 );
            return;
        }
    }

    //This is a real cludge.. (ryd)
    setProbMax( 0.0 );
}

void EvtBaryonPCR::init()
{
    //if (getNArg()!=0) {
    if ( getNArg() != 4 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBaryonPCR generator expected "
            << " 4 arguments but found:" << getNArg() << endl;
        //<< " 0 arguments but found:"<<getNArg()<<endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }

    if ( getNDaug() != 3 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Wrong number of daughters in EvtBaryonPCR.cc "
            << " 3 daughters expected but found: " << getNDaug() << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }

    //We expect the parent to be a scalar
    //and the daughters to be X lepton neutrino

    EvtSpinType::spintype parenttype = EvtPDL::getSpinType( getParentId() );
    EvtSpinType::spintype baryontype = EvtPDL::getSpinType( getDaug( 0 ) );
    EvtSpinType::spintype leptontype = EvtPDL::getSpinType( getDaug( 1 ) );
    EvtSpinType::spintype neutrinotype = EvtPDL::getSpinType( getDaug( 2 ) );

    if ( parenttype != EvtSpinType::DIRAC ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBaryonPCR generator expected "
            << " a DIRAC parent, found:" << EvtPDL::name( getParentId() )
            << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
    if ( leptontype != EvtSpinType::DIRAC ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBaryonPCR generator expected "
            << " a DIRAC 2nd daughter, found:" << EvtPDL::name( getDaug( 1 ) )
            << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
    if ( neutrinotype != EvtSpinType::NEUTRINO ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBaryonPCR generator expected "
            << " a NEUTRINO 3rd daughter, found:" << EvtPDL::name( getDaug( 2 ) )
            << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }

    baryonpcrffmodel = std::make_unique<EvtBaryonPCRFF>();

    if ( baryontype == EvtSpinType::DIRAC ||
         baryontype == EvtSpinType::RARITASCHWINGER ) {
        calcamp = std::make_unique<EvtSemiLeptonicBaryonAmp>();
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Wrong baryon spin type in EvtBaryonPCR.cc "
            << "Expected spin type " << EvtSpinType::DIRAC
            << ", found spin type " << baryontype << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << endl;
        ::abort();
    }
}
