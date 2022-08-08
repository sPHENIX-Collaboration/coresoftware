
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

#include "EvtGenModels/EvtBToDiBaryonlnupQCD.hh"

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtVector4R.hh"

std::string EvtBToDiBaryonlnupQCD::getName()
{
    return "BToDiBaryonlnupQCD";
}

EvtDecayBase* EvtBToDiBaryonlnupQCD::clone()
{
    return new EvtBToDiBaryonlnupQCD;
}

void EvtBToDiBaryonlnupQCD::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs(), true );

    calcAmp_->CalcAmp( p, _amp2 );
}

void EvtBToDiBaryonlnupQCD::init()
{
    if ( !( getNArg() == 6 || getNArg() == 7 ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBToDiBaryonlnupQCD model expected "
            << " 6 or 7 arguments but found:" << getNArg() << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << std::endl;
        ::abort();
    }

    if ( getNDaug() != 4 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Wrong number of daughters in EvtBToDiBaryonlnupQCD model: "
            << "4 daughters expected but found: " << getNDaug() << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << std::endl;
        ::abort();
    }

    // We expect B -> baryon baryon lepton neutrino
    EvtSpinType::spintype parentType = EvtPDL::getSpinType( getParentId() );
    EvtSpinType::spintype leptonType = EvtPDL::getSpinType( getDaug( 2 ) );
    EvtSpinType::spintype neutrinoType = EvtPDL::getSpinType( getDaug( 3 ) );

    if ( parentType != EvtSpinType::SCALAR ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBToDiBaryonlnupQCD model expected "
            << " a SCALAR parent, found:" << EvtPDL::name( getParentId() )
            << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << std::endl;
        ::abort();
    }

    if ( leptonType != EvtSpinType::DIRAC ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBToDiBaryonlnupQCD model expected "
            << " a DIRAC 3rd daughter, found:" << EvtPDL::name( getDaug( 2 ) )
            << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << std::endl;
        ::abort();
    }

    if ( neutrinoType != EvtSpinType::NEUTRINO ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBToDiBaryonlnupQCD model expected "
            << " a NEUTRINO 4th daughter, found:" << EvtPDL::name( getDaug( 3 ) )
            << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << std::endl;
        ::abort();
    }

    // Get the 6 form factor D parameters from model arguments in the decay file
    std::vector<double> DPars( 6 );
    for ( int i = 0; i < 6; i++ ) {
        DPars[i] = getArg( i );
    }

    // Form factor model
    ffModel_ = std::make_unique<EvtBToDiBaryonlnupQCDFF>( DPars );

    // Set amplitude calculation pointer.
    // Accomodate for spin 1/2 (DIRAC) or 3/2 (RARITASCHWINGER) baryons
    EvtSpinType::spintype baryon1Type = EvtPDL::getSpinType( getDaug( 0 ) );
    EvtSpinType::spintype baryon2Type = EvtPDL::getSpinType( getDaug( 1 ) );

    if ( ( baryon1Type == EvtSpinType::DIRAC &&
           baryon2Type == EvtSpinType::RARITASCHWINGER ) ||
         ( baryon1Type == EvtSpinType::RARITASCHWINGER &&
           baryon2Type == EvtSpinType::DIRAC ) ||
         ( baryon1Type == EvtSpinType::DIRAC &&
           baryon2Type == EvtSpinType::DIRAC ) ) {
        calcAmp_ = std::make_unique<EvtSLDiBaryonAmp>( *ffModel_ );

    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Wrong baryon spin type in EvtBToDiBaryonlnupQCD model. "
            << "Expected spin type " << EvtSpinType::DIRAC << " or "
            << EvtSpinType::RARITASCHWINGER << ", found spin types "
            << baryon1Type << " and " << baryon2Type << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << std::endl;
        ::abort();
    }
}

void EvtBToDiBaryonlnupQCD::initProbMax()
{
    // Set maximum prob using dec file parameter if present
    if ( getNArg() == 7 ) {
        setProbMax( getArg( 6 ) );

    } else {
        // Default probability for the B -> p p l nu mode, where l = e, mu or tau
        setProbMax( 3.0e6 );

        // Specific decay modes, where we have one proton plus a second
        // baryon that can be any (excited) state. They all have lower
        // maximum probabilities compared to the default pp mode in order
        // to improve accept/reject generation efficiency
        static EvtIdSet BMesons( "B-", "B+" );
        static EvtIdSet Delta( "Delta+", "anti-Delta-" );
        static EvtIdSet LambdaC( "Lambda_c+", "anti-Lambda_c-" );
        static EvtIdSet LambdaC1( "Lambda_c(2593)+", "anti-Lambda_c(2593)-" );
        static EvtIdSet LambdaC2( "Lambda_c(2625)+", "anti-Lambda_c(2625)-" );
        static EvtIdSet N1440( "N(1440)+", "anti-N(1440)-" );
        static EvtIdSet N1520( "N(1520)+", "anti-N(1520)-" );
        static EvtIdSet N1535( "N(1535)+", "anti-N(1535)-" );
        static EvtIdSet N1650( "N(1650)+", "anti-N(1650)-" );
        static EvtIdSet N1700( "N(1700)+", "anti-N(1700)-" );
        static EvtIdSet N1710( "N(1710)+", "anti-N(1710)-" );
        static EvtIdSet N1720( "N(1720)+", "anti-N(1720)-" );

        EvtId parId = getParentId();
        EvtId bar1Id = getDaug( 0 );
        EvtId bar2Id = getDaug( 1 );

        // These probabilties are sensitive to the sub-decay modes of the excited baryon states,
        // which limit the available phase space and allows for events to be generated within the
        // 10,000 event trial limit. Otherwise the amplitude varies too much (by more than a factor
        // of a million) and events fail to be generated correctly. In case of problems, specify
        // the maximum probability by passing an extra 7th model parameter
        if ( BMesons.contains( parId ) ) {
            if ( Delta.contains( bar1Id ) || Delta.contains( bar2Id ) ) {
                // Delta
                setProbMax( 1e7 );

            } else if ( LambdaC.contains( bar1Id ) ||
                        LambdaC.contains( bar2Id ) ) {
                // Lambda_c+
                setProbMax( 1000.0 );

            } else if ( LambdaC1.contains( bar1Id ) ||
                        LambdaC1.contains( bar2Id ) ) {
                // Lambda_c+(2593)
                setProbMax( 200.0 );

            } else if ( LambdaC2.contains( bar1Id ) ||
                        LambdaC2.contains( bar2Id ) ) {
                // Lambda_c+(2625)
                setProbMax( 500.0 );

            } else if ( N1440.contains( bar1Id ) || N1440.contains( bar2Id ) ) {
                // N(1440)
                setProbMax( 8e5 );

            } else if ( N1520.contains( bar1Id ) || N1520.contains( bar2Id ) ) {
                // N(1520)
                setProbMax( 8e6 );

            } else if ( N1535.contains( bar1Id ) || N1535.contains( bar2Id ) ) {
                // N(1535)
                setProbMax( 8e5 );

            } else if ( N1650.contains( bar1Id ) || N1650.contains( bar2Id ) ) {
                // N(1650)
                setProbMax( 8e5 );

            } else if ( N1700.contains( bar1Id ) || N1700.contains( bar2Id ) ) {
                // N(1700)
                setProbMax( 4e6 );

            } else if ( N1710.contains( bar1Id ) || N1710.contains( bar2Id ) ) {
                // N(1710)
                setProbMax( 5e5 );

            } else if ( N1720.contains( bar1Id ) || N1720.contains( bar2Id ) ) {
                // N(1720)
                setProbMax( 4e6 );

            }    // Baryon combinations

        }    // B parent

    }    // Specific modes
}
