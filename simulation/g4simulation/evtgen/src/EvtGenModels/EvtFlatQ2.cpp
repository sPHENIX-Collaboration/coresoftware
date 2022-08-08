
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

#include "EvtGenModels/EvtFlatQ2.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"

#include <fstream>
#include <string>

double lambda( double q, double m1, double m2 )
{
    double L( 1.0 );
    double mSum = m1 + m2;
    double mDiff = m1 - m2;
    double qSq = q * q;

    if ( qSq > 0.0 ) {
        double prodTerm = ( qSq - mSum * mSum ) * ( qSq - mDiff * mDiff );

        if ( prodTerm > 0.0 ) {
            L = sqrt( prodTerm ) / qSq;
        }
    }

    return L;
}

std::string EvtFlatQ2::getName()
{
    return "FLATQ2";
}

EvtDecayBase* EvtFlatQ2::clone()
{
    return new EvtFlatQ2;
}

void EvtFlatQ2::initProbMax()
{
    setProbMax( 100 );
}

void EvtFlatQ2::init()
{
    // check that there are 3 daughters
    checkNDaug( 3 );

    // We expect B -> X lepton lepton events
    checkSpinParent( EvtSpinType::SCALAR );

    EvtSpinType::spintype d1type = EvtPDL::getSpinType( getDaug( 1 ) );
    EvtSpinType::spintype d2type = EvtPDL::getSpinType( getDaug( 2 ) );

    if ( !( d1type == EvtSpinType::DIRAC || d1type == EvtSpinType::NEUTRINO ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtFlatQ2 expects 2nd daughter to "
            << "be a lepton" << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << std::endl;
        ::abort();
    }

    if ( !( d2type == EvtSpinType::DIRAC || d2type == EvtSpinType::NEUTRINO ) ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtFlatQ2 expects 3rd daughter to "
            << "be a lepton" << std::endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Will terminate execution!" << std::endl;
        ::abort();
    }

    // Specify if we want to use the phase space factor
    _usePhsp = false;
    if ( getNArg() > 0 ) {
        if ( getArg( 0 ) != 0 ) {
            _usePhsp = true;
        }
    }

    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "EvtFlatQ2 usePhsp = " << int( _usePhsp ) << std::endl;
}

void EvtFlatQ2::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtVector4R p4Xu = p->getDaug( 0 )->getP4();

    EvtVector4R p4ell1 = p->getDaug( 1 )->getP4();
    EvtVector4R p4ell2 = p->getDaug( 2 )->getP4();

    double pXu_x2 = p4Xu.get( 1 ) * p4Xu.get( 1 );
    double pXu_y2 = p4Xu.get( 2 ) * p4Xu.get( 2 );
    double pXu_z2 = p4Xu.get( 3 ) * p4Xu.get( 3 );
    double pXu = sqrt( pXu_x2 + pXu_y2 + pXu_z2 );
    double prob( 0.0 );
    if ( fabs( pXu ) > 0.0 ) {
        prob = 1 / pXu;
    }

    // Include the phase space factor if requested
    if ( _usePhsp ) {
        // Invariant mass of lepton pair
        double q = ( p4ell1 + p4ell2 ).mass();
        // Rest masses of the leptons
        double m1 = p4ell1.mass();
        double m2 = p4ell2.mass();
        // Phase space factor, which includes the various square roots
        double Lambda = lambda( q, m1, m2 );
        if ( Lambda > 0.0 ) {
            prob = prob / Lambda;
        }
    }

    if ( pXu > 0.01 ) {
        setProb( prob );
    }

    return;
}
