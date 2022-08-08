
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

#include "EvtGenModels/EvtBcVNpi.hh"

#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParser.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include "EvtGenModels/EvtTauHadnu.hh"
#include "EvtGenModels/EvtWnPi.hh"

#include <ctype.h>
#include <stdlib.h>

std::string EvtBcVNpi::getName()
{
    return "BC_VNPI";
}

EvtDecayBase* EvtBcVNpi::clone()
{
    return new EvtBcVNpi;
}

//======================================================
void EvtBcVNpi::init()
{
    //cout<<"BcVNpi::init()"<<endl;

    checkNArg( 1 );
    checkSpinParent( EvtSpinType::SCALAR );
    checkSpinDaughter( 0, EvtSpinType::VECTOR );
    for ( int i = 1; i <= ( getNDaug() - 1 ); i++ ) {
        checkSpinDaughter( i, EvtSpinType::SCALAR );
    };

    if ( getNDaug() < 2 || getNDaug() > 6 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Have not yet implemented this final state in BcVNpi model"
            << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Ndaug=" << getNDaug() << endl;
        for ( int id = 0; id < ( getNDaug() - 1 ); id++ )
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Daug " << id << " " << EvtPDL::name( getDaug( id ) ).c_str()
                << endl;
        return;
    }

    //     for(int i=0; i<getNDaug(); i++)
    //       cout<<"BcVNpi::init \t\t daughter "<<i<<" : "<<getDaug(i).getId()<<"   "<<EvtPDL::name(getDaug(i)).c_str()<<endl;

    idVector = getDaug( 0 ).getId();
    whichfit = int( getArg( 0 ) + 0.1 );
    //     cout<<"BcVNpi: whichfit ="<<whichfit<<"  idVector="<<idVector<<endl;
    ffmodel = std::make_unique<EvtBCVFF>( idVector, whichfit );

    wcurr = std::make_unique<EvtWnPi>();

    nCall = 0;
}

//======================================================
void EvtBcVNpi::initProbMax()
{
    //     cout<<"BcVNpi::initProbMax()"<<endl;
    if ( idVector == EvtPDL::getId( "J/psi" ).getId() && whichfit == 1 &&
         getNDaug() == 6 )
        setProbMax( 720000. );
    else if ( idVector == EvtPDL::getId( "J/psi" ).getId() && whichfit == 2 &&
              getNDaug() == 6 )
        setProbMax( 471817. );
    else if ( idVector == EvtPDL::getId( "J/psi" ).getId() && whichfit == 1 &&
              getNDaug() == 4 )
        setProbMax( 42000. );
    else if ( idVector == EvtPDL::getId( "J/psi" ).getId() && whichfit == 2 &&
              getNDaug() == 4 )
        setProbMax( 16000. );

    else if ( idVector == EvtPDL::getId( "psi(2S)" ).getId() && whichfit == 1 &&
              getNDaug() == 4 )
        setProbMax( 1200. );
    else if ( idVector == EvtPDL::getId( "psi(2S)" ).getId() && whichfit == 2 &&
              getNDaug() == 4 )
        setProbMax( 2600. );
    else if ( idVector == EvtPDL::getId( "psi(2S)" ).getId() && whichfit == 1 &&
              getNDaug() == 6 )
        setProbMax( 40000. );
    else if ( idVector == EvtPDL::getId( "psi(2S)" ).getId() && whichfit == 2 &&
              getNDaug() == 6 )
        setProbMax( 30000. );
}

//======================================================
void EvtBcVNpi::decay( EvtParticle* root_particle )
{
    ++nCall;
    //     cout<<"BcVNpi::decay()"<<endl;
    root_particle->initializePhaseSpace( getNDaug(), getDaugs() );

    EvtVector4R p4b( root_particle->mass(), 0., 0., 0. ),    // Bc momentum
        p4meson = root_particle->getDaug( 0 )->getP4(),      // J/psi momenta
        Q = p4b - p4meson;
    double Q2 = Q.mass2();

    // check pi-mesons and calculate hadronic current
    EvtVector4C hardCur;
    //     bool foundHadCurr=false;
    if ( getNDaug() == 2 ) {
        hardCur = wcurr->WCurrent( root_particle->getDaug( 1 )->getP4() );
        //       foundHadCurr=true;
    } else if ( getNDaug() == 3 ) {
        hardCur = wcurr->WCurrent( root_particle->getDaug( 1 )->getP4(),
                                   root_particle->getDaug( 2 )->getP4() );
        //       foundHadCurr=true;
    } else if ( getNDaug() == 4 ) {
        hardCur = wcurr->WCurrent( root_particle->getDaug( 1 )->getP4(),
                                   root_particle->getDaug( 2 )->getP4(),
                                   root_particle->getDaug( 3 )->getP4() );
        //       foundHadCurr=true;
    } else if ( getNDaug() ==
                6 )    // Bc -> psi pi+ pi+ pi- pi- pi+ from [Kuhn, Was, hep-ph/0602162
    {
        hardCur = wcurr->WCurrent( root_particle->getDaug( 1 )->getP4(),
                                   root_particle->getDaug( 2 )->getP4(),
                                   root_particle->getDaug( 3 )->getP4(),
                                   root_particle->getDaug( 4 )->getP4(),
                                   root_particle->getDaug( 5 )->getP4() );
        // 		foundHadCurr=true;
    } else {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "Have not yet implemented this final state in BCNPI model"
            << endl;
        EvtGenReport( EVTGEN_ERROR, "EvtGen" ) << "Ndaug=" << getNDaug() << endl;
        int id;
        for ( id = 0; id < ( getNDaug() - 1 ); id++ )
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Daug " << id << " " << EvtPDL::name( getDaug( id ) ).c_str()
                << endl;
        ::abort();
    };

    // calculate Bc -> V W form-factors
    double a1f, a2f, vf, a0f;
    double m_meson = root_particle->getDaug( 0 )->mass();
    double m_b = root_particle->mass();
    ffmodel->getvectorff( root_particle->getId(),
                          root_particle->getDaug( 0 )->getId(), Q2, m_meson,
                          &a1f, &a2f, &vf, &a0f );
    double a3f = ( ( m_b + m_meson ) / ( 2.0 * m_meson ) ) * a1f -
                 ( ( m_b - m_meson ) / ( 2.0 * m_meson ) ) * a2f;

    // calculate Bc -> V W current
    EvtTensor4C H;
    H = a1f * ( m_b + m_meson ) * EvtTensor4C::g();
    H.addDirProd( ( -a2f / ( m_b + m_meson ) ) * p4b, p4b + p4meson );
    H += EvtComplex( 0.0, vf / ( m_b + m_meson ) ) *
         dual( EvtGenFunctions::directProd( p4meson + p4b, p4b - p4meson ) );
    H.addDirProd( ( a0f - a3f ) * 2.0 * ( m_meson / Q2 ) * p4b, p4b - p4meson );
    EvtVector4C Heps = H.cont2( hardCur );

    for ( int i = 0; i < 4; i++ ) {
        EvtVector4C eps = root_particle->getDaug( 0 )
                              ->epsParent( i )
                              .conj();    // psi-meson polarization vector
        EvtComplex amp = eps * Heps;
        vertex( i, amp );
    };
}
