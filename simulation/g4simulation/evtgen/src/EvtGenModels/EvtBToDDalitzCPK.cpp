
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

#include "EvtGenModels/EvtBToDDalitzCPK.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtVector4C.hh"

//-----------------------------------------------------------------------------
// Implementation file for class : EvtBToDDalitzCPK
// Decay Model for B->D0K with D0->Ks pi+ pi-
// it is just a way to get the arguments...
// Works also for other B->D0K decays...
// 2003-12-08 : Patrick Robbe
//-----------------------------------------------------------------------------

//=============================================================================
// Name of the model
//=============================================================================
std::string EvtBToDDalitzCPK::getName()
{
    return "BTODDALITZCPK";
}
//=============================================================================
// Clone method
//=============================================================================
EvtBToDDalitzCPK* EvtBToDDalitzCPK::clone()
{
    return new EvtBToDDalitzCPK;
}
//=============================================================================
// Initialisation method
//=============================================================================
void EvtBToDDalitzCPK::init()
{
    // Check that there are 3 arguments
    checkNArg( 3 );
    // Check that there are 2 daughters
    checkNDaug( 2 );
    // Check that the particles of the decay are :
    // B+/-   -> D0/bar K+/-
    // B+/-   -> K+/- D0/bar
    // B0/bar -> K*0/bar D0/bar
    // and nothing else ...
    static EvtId BP = EvtPDL::getId( "B+" );
    static EvtId BM = EvtPDL::getId( "B-" );
    static EvtId B0 = EvtPDL::getId( "B0" );
    static EvtId B0B = EvtPDL::getId( "anti-B0" );
    static EvtId KP = EvtPDL::getId( "K+" );
    static EvtId KM = EvtPDL::getId( "K-" );
    static EvtId KS = EvtPDL::getId( "K*0" );
    static EvtId KSB = EvtPDL::getId( "anti-K*0" );
    static EvtId D0 = EvtPDL::getId( "D0" );
    static EvtId D0B = EvtPDL::getId( "anti-D0" );

    _flag = 0;

    EvtId parent = getParentId();
    EvtId d1 = getDaug( 0 );
    EvtId d2 = getDaug( 1 );

    if ( ( ( parent == BP ) || ( parent == BM ) ) &&
         ( ( d1 == D0 ) || ( d1 == D0B ) ) && ( ( d2 == KP ) || ( d2 == KM ) ) ) {
        _flag = 1;
        // PHSP Decay
    } else if ( ( ( parent == BP ) || ( parent == BM ) ) &&
                ( ( d1 == KP ) || ( d1 == KM ) ) &&
                ( ( d2 == D0 ) || ( d2 == D0B ) ) ) {
        _flag = 1;
        // also PHSP decay
    } else if ( ( ( parent == B0 ) || ( parent == B0B ) ) &&
                ( ( d1 == KS ) || ( d1 == KSB ) ) &&
                ( ( d2 == D0 ) || ( d2 == D0B ) ) ) {
        _flag = 2;
        // SVS Decay
    }

    if ( _flag == 0 ) {
        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
            << "EvtBToDDalitzCPK : Invalid mode." << std::endl;
        assert( 0 );
    }
}
//=============================================================================
// Set prob max
//=============================================================================
void EvtBToDDalitzCPK::initProbMax()
{
    if ( _flag == 1 ) {
        // PHSP
        setProbMax( 0. );
    } else if ( _flag == 2 ) {
        // SVS
        setProbMax( 1.0 );
    }
}
//=============================================================================
// decay particle
//=============================================================================
void EvtBToDDalitzCPK::decay( EvtParticle* p )
{
    if ( _flag == 1 ) {
        // PHSP
        p->initializePhaseSpace( getNDaug(), getDaugs() );
        vertex( 0. );
    } else if ( _flag == 2 ) {
        // SVS
        p->initializePhaseSpace( getNDaug(), getDaugs() );

        EvtParticle* v;
        v = p->getDaug( 0 );
        double massv = v->mass();
        EvtVector4R momv = v->getP4();
        EvtVector4R moms = p->getDaug( 1 )->getP4();
        double m_parent = p->mass();
        EvtVector4R p4_parent = momv + moms;

        double norm = massv / ( momv.d3mag() * m_parent );
        p4_parent = norm * p4_parent;
        vertex( 0, p4_parent * ( v->epsParent( 0 ) ) );
        vertex( 1, p4_parent * ( v->epsParent( 1 ) ) );
        vertex( 2, p4_parent * ( v->epsParent( 2 ) ) );
    }
}
