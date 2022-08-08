
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

#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"

#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"
using std::endl;

void EvtSemiLeptonicVectorAmp::CalcAmp( EvtParticle* parent, EvtAmp& amp,
                                        EvtSemiLeptonicFF* FormFactors )
{
    static EvtId EM = EvtPDL::getId( "e-" );
    static EvtId MUM = EvtPDL::getId( "mu-" );
    static EvtId TAUM = EvtPDL::getId( "tau-" );
    static EvtId EP = EvtPDL::getId( "e+" );
    static EvtId MUP = EvtPDL::getId( "mu+" );
    static EvtId TAUP = EvtPDL::getId( "tau+" );

    static EvtId D0 = EvtPDL::getId( "D0" );
    static EvtId D0B = EvtPDL::getId( "anti-D0" );
    static EvtId DP = EvtPDL::getId( "D+" );
    static EvtId DM = EvtPDL::getId( "D-" );
    static EvtId DSM = EvtPDL::getId( "D_s-" );
    static EvtId DSP = EvtPDL::getId( "D_s+" );

    //Add the lepton and neutrino 4 momenta to find q2

    EvtVector4R q = parent->getDaug( 1 )->getP4() + parent->getDaug( 2 )->getP4();
    double q2 = ( q.mass2() );

    double a1f, a2f, vf, a0f, a3f;
    double m_meson = parent->getDaug( 0 )->mass();

    FormFactors->getvectorff( parent->getId(), parent->getDaug( 0 )->getId(),
                              q2, m_meson, &a1f, &a2f, &vf, &a0f );

    double costhl_flag = 1.0;

    if ( parent->getId() == D0 || parent->getId() == D0B ||
         parent->getId() == DP || parent->getId() == DM ) {
        costhl_flag = -1.0;
    }
    if ( parent->getId() == DSP || parent->getId() == DSM ) {
        costhl_flag = -1.0;
    }
    vf = vf * costhl_flag;

    EvtVector4R p4b;
    p4b.set( parent->mass(), 0.0, 0.0, 0.0 );

    EvtVector4R p4meson = parent->getDaug( 0 )->getP4();

    EvtVector4C l1, l2;

    EvtId l_num = parent->getDaug( 1 )->getId();
    double m_b = parent->mass();

    a3f = ( ( m_b + m_meson ) / ( 2.0 * m_meson ) ) * a1f -
          ( ( m_b - m_meson ) / ( 2.0 * m_meson ) ) * a2f;

    EvtTensor4C tds;
    if ( l_num == EM || l_num == MUM || l_num == TAUM ) {
        tds = a1f * ( m_b + m_meson ) * EvtTensor4C::g();
        tds.addDirProd( ( -a2f / ( m_b + m_meson ) ) * p4b, p4b + p4meson );
        tds += EvtComplex( 0.0, vf / ( m_b + m_meson ) ) *
               dual( EvtGenFunctions::directProd( p4meson + p4b, p4b - p4meson ) );
        tds.addDirProd( ( a0f - a3f ) * 2.0 * ( m_meson / q2 ) * p4b,
                        p4b - p4meson );

        l1 = EvtLeptonVACurrent( parent->getDaug( 1 )->spParent( 0 ),
                                 parent->getDaug( 2 )->spParentNeutrino() );
        l2 = EvtLeptonVACurrent( parent->getDaug( 1 )->spParent( 1 ),
                                 parent->getDaug( 2 )->spParentNeutrino() );
    } else {
        if ( l_num == EP || l_num == MUP || l_num == TAUP ) {
            tds = a1f * ( m_b + m_meson ) * EvtTensor4C::g();
            tds.addDirProd( ( -a2f / ( m_b + m_meson ) ) * p4b, p4b + p4meson );
            tds -= EvtComplex( 0.0, vf / ( m_b + m_meson ) ) *
                   dual( EvtGenFunctions::directProd( p4meson + p4b,
                                                      p4b - p4meson ) );
            tds.addDirProd( ( a0f - a3f ) * 2.0 * ( m_meson / q2 ) * p4b,
                            p4b - p4meson );

            l1 = EvtLeptonVACurrent( parent->getDaug( 2 )->spParentNeutrino(),
                                     parent->getDaug( 1 )->spParent( 0 ) );
            l2 = EvtLeptonVACurrent( parent->getDaug( 2 )->spParentNeutrino(),
                                     parent->getDaug( 1 )->spParent( 1 ) );
        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Wrong lepton number" << endl;
        }
    }

    EvtVector4C et0 = tds.cont1( parent->getDaug( 0 )->epsParent( 0 ).conj() );
    EvtVector4C et1 = tds.cont1( parent->getDaug( 0 )->epsParent( 1 ).conj() );
    EvtVector4C et2 = tds.cont1( parent->getDaug( 0 )->epsParent( 2 ).conj() );

    amp.vertex( 0, 0, l1.cont( et0 ) );
    amp.vertex( 0, 1, l2.cont( et0 ) );

    amp.vertex( 1, 0, l1.cont( et1 ) );
    amp.vertex( 1, 1, l2.cont( et1 ) );

    amp.vertex( 2, 0, l1.cont( et2 ) );
    amp.vertex( 2, 1, l2.cont( et2 ) );

    return;
}
