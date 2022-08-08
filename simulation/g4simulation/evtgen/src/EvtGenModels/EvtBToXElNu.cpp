
/***********************************************************************
* Copyright 1998-2021 CERN for the benefit of the EvtGen authors       *
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

#include "EvtGenModels/EvtBToXElNu.hh"

#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicVectorAmp.hh"

#include "EvtGenModels/EvtBCLFF.hh"
#include "EvtGenModels/EvtBGLFF.hh"

#include <cstdlib>
#include <string>

using std::endl;

std::string EvtBToXElNu::getName()
{
    return "BTOXELNU";
}

EvtDecayBase* EvtBToXElNu::clone()
{
    return new EvtBToXElNu;
}

void EvtBToXElNu::decay( EvtParticle* p )
{
    p->initializePhaseSpace( getNDaug(), getDaugs() );
    m_calcamp->CalcAmp( p, _amp2, m_ffmodel.get() );
}

void EvtBToXElNu::initProbMax()
{
    EvtId parnum, mesnum, lnum, nunum;

    parnum = getParentId();
    mesnum = getDaug( 0 );
    lnum = getDaug( 1 );
    nunum = getDaug( 2 );

    const double mymaxprob = m_calcamp->CalcMaxProb( parnum, mesnum, lnum,
                                                     nunum, m_ffmodel.get() );

    setProbMax( mymaxprob );
}

void EvtBToXElNu::init()
{
    checkNDaug( 3 );

    // We expect the parent to be a scalar
    // and the daughters to be X lepton neutrino
    checkSpinParent( EvtSpinType::SCALAR );

    checkSpinDaughter( 1, EvtSpinType::DIRAC );
    checkSpinDaughter( 2, EvtSpinType::NEUTRINO );

    EvtSpinType::spintype d1type = EvtPDL::getSpinType( getDaug( 0 ) );
    const std::string model = getArgStr( 0 );
    EvtGenReport( EVTGEN_INFO, "EvtGen" )
        << "Using " << model << " EvtGen FF model" << endl;

    // -- BGL -- BGL -- BGL -- BGL -- BGL -- BGL -- BGL -- BGL -- BGL -- BGL -- BGL -- BGL -- BGL
    if ( model == "BGL" ) {
        const EvtId lnum = getDaug( 1 );
        if ( lnum == EvtPDL::getId( "tau-" ) || lnum == EvtPDL::getId( "tau+" ) ) {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "The current BGL model should not be used for taus" << endl;
            ::abort();
        }

        if ( d1type == EvtSpinType::SCALAR ) {
            if ( getNArg() == 9 ) {
                m_ffmodel = std::make_unique<EvtBGLFF>(
                    getArg( 1 ), getArg( 2 ), getArg( 3 ), getArg( 4 ),
                    getArg( 5 ), getArg( 6 ), getArg( 7 ), getArg( 8 ) );
                m_calcamp = std::make_unique<EvtSemiLeptonicScalarAmp>();
            } else {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "BGL (N=3) model for scalar meson daughters needs 8 arguments. Sorry."
                    << endl;
                ::abort();
            }
        } else if ( d1type == EvtSpinType::VECTOR ) {
            if ( getNArg() == 7 ) {
                m_ffmodel = std::make_unique<EvtBGLFF>( getArg( 1 ),
                                                        getArg( 2 ), getArg( 3 ),
                                                        getArg( 4 ), getArg( 5 ),
                                                        getArg( 6 ) );
                m_calcamp = std::make_unique<EvtSemiLeptonicVectorAmp>();
            } else {
                EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                    << "BGL (N=3) model for vector meson daughters needs 6 arguments. Sorry."
                    << endl;
                ::abort();
            }
        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "Only Scalar and Vector models implemented. Sorry." << endl;
            ::abort();
        }
        // -- BCL -- BCL -- BCL -- BCL -- BCL -- BCL -- BCL -- BCL -- BCL -- BCL -- BCL -- BCL -- BCL
    } else if ( model == "BCL" ) {
        // -- Need to subtract 1 as the first argument is the model name (BCL)
        const int numArgs = getNArg() - 1;
        double* args = getArgs();
        args = &args[1];

        m_ffmodel = std::make_unique<EvtBCLFF>( numArgs, args );
        if ( d1type == EvtSpinType::SCALAR ) {
            m_calcamp = std::make_unique<EvtSemiLeptonicScalarAmp>();
        } else if ( d1type == EvtSpinType::VECTOR ) {
            m_calcamp = std::make_unique<EvtSemiLeptonicVectorAmp>();
        } else {
            EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                << "BCL model handles currently only scalar and vector meson daughters. Sorry."
                << endl;
            ::abort();
        }
    } else {
        EvtGenReport( EVTGEN_INFO, "EvtGen" )
            << "  Unknown form-factor model, valid options are BGL, BCL"
            << std::endl;
        ::abort();
    }
}
