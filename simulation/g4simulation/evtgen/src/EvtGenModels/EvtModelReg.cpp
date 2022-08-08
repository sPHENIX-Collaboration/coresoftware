
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

#include "EvtGenModels/EvtModelReg.hh"

#include "EvtGenBase/EvtModel.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenModels/EvtBBScalar.hh"
#include "EvtGenModels/EvtBHadronic.hh"
#include "EvtGenModels/EvtBLLNuL.hh"
#include "EvtGenModels/EvtBTo3piCP.hh"
#include "EvtGenModels/EvtBTo4piCP.hh"
#include "EvtGenModels/EvtBToDDalitzCPK.hh"
#include "EvtGenModels/EvtBToDiBaryonlnupQCD.hh"
#include "EvtGenModels/EvtBToKpipiCP.hh"
#include "EvtGenModels/EvtBToPlnuBK.hh"
#include "EvtGenModels/EvtBToVlnuBall.hh"
#include "EvtGenModels/EvtBToXElNu.hh"
#include "EvtGenModels/EvtBaryonPCR.hh"
#include "EvtGenModels/EvtBcBsNPi.hh"
#include "EvtGenModels/EvtBcBsStarNPi.hh"
#include "EvtGenModels/EvtBcPsiNPi.hh"
#include "EvtGenModels/EvtBcSMuNu.hh"
#include "EvtGenModels/EvtBcTMuNu.hh"
#include "EvtGenModels/EvtBcToNPi.hh"
#include "EvtGenModels/EvtBcVHad.hh"
#include "EvtGenModels/EvtBcVMuNu.hh"
#include "EvtGenModels/EvtBcVNpi.hh"
#include "EvtGenModels/EvtBsMuMuKK.hh"
#include "EvtGenModels/EvtBsquark.hh"
#include "EvtGenModels/EvtBto2piCPiso.hh"
#include "EvtGenModels/EvtBtoKD3P.hh"
#include "EvtGenModels/EvtBtoKpiCPiso.hh"
#include "EvtGenModels/EvtBtoXsEtap.hh"
#include "EvtGenModels/EvtBtoXsgamma.hh"
#include "EvtGenModels/EvtBtoXsll.hh"
#include "EvtGenModels/EvtCBTo3piMPP.hh"
#include "EvtGenModels/EvtCBTo3piP00.hh"
#include "EvtGenModels/EvtD0gammaDalitz.hh"
#include "EvtGenModels/EvtD0mixDalitz.hh"
#include "EvtGenModels/EvtDDalitz.hh"
#include "EvtGenModels/EvtDMix.hh"
#include "EvtGenModels/EvtDToKpienu.hh"
#include "EvtGenModels/EvtEtaDalitz.hh"
#include "EvtGenModels/EvtEtaLLPiPi.hh"
#include "EvtGenModels/EvtFlatQ2.hh"
#include "EvtGenModels/EvtFlatSqDalitz.hh"
#include "EvtGenModels/EvtFourBodyPhsp.hh"
#include "EvtGenModels/EvtGenericDalitz.hh"
#include "EvtGenModels/EvtGoityRoberts.hh"
#include "EvtGenModels/EvtHQET.hh"
#include "EvtGenModels/EvtHQET2.hh"
#include "EvtGenModels/EvtHelAmp.hh"
#include "EvtGenModels/EvtHypNonLepton.hh"
#include "EvtGenModels/EvtISGW.hh"
#include "EvtGenModels/EvtISGW2.hh"
#include "EvtGenModels/EvtKKLambdaC.hh"
#include "EvtGenModels/EvtKStopizmumu.hh"
#include "EvtGenModels/EvtKstarnunu.hh"
#include "EvtGenModels/EvtKstarstargamma.hh"
#include "EvtGenModels/EvtLNuGamma.hh"
#include "EvtGenModels/EvtLambdaB2LambdaV.hh"
#include "EvtGenModels/EvtLambdaP_BarGamma.hh"
#include "EvtGenModels/EvtLambdacPHH.hh"
#include "EvtGenModels/EvtLb2Baryonlnu.hh"
#include "EvtGenModels/EvtLb2Lll.hh"
#include "EvtGenModels/EvtLb2plnuLCSR.hh"
#include "EvtGenModels/EvtLb2plnuLQCD.hh"
#include "EvtGenModels/EvtMelikhov.hh"
#include "EvtGenModels/EvtMultibody.hh"
#include "EvtGenModels/EvtOmegaDalitz.hh"
#include "EvtGenModels/EvtPVVCPLH.hh"
#include "EvtGenModels/EvtPartWave.hh"
#include "EvtGenModels/EvtPhiDalitz.hh"
#include "EvtGenModels/EvtPhsp.hh"
#include "EvtGenModels/EvtPhspDecaytimeCut.hh"
#include "EvtGenModels/EvtPhspFlatLifetime.hh"
#include "EvtGenModels/EvtPi0Dalitz.hh"
#include "EvtGenModels/EvtPropSLPole.hh"
#include "EvtGenModels/EvtPsi2JpsiPiPi.hh"
#include "EvtGenModels/EvtPto3P.hh"
#include "EvtGenModels/EvtRareLbToLll.hh"
#include "EvtGenModels/EvtSLBKPole.hh"
#include "EvtGenModels/EvtSLN.hh"
#include "EvtGenModels/EvtSLPole.hh"
#include "EvtGenModels/EvtSSDCP.hh"
#include "EvtGenModels/EvtSSD_DirectCP.hh"
#include "EvtGenModels/EvtSSSCP.hh"
#include "EvtGenModels/EvtSSSCPT.hh"
#include "EvtGenModels/EvtSSSCPpng.hh"
#include "EvtGenModels/EvtSTS.hh"
#include "EvtGenModels/EvtSTSCP.hh"
#include "EvtGenModels/EvtSVP.hh"
#include "EvtGenModels/EvtSVPCP.hh"
#include "EvtGenModels/EvtSVPHelAmp.hh"
#include "EvtGenModels/EvtSVPHelCPMix.hh"
#include "EvtGenModels/EvtSVS.hh"
#include "EvtGenModels/EvtSVSCP.hh"
#include "EvtGenModels/EvtSVSCPLH.hh"
#include "EvtGenModels/EvtSVSCPiso.hh"
#include "EvtGenModels/EvtSVSNONCPEIGEN.hh"
#include "EvtGenModels/EvtSVVCP.hh"
#include "EvtGenModels/EvtSVVCPLH.hh"
#include "EvtGenModels/EvtSVVHelAmp.hh"
#include "EvtGenModels/EvtSVVHelCPMix.hh"
#include "EvtGenModels/EvtSVVNONCPEIGEN.hh"
#include "EvtGenModels/EvtSingleParticle.hh"
#include "EvtGenModels/EvtSll.hh"
#include "EvtGenModels/EvtTSS.hh"
#include "EvtGenModels/EvtTVP.hh"
#include "EvtGenModels/EvtTVSPwave.hh"
#include "EvtGenModels/EvtTauHadnu.hh"
#include "EvtGenModels/EvtTauScalarnu.hh"
#include "EvtGenModels/EvtTauVectornu.hh"
#include "EvtGenModels/EvtTaulnunu.hh"
#include "EvtGenModels/EvtThreeBodyPhsp.hh"
#include "EvtGenModels/EvtVPHOtoVISRHi.hh"
#include "EvtGenModels/EvtVSPPwave.hh"
#include "EvtGenModels/EvtVSS.hh"
#include "EvtGenModels/EvtVSSBMixCPT.hh"
#include "EvtGenModels/EvtVSSMix.hh"
#include "EvtGenModels/EvtVVP.hh"
#include "EvtGenModels/EvtVVPIPI_WEIGHTED.hh"
#include "EvtGenModels/EvtVVSPwave.hh"
#include "EvtGenModels/EvtVVpipi.hh"
#include "EvtGenModels/EvtVectorIsr.hh"
#include "EvtGenModels/EvtVll.hh"
#include "EvtGenModels/EvtVtoSll.hh"
#include "EvtGenModels/EvtVub.hh"
#include "EvtGenModels/EvtVubBLNP.hh"
#include "EvtGenModels/EvtVubBLNPHybrid.hh"
#include "EvtGenModels/EvtVubHybrid.hh"
#include "EvtGenModels/EvtVubNLO.hh"
#include "EvtGenModels/EvtXPsiGamma.hh"
#include "EvtGenModels/EvtY3SToY1SpipiMoxhay.hh"
#include "EvtGenModels/EvtYmSToYnSpipiCLEO.hh"
#include "EvtGenModels/EvtbTosllAli.hh"
#include "EvtGenModels/EvtbTosllBall.hh"
#include "EvtGenModels/EvtbTosllMS.hh"
#include "EvtGenModels/EvtbTosllMSExt.hh"
#include "EvtGenModels/Evtbs2llGammaISRFSR.hh"
#include "EvtGenModels/Evtbs2llGammaMNT.hh"
#include "EvtGenModels/EvtbsToLLLL.hh"
#include "EvtGenModels/EvtbsToLLLLHyperCP.hh"

#include <assert.h>
#include <ctype.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>

using std::cout;
using std::endl;
using std::fstream;

EvtModelReg::EvtModelReg( const std::list<EvtDecayBase*>* extraModels )
{
    EvtModel& modelist = EvtModel::instance();

    if ( extraModels ) {
        for ( std::list<EvtDecayBase*>::const_iterator it = extraModels->begin();
              it != extraModels->end(); ++it ) {
            modelist.registerModel( *it );
        }
    }

    modelist.registerModel( new EvtBBScalar );
    modelist.registerModel( new EvtLambdaP_BarGamma );
    modelist.registerModel( new EvtFlatQ2 );
    modelist.registerModel( new EvtTauHadnu );
    modelist.registerModel( new EvtTauVectornu );
    modelist.registerModel( new EvtVVP );
    modelist.registerModel( new EvtSLN );
    modelist.registerModel( new EvtISGW2 );
    modelist.registerModel( new EvtMelikhov );
    modelist.registerModel( new EvtSLPole );
    modelist.registerModel( new EvtPropSLPole );
    modelist.registerModel( new EvtSLBKPole );
    modelist.registerModel( new EvtHQET );
    modelist.registerModel( new EvtHQET2 );
    modelist.registerModel( new EvtISGW );
    modelist.registerModel( new EvtBHadronic );
    modelist.registerModel( new EvtVSS );
    modelist.registerModel( new EvtVSSMix );
    modelist.registerModel( new EvtVSSBMixCPT );
    modelist.registerModel( new EvtVSPPwave );
    modelist.registerModel( new EvtGoityRoberts );
    modelist.registerModel( new EvtSVS );
    modelist.registerModel( new EvtTSS );
    modelist.registerModel( new EvtTVSPwave );
    modelist.registerModel( new EvtSVVHelAmp );
    modelist.registerModel( new EvtSVPHelAmp );
    modelist.registerModel( new EvtSVPCP );
    modelist.registerModel( new EvtVVSPwave );
    modelist.registerModel( new EvtDDalitz );
    modelist.registerModel( new EvtOmegaDalitz );
    modelist.registerModel( new EvtEtaDalitz );
    modelist.registerModel( new EvtPhsp );
    modelist.registerModel( new EvtPhspDecaytimeCut );
    modelist.registerModel( new EvtBtoXsgamma );
    modelist.registerModel( new EvtBtoXsll );
    modelist.registerModel( new EvtBtoXsEtap );
    modelist.registerModel( new EvtSSSCP );
    modelist.registerModel( new EvtSSSCPpng );
    modelist.registerModel( new EvtSTSCP );
    modelist.registerModel( new EvtSTS );
    modelist.registerModel( new EvtSSSCPT );
    modelist.registerModel( new EvtSVSCP );
    modelist.registerModel( new EvtSSDCP );
    modelist.registerModel( new EvtSVSNONCPEIGEN );
    modelist.registerModel( new EvtSVVNONCPEIGEN );
    modelist.registerModel( new EvtSVVCP );
    modelist.registerModel( new EvtSVVCPLH );
    modelist.registerModel( new EvtSVSCPLH );
    modelist.registerModel( new EvtSll );
    modelist.registerModel( new EvtVll );
    modelist.registerModel( new EvtTaulnunu );
    modelist.registerModel( new EvtTauScalarnu );
    modelist.registerModel( new EvtKstarnunu );
    modelist.registerModel( new EvtbTosllBall );
    modelist.registerModel( new EvtBto2piCPiso );
    modelist.registerModel( new EvtBtoKpiCPiso );
    modelist.registerModel( new EvtSVSCPiso );
    modelist.registerModel( new EvtSingleParticle );
    modelist.registerModel( new EvtVectorIsr );
    modelist.registerModel( new EvtPi0Dalitz );
    modelist.registerModel( new EvtHelAmp );
    modelist.registerModel( new EvtPartWave );
    modelist.registerModel( new EvtVVpipi );
    modelist.registerModel( new EvtY3SToY1SpipiMoxhay );
    modelist.registerModel( new EvtYmSToYnSpipiCLEO );
    modelist.registerModel( new EvtBsquark );
    modelist.registerModel( new EvtPhiDalitz );
    modelist.registerModel( new EvtBToPlnuBK );
    modelist.registerModel( new EvtBToVlnuBall );
    modelist.registerModel( new EvtVVPIPI_WEIGHTED );
    modelist.registerModel( new EvtVPHOtoVISRHi );

    modelist.registerModel( new EvtBTo4piCP );
    modelist.registerModel( new EvtBTo3piCP );
    modelist.registerModel( new EvtCBTo3piP00 );
    modelist.registerModel( new EvtCBTo3piMPP );
    modelist.registerModel( new EvtBToKpipiCP );

    modelist.registerModel( new EvtLb2Lll );
    modelist.registerModel( new EvtRareLbToLll );
    modelist.registerModel( new EvtHypNonLepton );
    modelist.registerModel( new EvtSVVHelCPMix );
    modelist.registerModel( new EvtSVPHelCPMix );

    modelist.registerModel( new EvtLNuGamma );
    modelist.registerModel( new EvtKstarstargamma );

    modelist.registerModel( new EvtVub );

    modelist.registerModel( new EvtVubHybrid );
    modelist.registerModel( new EvtVubNLO );
    modelist.registerModel( new EvtVubBLNP );
    modelist.registerModel( new EvtVubBLNPHybrid );

    modelist.registerModel( new EvtPto3P );
    modelist.registerModel( new EvtBtoKD3P );
    modelist.registerModel( new EvtKKLambdaC );
    modelist.registerModel( new EvtMultibody );
    modelist.registerModel( new EvtDMix );
    modelist.registerModel( new EvtD0mixDalitz );
    modelist.registerModel( new EvtD0gammaDalitz );

    modelist.registerModel( new EvtbTosllAli );
    modelist.registerModel( new EvtBaryonPCR );

    modelist.registerModel( new EvtBToDDalitzCPK );
    modelist.registerModel( new EvtLambdaB2LambdaV );
    modelist.registerModel( new EvtLambda2PPiForLambdaB2LambdaV );
    modelist.registerModel( new EvtV2VpVmForLambdaB2LambdaV );
    modelist.registerModel( new EvtPVVCPLH );
    modelist.registerModel( new EvtSSD_DirectCP );

    modelist.registerModel( new EvtBcToNPi( true ) );    // true = print author info
    modelist.registerModel( new EvtBcPsiNPi );
    modelist.registerModel( new EvtBcBsNPi );
    modelist.registerModel( new EvtBcBsStarNPi );

    modelist.registerModel( new EvtBcSMuNu );
    modelist.registerModel( new EvtBcVMuNu );
    modelist.registerModel( new EvtBcTMuNu );
    modelist.registerModel( new EvtBcVNpi );
    modelist.registerModel( new EvtSVP );
    modelist.registerModel( new EvtTVP );
    modelist.registerModel( new EvtXPsiGamma );

    modelist.registerModel( new EvtbsToLLLL );
    modelist.registerModel( new EvtbsToLLLLHyperCP );
    modelist.registerModel( new EvtBLLNuL );

    modelist.registerModel( new EvtKStopizmumu );
    modelist.registerModel( new EvtVtoSll );

    modelist.registerModel( new EvtBsMuMuKK );
    modelist.registerModel( new EvtGenericDalitz() );

    modelist.registerModel( new EvtBcVHad );

    modelist.registerModel( new Evtbs2llGammaMNT );
    modelist.registerModel( new Evtbs2llGammaISRFSR );
    modelist.registerModel( new EvtbTosllMS );
    modelist.registerModel( new EvtbTosllMSExt );

    modelist.registerModel( new EvtLb2plnuLQCD );
    modelist.registerModel( new EvtLb2plnuLCSR );
    modelist.registerModel( new EvtLb2Baryonlnu );

    modelist.registerModel( new EvtBToDiBaryonlnupQCD );

    modelist.registerModel( new EvtFlatSqDalitz );
    modelist.registerModel( new EvtPhspFlatLifetime );

    modelist.registerModel( new EvtLambdacPHH );

    modelist.registerModel( new EvtDToKpienu );
    modelist.registerModel( new EvtPsi2JpsiPiPi );

    modelist.registerModel( new EvtThreeBodyPhsp );
    modelist.registerModel( new EvtFourBodyPhsp );
    modelist.registerModel( new EvtEtaLLPiPi );

    modelist.registerModel( new EvtBToXElNu );
}
