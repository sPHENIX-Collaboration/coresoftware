
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

#include "EvtGenModels/EvtKstarstargamma.hh"

#include "EvtGenBase/EvtEvalHelAmp.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPhotonParticle.hh"
#include "EvtGenBase/EvtPropBreitWignerRel.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtTensorParticle.hh"
#include "EvtGenBase/EvtTwoBodyVertex.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVectorParticle.hh"

#include <stdlib.h>
#include <string>

std::string EvtKstarstargamma::getName()
{
    return "KSTARSTARGAMMA";
}

EvtDecayBase* EvtKstarstargamma::clone()
{
    return new EvtKstarstargamma;
}

void EvtKstarstargamma::init()
{
    // check that there are 0 arguments
    checkNArg( 0 );

    // check that there are 3 daughters
    checkNDaug( 3 );

    // check the parent and daughter spins
    checkSpinParent( EvtSpinType::SCALAR );
    checkSpinDaughter( 0, EvtSpinType::SCALAR );
    checkSpinDaughter( 1, EvtSpinType::SCALAR );
    checkSpinDaughter( 2, EvtSpinType::PHOTON );
}

void EvtKstarstargamma::initProbMax()
{
    //setProbMax(1.0);
}

void EvtKstarstargamma::decay( EvtParticle* /*p*/ )
{
    /*

  The EvtEvalHelAmp is completely broken...

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtParticle* kaon = p->getDaug(0);
  EvtParticle* pion = p->getDaug(1);
  EvtParticle* photon = p->getDaug(2);


  EvtComplexPtrPtr Hd1=new EvtComplexPtr[5];
  Hd1[0]=new EvtComplex[2];
  Hd1[1]=new EvtComplex[2];
  Hd1[2]=new EvtComplex[2];
  Hd1[3]=new EvtComplex[2];
  Hd1[4]=new EvtComplex[2];

  Hd1[0][0]=0.0;
  Hd1[0][1]=0.0;
  Hd1[1][0]=0.0;
  Hd1[1][1]=0.0;
  Hd1[2][0]=0.0;
  Hd1[2][1]=0.0;
  Hd1[3][0]=0.0;
  Hd1[3][1]=1.0;
  Hd1[4][0]=0.0;
  Hd1[4][1]=0.0;

  EvtEvalHelAmp d1(EvtSpinType::SCALAR,EvtSpinType::TENSOR,
		   EvtSpinType::PHOTON,Hd1);

  EvtScalarParticle theB;

  theB.init(p->getId(),p->getP4Restframe());

  EvtVector4R theKstarP4=kaon->getP4()+pion->getP4();

  EvtTensorParticle theKstar;
  theKstar.init(EvtPDL::getId(std::string("K_2*0")),theKstarP4);

  EvtPhotonParticle thePhoton;
  thePhoton.init(EvtPDL::getId(std::string("K_2*0")),photon->getP4());

  theKstar.addDaug(&theB);
  thePhoton.addDaug(&theB);

  EvtAmp amp1;

  d1.evalAmp(&theB,amp1);

  EvtComplexPtrPtr Hd2=new EvtComplexPtr[1];
  Hd2[0]=new EvtComplex[1];

  Hd2[0][0]=1.0;


  EvtEvalHelAmp d2(EvtSpinType::TENSOR,EvtSpinType::SCALAR,
		   EvtSpinType::SCALAR,Hd2);


  EvtVector4R theKstarP4boost(theKstarP4.get(0),-theKstarP4.get(1),-theKstarP4.get(2),-theKstarP4.get(3));

  EvtScalarParticle theKaon;
  theKaon.init(EvtPDL::getId(std::string("K+")),boostTo(kaon->getP4(),theKstarP4boost));

  EvtScalarParticle thePion;
  thePion.init(EvtPDL::getId(std::string("pi+")),boostTo(pion->getP4(),theKstarP4boost));

  theKaon.addDaug(&theKstar);
  thePion.addDaug(&theKstar);

  // Calculate the propagator

  double m = theKstarP4.mass();
  EvtTwoBodyVertex v(0.5,0.14,1.4,2);
  EvtTwoBodyKine v1(0.5,0.14,m);
  EvtPropBreitWignerRel prop(1.4,0.2);

  // Mass-dependent width correction and amplitude calculation

  double width = prop.g0() * v.widthFactor(v1);
  prop.set_g0(width);
  EvtComplex bwamp = prop.evaluate(m);


  EvtAmp amp2;

  d2.evalAmp(&theKstar,amp2);

  vertex(0,bwamp*(amp1._amp[0]*amp2._amp[0]+
	   amp1._amp[1]*amp2._amp[1]+
	   amp1._amp[2]*amp2._amp[2]+
	   amp1._amp[3]*amp2._amp[3]+
           amp1._amp[4]*amp2._amp[4]));

  vertex(1,bwamp*(amp1._amp[5]*amp2._amp[0]+
	   amp1._amp[6]*amp2._amp[1]+
	   amp1._amp[7]*amp2._amp[2]+
	   amp1._amp[8]*amp2._amp[3]+
           amp1._amp[9]*amp2._amp[4]));

*/

    return;
}
