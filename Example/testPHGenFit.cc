#include <vector>

#include <TVector3.h>
#include <TMatrixDSym.h>

#include <GenFit/AbsTrackRep.h>
#include <GenFit/RKTrackRep.h>

#include <phgenfit/Fitter.h>
#include <phgenfit/Track.h>
#include <phgenfit/Measurement.h>
#include <phgenfit/PlanarMeasurement.h>


int main(int argc, char**argv)
{
	PHGenFit::Fitter* fitter = new PHGenFit::Fitter("sPHENIX_TGeo.root","sPHENIX.2d.root",1.4/1.5);

	int pid_assumption = -13; //mu+
	genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pid_assumption);

	TMatrixDSym cov(6);
	PHGenFit::Track* track = new PHGenFit::Track(rep, TVector3(0,0,0), TVector3(0,0,0), cov);

	std::vector<PHGenFit::Measurement*> measurements;
	for(int imeas = 0; imeas < 3; imeas++)
	{
		PHGenFit::Measurement* meas = PHGenFit::PlanarMeasurement(TVector3(0,0,0),TVector3(1,0,0),TVector3(0,1,0),1, 1);
		measurements.push_back(meas);
	}


	track->addMeasurements(measurements);

	fitter->processTrack(track);


	return 0;
}
