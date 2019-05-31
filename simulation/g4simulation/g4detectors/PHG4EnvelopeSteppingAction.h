#ifndef G4DETECTORS_PHG4ENVELOPESTEPPINGACTION_H
#define G4DETECTORS_PHG4ENVELOPESTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class PHCompositeNode;
class PHG4EnvelopeDetector;
class PHG4Hit;
class PHG4HitContainer;

class PHG4EnvelopeSteppingAction: public PHG4SteppingAction
{
	public:
		//Constructor
		PHG4EnvelopeSteppingAction( PHG4EnvelopeDetector* );
	
		//Destructor
		virtual ~PHG4EnvelopeSteppingAction()
		{}
	
		//Stepping Action
		virtual bool UserSteppingAction( const G4Step*, bool);
	
		//reimplemented from base class
		virtual void SetInterfacePointers( PHCompositeNode* );
	
	private:
		
		//pointer to the detector
		PHG4EnvelopeDetector* detector_;
	
		//pointer to hit container
		PHG4HitContainer * hits_;
		PHG4Hit *hit;	
};

#endif //PHG4EnvelopeSteppingAction_h
