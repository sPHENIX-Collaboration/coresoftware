#ifndef G4DETECTORS_PHG4ENVELOPESUBSYSTEM_H
#define G4DETECTORS_PHG4ENVELOPESUBSYSTEM_H

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4String.hh>

class PHG4Detector;
class PHG4EnvelopeDetector;
class PHG4EnvelopeSteppingAction;
class PHG4EventAction;
class PHG4SteppingAction;
class PHCompositeNode;

class PHG4EnvelopeSubsystem: public PHG4Subsystem
{
	public:
		//Constructor
		PHG4EnvelopeSubsystem ( const std::string &name = "ENVELOPE_DEFAULT", const int layer = 0 );
		
		//Destructor
		virtual ~PHG4EnvelopeSubsystem ( void )
		{}
			
		/*
			Creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
			Creates the stepping action and place it on the node tree, under "ACTIONS" node
			Creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
		*/
	
		int Init(PHCompositeNode *);
	
		//Event Processing
		int process_event(PHCompositeNode *);
	
		//Accessors (reimplemented)
		virtual PHG4Detector* GetDetector( void ) const;
		virtual PHG4SteppingAction* GetSteppingAction( void ) const;
	
	private:
		//Pointer to Geant4 implementation of detector
		PHG4EnvelopeDetector* detector_;
	
		//Stepping Action
		PHG4EnvelopeSteppingAction* steppingAction_;
		PHG4EventAction *eventAction_;
			
		G4String material;
		int active;
	
		std::string detector_type;
			
};

#endif
