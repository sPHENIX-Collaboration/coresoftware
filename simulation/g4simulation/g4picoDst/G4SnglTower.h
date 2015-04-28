/*!
 *** \\ class G4SnglTower
 *** \\ author Dr. Liang Xue from Georgia State University
 */

#ifndef G4SNGLTOWER_H__
#define G4SNGLTOWER_H__

#include "phool/PHObject.h"
#include <iostream>

class G4SnglTower: public PHObject
{

	public:
		//! constructor
		G4SnglTower();

		//! destructor
		virtual ~G4SnglTower();

		void Reset();

		int get_detid() const { return detid; }
		float get_ieta() const { return ieta; }
		float get_iphi() const { return iphi; }
		float get_eta() const { return eta; }
		float get_phi() const { return phi; }
		float get_edep() const { return edep; }

		void set_detid( int i ) { detid = i; }
		void set_ieta( float f ) { ieta = f; }
		void set_iphi( float f ) { iphi = f; }
		void set_eta( float f ) { eta = f; }
		void set_phi( float f ) { phi = f; }
		void set_edep( float f ) { edep = f; }

		virtual void print() const;

	public:
		int detid;
		float ieta;
		float iphi;
		float eta;
		float phi;
		float edep;

		ClassDef(G4SnglTower,1)
};

#endif 
