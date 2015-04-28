/*!
 *** \\ class G4SnglCluster
 *** \\ author Dr. Liang Xue from Georgia State University
 */

#ifndef G4SNGLCLUSTER_H__
#define G4SNGLCLUSTER_H__

#include "phool/PHObject.h"
#include <iostream>

class G4SnglCluster: public PHObject
{

	public:
		//! constructor
		G4SnglCluster();

		//! destructor
		virtual ~G4SnglCluster();

		void Reset();

		int get_detid() const { return detid; }
		int get_ntowers() const { return ntowers; }
		float get_eta() const { return eta; }
		float get_phi() const { return phi; }
		float get_edep() const { return edep; }

		void set_detid( int i ) { detid = i; }
		void set_ntowers( int i ) { ntowers = i; }
		void set_eta( float f ) { eta = f; }
		void set_phi( float f ) { phi = f; }
		void set_edep( float f ) { edep = f; }

		virtual void print() const;

	public:
		int detid;
		int ntowers;
		float eta;
		float phi;
		float edep;

		ClassDef(G4SnglCluster,1)
};

#endif 
