/*!
 *** \\ class G4SnglHit
 *** \\ author Dr. Liang Xue from Georgia State University
 */

#ifndef G4SNGLHIT_H__
#define G4SNGLHIT_H__

#include "phool/PHObject.h"
#include <iostream>

class G4SnglHit: public PHObject
{

	public:
		//! constructor
		G4SnglHit();

		//! destructor
		virtual ~G4SnglHit();

		void Reset();

		int get_detid() const { return detid; }
		int get_layer() const { return layer; }
		int get_scintid() const { return scintid; }
		int get_trackid() const { return trackid; }
		float get_x0() const { return x0; }
		float get_y0() const { return y0; }
		float get_z0() const { return z0; }
		float get_x1() const { return x1; }
		float get_y1() const { return y1; }
		float get_z1() const { return z1; }
		float get_edep() const { return edep; }

		void set_detid( int i ) { detid = i; }
		void set_layer( int i ) { layer = i; }
		void set_scintid( int i ) { scintid = i; }
		void set_trackid( int i ) { trackid = i; }
		void set_x0( float f ) { x0 = f; }
		void set_y0( float f ) { y0 = f; }
		void set_z0( float f ) { z0 = f; }
		void set_x1( float f ) { x1 = f; }
		void set_y1( float f ) { y1 = f; }
		void set_z1( float f ) { z1 = f; }
		void set_edep( float f ) { edep = f; }

		virtual void print() const;

	public:
		int detid;
		int layer;
		int scintid;
		int trackid;
		float x0;
		float y0;
		float z0;
		float x1;
		float y1;
		float z1;
		float edep;

		ClassDef(G4SnglHit,1)
};

#endif 
