/*!
 *** \\ class G4SnglHTCContainer, container for G4 Hits, Towers, and Clusters
 *** \\ author Dr. Liang Xue from Georgia State University
 */

#ifndef G4SNGLHTCCONTAINER_H__
#define G4SNGLHTCCONTAINER_H__

#include <iostream>
#include "phool/PHObject.h"
#include "TClonesArray.h"

class G4SnglHit;
class G4SnglTower;
class G4SnglCluster;

class G4SnglHTCContainer: public PHObject
{
	public:
		//! constructor
		G4SnglHTCContainer();

		//! destructor
		virtual ~G4SnglHTCContainer();

		virtual unsigned int get_nG4SnglHits() const {return nG4SnglHits;}
		virtual unsigned int get_nG4SnglTowers() const {return nG4SnglTowers;}
		virtual unsigned int get_nG4SnglClusters() const {return nG4SnglClusters;}

		void Reset();

		G4SnglHit* AddG4SnglHit(const G4SnglHit &snglhit);
		G4SnglTower* AddG4SnglTower(const G4SnglTower &sngltwr);
		G4SnglCluster* AddG4SnglCluster(const G4SnglCluster &snglclr);

		TClonesArray*   GetG4SnglHitContainer() {return G4SnglHits;}
		TClonesArray*   GetG4SnglTowerContainer() {return G4SnglTowers;}
		TClonesArray*   GetG4SnglClusterContainer() {return G4SnglClusters;}

		int get_PID() const { return PID; }
		float get_Energy() const { return Energy; }
		float get_Theta() const { return Theta; }
		float get_Phi() const { return Phi; }
		float get_Px() const { return Px; }
		float get_Py() const { return Py; }
		float get_Pz() const { return Pz; }

		void set_PID(int i) { PID = i; }
		void set_Energy(float f) { Energy = f; }
		void set_Phi(float f) { Phi = f; }
		void set_Theta(float f) { Theta = f; }
		void set_Px(float f) { Px = f; }
		void set_Py(float f) { Py = f; }
		void set_Pz(float f) { Pz = f; }

	public:

		unsigned int nG4SnglHits;
		TClonesArray *G4SnglHits;
		unsigned int nG4SnglTowers;
		TClonesArray *G4SnglTowers;
		unsigned int nG4SnglClusters;
		TClonesArray *G4SnglClusters;

		int PID;
		float Energy;
		float Phi;
		float Theta;
		float Px;
		float Py;
		float Pz;


		ClassDef(G4SnglHTCContainer,1)
};

#endif 
