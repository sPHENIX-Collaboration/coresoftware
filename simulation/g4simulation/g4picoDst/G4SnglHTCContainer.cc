#include "G4SnglHit.h"
#include "G4SnglTower.h"
#include "G4SnglCluster.h"
#include "G4SnglHTCContainer.h"

#include <phool/phool.h>

#include <iostream>

using namespace std;

static const unsigned int NMAX = 100000;

G4SnglHTCContainer::G4SnglHTCContainer()
{
	nG4SnglHits = 0;
        G4SnglHits = new TClonesArray("G4SnglHit",NMAX);
	nG4SnglTowers = 0;
	G4SnglTowers = new TClonesArray("G4SnglTower",NMAX);
	nG4SnglClusters = 0;
	G4SnglClusters = new TClonesArray("G4SnglCluster",NMAX);
	Reset();
}

G4SnglHTCContainer::~G4SnglHTCContainer()
{
	G4SnglHits->Clear();
	G4SnglTowers->Clear();
	G4SnglClusters->Clear();
	delete G4SnglHits;
	delete G4SnglTowers;
	delete G4SnglClusters;
	return;
}


G4SnglHit*
G4SnglHTCContainer::AddG4SnglHit(const G4SnglHit &snglhit)
{
	// The code below needs a bit of explanation.
	// Before this code is called, TClonesArray pointed to by
	// G4SnglHits has some size.
	// When ExpandCreate() is called, a new G4SnglHit is 
	// created on the TClonesArray at the last positon. If there
	// is not sufficient space in the TClonesArray, the TClonesArray
	// is expaneded to store the new G4SnglHit object. Since
	// the creation is done by default constructor of G4SnglHit, the
	// data of "snglhit" have to be copied to the newly created
	// object.

	if( !G4SnglHits ) return NULL;

	// this is a bit ugly but GetLast returns the index-1.
	int nnew = G4SnglHits->GetLast()+2;
	G4SnglHits->ExpandCreate(nnew);

	// ExpandCreate calls default constrcutor. So, the contents of snglhit
	// should be copied to newsnglmu

	G4SnglHit *newsnglhit = (G4SnglHit*) ( G4SnglHits->UncheckedAt(G4SnglHits->GetLast()) );
	*newsnglhit = snglhit;
	nG4SnglHits++;

	return newsnglhit;
}
	

G4SnglTower*
G4SnglHTCContainer::AddG4SnglTower(const G4SnglTower &sngltwr)
{
	// The code below needs a bit of explanation.
	// Before this code is called, TClonesArray pointed to by
	// G4SnglTowers has some size.
	// When ExpandCreate() is called, a new G4SnglTower is 
	// created on the TClonesArray at the last positon. If there
	// is not sufficient space in the TClonesArray, the TClonesArray
	// is expaneded to store the new G4SnglTower object. Since
	// the creation is done by default constructor of G4SnglTower, the
	// data of "sngltwr" have to be copied to the newly created
	// object.

	if( !G4SnglTowers ) return NULL;

	// this is a bit ugly but GetLast returns the index-1.
	int nnew = G4SnglTowers->GetLast()+2;
	G4SnglTowers->ExpandCreate(nnew);
	
	// ExpandCreate calls default constrcutor. So, the contents of sngltwr
	// should be copied to newsnglmu

	G4SnglTower *newsngltwr = (G4SnglTower*) ( G4SnglTowers->UncheckedAt(G4SnglTowers->GetLast()) );
	*newsngltwr = sngltwr;
	nG4SnglTowers++;

	return newsngltwr;
}

G4SnglCluster*
G4SnglHTCContainer::AddG4SnglCluster(const G4SnglCluster &snglclr)
{
	// The code below needs a bit of explanation.
	// Before this code is called, TClonesArray pointed to by
	// G4SnglClusters has some size.
	// When ExpandCreate() is called, a new G4SnglCluster is 
	// created on the TClonesArray at the last positon. If there
	// is not sufficient space in the TClonesArray, the TClonesArray
	// is expaneded to store the new G4SnglCluster object. Since
	// the creation is done by default constructor of G4SnglCluster, the
	// data of "snglclr" have to be copied to the newly created
	// object.

	if( !G4SnglClusters ) return NULL;

	// this is a bit ugly but GetLast returns the index-1.
	int nnew = G4SnglClusters->GetLast()+2;
	G4SnglClusters->ExpandCreate(nnew);
	
	// ExpandCreate calls default constrcutor. So, the contents of snglclr
	// should be copied to newsnglmu

	G4SnglCluster *newsnglclr = (G4SnglCluster*) ( G4SnglClusters->UncheckedAt(G4SnglClusters->GetLast()) );
	*newsnglclr = snglclr;
	nG4SnglClusters++;

	return newsnglclr;
}

void
G4SnglHTCContainer::Reset()
{
	PID = -9999;
	Energy = -9999.;
	Theta = -9999.;
	Phi = -9999.;
	Px = -9999.;
	Py = -9999.;
	Pz = -9999.;
	
	G4SnglHits->Clear();
	G4SnglTowers->Clear();
	G4SnglClusters->Clear();
	if( nG4SnglHits > NMAX) { G4SnglHits->Expand(NMAX); }
	if( nG4SnglTowers > NMAX) { G4SnglTowers->Expand(NMAX); }
	if( nG4SnglClusters > NMAX) { G4SnglClusters->Expand(NMAX); }
	nG4SnglHits = 0;
	nG4SnglTowers = 0;
	nG4SnglClusters = 0;
	return;
}
