#include <iostream>
#include <string.h>
#include <cstdlib>
#include <vector>
#include <math.h>
using namespace std;

struct My4Mom {	
	double pt;
	double mdelta;
	double phi;
	double y;
} ;

struct MyPart {
	double pt;
	double mdelta;
	double phi;
	double y;
	int    id;
} ;

struct MyDist {
	double dR;
	size_t ipart;
	size_t itherm;
} ;
    
int MyDistComp (const void * a, const void * b) {
	MyDist *dist1 = (MyDist *) a;
	MyDist *dist2 = (MyDist *) b;
	if ((*dist1).dR < (*dist2).dR) return -1;
	else if ((*dist1).dR > (*dist2).dR) return 1; 
	return 0;
}

double sqr(double x){
	return x*x;
}

vector<MyPart> readEvent(int (*ptrK)[23000], double (*ptrP)[23000], int pyjetsN, int lstart)
{
	//cout<<"Will start reading the event record with "<<pyjetsN<<" particles now.\n";
	vector<MyPart> parts;
	for (int i = lstart; i < pyjetsN; ++i) {
		if (ptrK[0][i] < 6 || ptrK[0][i] == 17) {
			//cout<<i<<": "<<ptrK[0][i]<<" - "<<ptrK[1][i]<<" -> ("<<ptrP[0][i]<<" , "<<ptrP[1][i]<<" , "<<ptrP[2][i]<<" ; "<<ptrP[3][i]<<") m = "<<ptrP[4][i]<<endl;
			if (ptrP[3][i] - fabs(ptrP[2][i]) > 0.) {
				MyPart part;
				double pT,mass2;
				pT = sqrt(sqr(ptrP[0][i]) + sqr(ptrP[1][i]));
				mass2 = sqr(ptrP[3][i]) - sqr(pT) - sqr(ptrP[2][i]);
				part.id = ptrK[1][i];
				part.pt = pT;
				part.mdelta = sqrt(mass2 + sqr(pT)) - pT;
				part.phi = atan2(ptrP[1][i],ptrP[0][i]);
				part.y = 0.5*log((ptrP[3][i]+ptrP[2][i])/(ptrP[3][i]-ptrP[2][i]));
				parts.push_back(part);
			}
			else {
				cout<<i<<": "<<ptrK[0][i]<<" - "<<ptrK[1][i]<<" -> ("<<ptrP[0][i]<<" , "<<ptrP[1][i]<<" , "<<ptrP[2][i]<<" ; "<<ptrP[3][i]<<") m = "<<ptrP[4][i]<<endl;
				//_subtrevent.push_back(p);
			}
		}
	}
	//cout<<"Done with reading the event record, have found "<<parts.size()<<" final state particles.\n";
	return parts;
}

vector<My4Mom> readThermals(double (*ptrscatcens)[23000], int nscatcen) {
	//cout<<"Will start reading "<<nscatcen<<" thermal momenta now.\n";
	vector<My4Mom> thermals;
	for (int i = 0; i < nscatcen; ++i) {
		if (ptrscatcens[3][i] - fabs(ptrscatcens[2][i]) > 0.) {
			My4Mom mom;
			double pT,mass2;
			pT = sqrt(sqr(ptrscatcens[0][i]) + sqr(ptrscatcens[1][i]));
			mass2 = sqr(ptrscatcens[3][i]) - sqr(pT) - sqr(ptrscatcens[2][i]);
			mom.pt = pT;
			mom.mdelta = sqrt(mass2 + sqr(pT)) - pT;
			mom.phi = atan2(ptrscatcens[1][i],ptrscatcens[0][i]);
			mom.y = 0.5*log((ptrscatcens[3][i]+ptrscatcens[2][i])/(ptrscatcens[3][i]-ptrscatcens[2][i]));
			thermals.push_back(mom);
		}
		else {
			cout<<i<<": "<<ptrscatcens[0][i]<<" - "<<ptrscatcens[1][i]<<" -> ("<<ptrscatcens[0][i]<<" , "<<ptrscatcens[1][i]<<" , "<<ptrscatcens[2][i]<<" ; "<<ptrscatcens[3][i]<<") m = "<<ptrscatcens[4][i]<<endl;
			//_subtrevent.push_back(p);
		}
	}
	//cout<<"Done with reading the thermal momenta.\n";	
	return thermals;
}

vector<MyDist> buildPairs(vector<MyPart> * myparticles, vector<My4Mom> * mythermals) {
    //create list with all particle-ghosts distances and sort it
	//cout<<"Will start creating the pairs now.\n";
	  vector<MyDist> dists;
	  for (size_t i=0; i < myparticles->size(); ++i) {
		  for (size_t j=0; j < mythermals->size(); ++j) {
			  //cout<<"i, j = "<<i<<" "<<j<<endl;
			  MyDist dist;
			  double Deltaphi(abs((*myparticles)[i].phi-(*mythermals)[j].phi));
			  if (Deltaphi > M_PI) Deltaphi = 2.*M_PI - Deltaphi;
			  dist.dR = sqrt(sqr(Deltaphi) + sqr((*myparticles)[i].y-(*mythermals)[j].y));
			  dist.ipart = i;
			  dist.itherm = j;
			  //cout<<"distance: "<<dist.dR<<endl;
			  dists.push_back(dist);
		  }
	  }
	  //dists.sort(MyDistComp);
	  qsort(&dists[0], dists.size(), sizeof(MyDist), MyDistComp);
	  //cout<<"Done with creating the pairs.\n";
	  return dists;
}
    
void doSubtraction(vector<MyDist> * dists, vector<MyPart> * myparticles, vector<My4Mom> * mythermals, double dRmax) {
	//cout<<"Will do the subtraction now with max distance "<<dRmax<<endl;
		  //go through all particle-ghost pairs and re-distribute momentum and mass
	  for (vector<MyDist>::iterator liter=dists->begin(); liter != dists->end(); ++liter) {
		  //cout<<"dealing with dist "<<liter->dR<<endl;
		  if (liter->dR > dRmax) break;
		  size_t pnum = liter->ipart;
		  size_t tnum = liter->itherm;
		  double ptp = (*myparticles)[pnum].pt;
		  double ptt = (*mythermals)[tnum].pt;
		  //cout<<"pts: "<<ptp<<" vs "<<ptg<<endl;
		  if (ptp > ptt) {
			  (*myparticles)[pnum].pt -= ptt;
			  (*mythermals)[tnum].pt = 0.;
		  }
		  else {
			  (*mythermals)[tnum].pt -= ptp;
			  (*myparticles)[pnum].pt = 0.;
		  }
		  double mdp = (*myparticles)[pnum].mdelta;
		  double mdt = (*mythermals)[tnum].mdelta;
		  //cout<<"masses: "<<mdp<<" vs "<<mdg<<endl;
		  if (mdp > mdt) {
			  (*myparticles)[pnum].mdelta -= mdt;
			  (*mythermals)[tnum].mdelta = 0.;
		  }
		  else {
			  (*mythermals)[tnum].mdelta -= mdp;
			  (*myparticles)[pnum].mdelta = 0.;
		  }
	  }
	  return;
}
    
void collectResult(vector<MyPart> * myparticles, vector<My4Mom> * mythermals, int (*ptrK)[23000], double (*ptrP)[23000], int * pyjetsN, int lstart) {
		  //collect resulting 4-momenta to get subtracted event
	int N(lstart);
	// lstart is an offset into PYJETS, NOT into the compacted myparticles vector,
	// which readEvent() already started at lstart.  Indexing it here dropped the
	// first lstart final-state particles from every subtracted event.
	for (size_t i=0; i < myparticles->size(); ++i) {
		if ((*myparticles)[i].pt > 0.) {
			N++;
			if (N == 23000) {
				*pyjetsN = 0;
				return;
			}
			ptrK[0][N-1] = 1;
			ptrK[1][N-1] = (*myparticles)[i].id;
			ptrK[2][N-1] = 0;
			ptrK[3][N-1] = 0;
			ptrK[4][N-1] = 0;
			ptrP[0][N-1] = (*myparticles)[i].pt*cos((*myparticles)[i].phi);
			ptrP[1][N-1] = (*myparticles)[i].pt*sin((*myparticles)[i].phi);
			ptrP[2][N-1] = ((*myparticles)[i].pt+(*myparticles)[i].mdelta)*sinh((*myparticles)[i].y);
			ptrP[3][N-1] = ((*myparticles)[i].pt+(*myparticles)[i].mdelta)*cosh((*myparticles)[i].y);
			ptrP[4][N-1] = sqrt(sqr((*myparticles)[i].pt+(*myparticles)[i].mdelta) - sqr((*myparticles)[i].pt));
		  }
	  }
	  for (size_t i=0; i < mythermals->size(); ++i) {
		if ((*mythermals)[i].pt > 0.) {
			int thermid;
			double R(double(rand())/double(RAND_MAX));
			if (R < 1./3.) thermid = 211;
			else if (R < 2./3.) thermid = -211;
			else thermid = 111;
			N++;
			if (N == 23000) {
				*pyjetsN = 0;
				return;
			}
			ptrK[0][N-1] = 1;
			ptrK[1][N-1] = thermid;
			ptrK[2][N-1] = 0;
			ptrK[3][N-1] = 0;
			ptrK[4][N-1] = 0;
			ptrP[0][N-1] = (*mythermals)[i].pt*cos((*mythermals)[i].phi);
			ptrP[1][N-1] = (*mythermals)[i].pt*sin((*mythermals)[i].phi);
			ptrP[2][N-1] = ((*mythermals)[i].pt+(*mythermals)[i].mdelta)*sinh((*mythermals)[i].y);
			ptrP[3][N-1] = ((*mythermals)[i].pt+(*mythermals)[i].mdelta)*cosh((*mythermals)[i].y);
			ptrP[4][N-1] = sqrt(sqr((*mythermals)[i].pt+(*mythermals)[i].mdelta) - sqr((*mythermals)[i].pt));
		}
	  }
	  *pyjetsN = N;
        return;
}


extern "C" {
    void subtract_thmom_(int pyjetsK[5][23000], double pyjetsP[5][23000], int* pyjetsN, double scatcens[5][23000], int* nscatcen, double * Rmax, char *proc);
}

void subtract_thmom_(int pyjetsK[5][23000], double pyjetsP[5][23000], int* pyjetsN, double scatcens[5][23000], int* nscatcen, double * Rmax, char *proc) 
{
	// proc is a Fortran character*4: NOT NUL-terminated, so compare 4 chars only.
	// strcmp()==0 means EQUAL -- the original tests were inverted, giving lstart=2
	// for every process.  Mapping below matches putparticlesonshell() in
	// jewel-2.6.0.f, which is the authors' own statement of the intent.
	int lstart(0);
	if (!strncmp(proc,"PPYJ",4) || !strncmp(proc,"PPYQ",4)
	    || !strncmp(proc,"PPYG",4)) lstart = 1;
	if (!strncmp(proc,"PPZJ",4) || !strncmp(proc,"PPWJ",4)) lstart = 2;
	//cout<<"Lstart = "<<lstart<<endl;
	//cout<<"Will do the subtraction now.\n";
	vector<MyPart> myparticles = readEvent(pyjetsK, pyjetsP, *pyjetsN, lstart);
	vector<My4Mom> mythermals = readThermals(scatcens, *nscatcen);
	vector<MyDist> dists = buildPairs(&myparticles, &mythermals);
	doSubtraction(&dists, &myparticles, &mythermals, *Rmax);
	collectResult(&myparticles, &mythermals, pyjetsK, pyjetsP, pyjetsN, lstart);
	//cout<<"Done with the subtraction, resulted in "<<*pyjetsN<<" particles.\n";
	return;
}

