#include <iostream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <cmath>

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
	if ((*dist1).dR < (*dist2).dR) { return -1;
	} if ((*dist1).dR > (*dist2).dR) {return 1; }
	return 0;
}

double sqr(double x){
	return x*x;
}

std::vector<MyPart> readEvent(int (*ptrK)[23000], double (*ptrP)[23000], int pyjetsN, int lstart)
{
	//std::cout<<"Will start reading the event record with "<<pyjetsN<<" particles now.\n";
	std::vector<MyPart> parts;
	for (int i = lstart; i < pyjetsN; ++i) {
		if (ptrK[0][i] < 6 || ptrK[0][i] == 17) {
			//std::cout<<i<<": "<<ptrK[0][i]<<" - "<<ptrK[1][i]<<" -> ("<<ptrP[0][i]<<" , "<<ptrP[1][i]<<" , "<<ptrP[2][i]<<" ; "<<ptrP[3][i]<<") m = "<<ptrP[4][i]<<std::endl;
			if (ptrP[3][i] - fabs(ptrP[2][i]) > 0.) {
				MyPart part{};
				double pT;
				double mass2;
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
				std::cout<<i<<": "<<ptrK[0][i]<<" - "<<ptrK[1][i]<<" -> ("<<ptrP[0][i]<<" , "<<ptrP[1][i]<<" , "<<ptrP[2][i]<<" ; "<<ptrP[3][i]<<") m = "<<ptrP[4][i]<<std::endl;
				//_subtrevent.push_back(p);
			}
		}
	}
	//std::cout<<"Done with reading the event record, have found "<<parts.size()<<" final state particles.\n";
	return parts;
}

std::vector<My4Mom> readThermals(double (*ptrscatcens)[23000], int nscatcen) {
	//std::cout<<"Will start reading "<<nscatcen<<" thermal momenta now.\n";
	std::vector<My4Mom> thermals;
	for (int i = 0; i < nscatcen; ++i) {
		if (ptrscatcens[3][i] - fabs(ptrscatcens[2][i]) > 0.) {
			My4Mom mom{};
			double pT;
			double mass2;
			pT = sqrt(sqr(ptrscatcens[0][i]) + sqr(ptrscatcens[1][i]));
			mass2 = sqr(ptrscatcens[3][i]) - sqr(pT) - sqr(ptrscatcens[2][i]);
			mom.pt = pT;
			mom.mdelta = sqrt(mass2 + sqr(pT)) - pT;
			mom.phi = atan2(ptrscatcens[1][i],ptrscatcens[0][i]);
			mom.y = 0.5*log((ptrscatcens[3][i]+ptrscatcens[2][i])/(ptrscatcens[3][i]-ptrscatcens[2][i]));
			thermals.push_back(mom);
		}
		else {
			std::cout<<i<<": "<<ptrscatcens[0][i]<<" - "<<ptrscatcens[1][i]<<" -> ("<<ptrscatcens[0][i]<<" , "<<ptrscatcens[1][i]<<" , "<<ptrscatcens[2][i]<<" ; "<<ptrscatcens[3][i]<<") m = "<<ptrscatcens[4][i]<<std::endl;
			//_subtrevent.push_back(p);
		}
	}
	//std::cout<<"Done with reading the thermal momenta.\n";	
	return thermals;
}

std::vector<MyDist> buildPairs(std::vector<MyPart> * myparticles, std::vector<My4Mom> * mythermals) {
    //create list with all particle-ghosts distances and sort it
	//std::cout<<"Will start creating the pairs now.\n";
	  std::vector<MyDist> dists;
	  for (size_t i=0; i < myparticles->size(); ++i) {
		  for (size_t j=0; j < mythermals->size(); ++j) {
			  //std::cout<<"i, j = "<<i<<" "<<j<<std::endl;
			  MyDist dist{};
			  double Deltaphi(abs((*myparticles)[i].phi-(*mythermals)[j].phi));
			  if (Deltaphi > M_PI) { Deltaphi = 2.*M_PI - Deltaphi;
}
			  dist.dR = sqrt(sqr(Deltaphi) + sqr((*myparticles)[i].y-(*mythermals)[j].y));
			  dist.ipart = i;
			  dist.itherm = j;
			  //std::cout<<"distance: "<<dist.dR<<std::endl;
			  dists.push_back(dist);
		  }
	  }
	  //dists.sort(MyDistComp);
	  qsort(dists.data(), dists.size(), sizeof(MyDist), MyDistComp);
	  //std::cout<<"Done with creating the pairs.\n";
	  return dists;
}
    
void doSubtraction(std::vector<MyDist> * dists, std::vector<MyPart> * myparticles, std::vector<My4Mom> * mythermals, double dRmax) {
	//std::cout<<"Will do the subtraction now with max distance "<<dRmax<<std::endl;
		  //go through all particle-ghost pairs and re-distribute momentum and mass
	  for (auto & dist : *dists) {
		  //std::cout<<"dealing with dist "<<liter->dR<<std::endl;
		  if (dist.dR > dRmax) { break;
}
		  size_t pnum = dist.ipart;
		  size_t tnum = dist.itherm;
		  double ptp = (*myparticles)[pnum].pt;
		  double ptt = (*mythermals)[tnum].pt;
		  //std::cout<<"pts: "<<ptp<<" vs "<<ptg<<std::endl;
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
		  //std::cout<<"masses: "<<mdp<<" vs "<<mdg<<std::endl;
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
    
void collectResult(std::vector<MyPart> * myparticles, std::vector<My4Mom> * mythermals, int (*ptrK)[23000], double (*ptrP)[23000], int * pyjetsN, int lstart) {
		  //collect resulting 4-momenta to get subtracted event
	int N(lstart);
	// lstart is an offset into PYJETS, NOT into the compacted myparticles vector,
	// which readEvent() already started at lstart.  Indexing it here dropped the
	// first lstart final-state particles from every subtracted event.
	for (auto & myparticle : *myparticles) {
		if (myparticle.pt > 0.) {
			N++;
			if (N == 23000) {
				*pyjetsN = 0;
				return;
			}
			ptrK[0][N-1] = 1;
			ptrK[1][N-1] = myparticle.id;
			ptrK[2][N-1] = 0;
			ptrK[3][N-1] = 0;
			ptrK[4][N-1] = 0;
			ptrP[0][N-1] = myparticle.pt*cos(myparticle.phi);
			ptrP[1][N-1] = myparticle.pt*sin(myparticle.phi);
			ptrP[2][N-1] = (myparticle.pt+myparticle.mdelta)*sinh(myparticle.y);
			ptrP[3][N-1] = (myparticle.pt+myparticle.mdelta)*cosh(myparticle.y);
			ptrP[4][N-1] = sqrt(sqr(myparticle.pt+myparticle.mdelta) - sqr(myparticle.pt));
		  }
	  }
	  for (auto & mythermal : *mythermals) {
		if (mythermal.pt > 0.) {
			int thermid;
			double R(double(rand())/double(RAND_MAX));
			if (R < 1./3.) { thermid = 211;
			} else if (R < 2./3.) { thermid = -211;
			} else { thermid = 111;
}
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
			ptrP[0][N-1] = mythermal.pt*cos(mythermal.phi);
			ptrP[1][N-1] = mythermal.pt*sin(mythermal.phi);
			ptrP[2][N-1] = (mythermal.pt+mythermal.mdelta)*sinh(mythermal.y);
			ptrP[3][N-1] = (mythermal.pt+mythermal.mdelta)*cosh(mythermal.y);
			ptrP[4][N-1] = sqrt(sqr(mythermal.pt+mythermal.mdelta) - sqr(mythermal.pt));
		}
	  }
	  *pyjetsN = N;
        return;
}


extern "C" {
    void subtract_thmom_(int pyjetsK[5][23000], double pyjetsP[5][23000], int* pyjetsN, double scatcens[5][23000], const int* nscatcen, const double * Rmax, char *proc);
}

void subtract_thmom_(int pyjetsK[5][23000], double pyjetsP[5][23000], int* pyjetsN, double scatcens[5][23000], const int* nscatcen, const double * Rmax, char *proc) 
{
	// proc is a Fortran character*4: NOT NUL-terminated, so compare 4 chars only.
	// strcmp()==0 means EQUAL -- the original tests were inverted, giving lstart=2
	// for every process.  Mapping below matches putparticlesonshell() in
	// jewel-2.6.0.f, which is the authors' own statement of the intent.
	int lstart(0);
	if (!strncmp(proc,"PPYJ",4) || !strncmp(proc,"PPYQ",4)
	    || !strncmp(proc,"PPYG",4)) { lstart = 1;
}
	if (!strncmp(proc,"PPZJ",4) || !strncmp(proc,"PPWJ",4)) { lstart = 2;
}
	//std::cout<<"Lstart = "<<lstart<<std::endl;
	//std::cout<<"Will do the subtraction now.\n";
	std::vector<MyPart> myparticles = readEvent(pyjetsK, pyjetsP, *pyjetsN, lstart);
	std::vector<My4Mom> mythermals = readThermals(scatcens, *nscatcen);
	std::vector<MyDist> dists = buildPairs(&myparticles, &mythermals);
	doSubtraction(&dists, &myparticles, &mythermals, *Rmax);
	collectResult(&myparticles, &mythermals, pyjetsK, pyjetsP, pyjetsN, lstart);
	//std::cout<<"Done with the subtraction, resulted in "<<*pyjetsN<<" particles.\n";
	return;
}

