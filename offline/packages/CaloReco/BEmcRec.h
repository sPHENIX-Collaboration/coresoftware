#ifndef BEMCREC_H
#define BEMCREC_H

// Name: BEmcRec.h
// Author: A. Bazilevsky, Apr 2012
// Modified from EmcSectorRec.h and EmcScSectorRec.h

//#include "PHMatrix.h"
//#include "PHVector.h"
#include <vector>


class EmcCluster;
class EmcModule;

typedef struct SecGeom{

  short nx;          // Number of cells in X dir
  short ny;          // Number of cells in Y dir
  float Tower_xSize; // Tower size in X dir
  float Tower_ySize; // Tower size in Y dir

} SecGeom;

// ///////////////////////////////////////////////////////////////////////////

/** ABC of a clusterizer for one EMCAL sector. 

@ingroup clustering

 */

class BEmcRec
{
  
 public:

  BEmcRec();
  ~BEmcRec();

  void SetVertex(float *vv){  fVx = vv[0]; fVy = vv[1]; fVz = vv[2]; }
  void SetGeometry(int nx, int ny, float txsz, float tysz);
  //  void SetGeometry(SecGeom const &geom, PHMatrix * rm, PHVector * tr );
  void SetConf(int nx, int ny) { SetGeometry(nx, ny, 1., 1.); }

  void SetPlaneGeometry() { bCYL=false; }
  void SetCylindricalGeometry() { bCYL=true; }
  bool isCylindrical() const { return bCYL; }

  int GetNx() const { return fNx; }
  int GetNy() const { return fNy; }
  float GetModSizex() const { return fModSizex; }
  float GetModSizey() const { return fModSizey; }
  float GetVx() const { return fVx; }
  float GetVy() const { return fVy; }
  float GetVz() const { return fVz; }
  void SetPeakThreshold(float Thresh) { fgMinPeakEnergy = Thresh; }
  float GetPeakThreshold() { return fgMinPeakEnergy; }
  float GetTowerThreshold() { return fgTowerThresh; }
  
#ifndef __CINT__
  void SetModules(std::vector<EmcModule> const *modules);
  std::vector<EmcModule> *GetModules(){ return fModules; }
  std::vector<EmcCluster> *GetClusters(){ return fClusters; }
#endif

  int iTowerDist(int ix1, int ix2);
  float fTowerDist(float x1, float x2);

  int FindClusters();
  void GetImpactAngle(float x, float y, float *sinT );
  void GlobalToSector(float, float, float, float*, float*, float*);
  void SectorToGlobal(float xsec, float ysec, float zsec, float* px,
		      float* py, float* pz );
  void SectorToGlobalErr( float dxsec, float dysec, float dzsec, float* pdx,
			  float* pdy, float* pdz );
  /// Converts coordinates in units of towers into cm's (Local coord. system)
  void TowersToSector(float, float, float &, float &);
  /// Returns  coordinates of the tower centers in cm's (Local coord. system)
  void TowersToSector(int,   int,   float &, float &);
  /// Converts Local Sector coordinates in cm into integer tower numbers
  void SectorToTowers(float, float, int &,   int &);

  void Gamma(int, EmcModule*, float*, float*, float*, float*, float*,
		     float*, float*, float*,
		     int &ndf); // ndf added MV 28.01.00
  void Momenta(int, EmcModule*, float*, float*, float*, float*, float*,
	       float* );
  int ShiftX(int ishift, int nh, EmcModule* phit0, EmcModule* phit1);
  EmcModule ShiftX(int ish, EmcModule ehit);

  void SetTowerThreshold(float Thresh);
  void SetProfileParameters(int, float, float, float);
  void SetChi2Limit(int lim);
  int GetTowerID( int iy, int iz, int nn, int* iyy, int* izz, float* ee );
  float GetProb(std::vector<EmcModule> HitList, float &chi2, int &ndf);
  float ClusterChisq(int, EmcModule*, float, float, float,
			     int &ndf); // ndf added MV 28.01.00
  float Chi2Correct(float chi2,int ndf);
  void CorrectPosition(float energy, float x, float y, float *xcorr,
				float *ycorr, bool callSetPar=true);
  void CorrectEnergy(float energy, float x, float y, float *ecorr);
  void CorrectECore(float ecore, float x, float y, float *ecorecorr);
  void TwoGamma(int, EmcModule*, float*, float*, float*, float*,
			float*, float*, float*);
  float Chi2Limit(int ndf);
  float PredictEnergy(float, float, float);
  void CalculateErrors(float e, float x, float y, float* pde,
			       float* pdx, float* pdy, float* pdz);
  void getTowerPos(int ix, int iy, float &x, float & y);

  // Auxiliary static functions
  static int HitNCompare(const void*, const void*);
  static int HitACompare(const void*, const void*);
  static void CopyVector(int*, int*, int);
  static void CopyVector(EmcModule*, EmcModule*, int);
  static void ZeroVector(int*, int);
  static void ZeroVector(float*, int);
  static void ZeroVector(EmcModule*, int);
  static void ResizeVector(int* , int, int);
  static void c3to5(float, float, float, float, float, float, float*, float*,
		    float*,  float*, float*, float*);

 protected:

  // geometry
  bool bCYL; // Cylindrical? 
  int fNx; // length in X direction
  int fNy; // length in Y direction
  float fModSizex; // module size in X direction (cm)
  float fModSizey; // module size in Y direction
  float fVx; // vertex position (cm)
  float fVy;
  float fVz;
  //  Sector - to - Global PHENIX transformation matrix
  //  PHMatrix emcrm;
  //  Sector - to - Global PHENIX translation vector
  //  PHVector emctr;
  //  Inverse transforamtion matrix and vector
  //  PHMatrix invemcrm;
  //  PHVector invemctr;

  // Tabulated values of Chi2 corresponding to 1% and 2% CL
  // for 50 values ndf (1-50)
  static float fgChi2Level1[];
  static float fgChi2Level2[];
  static float fgChi2Level[]; // actual level, chosen with SetChi2Limit().

#ifndef __CINT__
  std::vector<EmcModule> *fModules;
  std::vector<EmcCluster> *fClusters;
#endif

  float fgTowerThresh;
  float fgMinPeakEnergy;
  static float const fgMinShowerEnergy;
  static int const fgMaxLen;

  // From EmcScSectorRec
  //

  // Parameters for sigma calculation in Chi2 analysis
  static float fgEpar00;
  static float fgEpar0;
  static float fgEpar1;
  static float fgEpar2;
  static float fgEpar3;
  static float fgEpar4;

  // Parameters for shower shape and Chi2 calculation
  // Set by SetProfileParameters()
  float fSin4T;
  float fSinTx;
  float fSinTy;
  float fPpar1;
  float fPpar2;
  float fPpar3;
  float fPpar4;
  float fPshiftx;
  float fPshifty;
};

#endif // #ifndef BEMCREC_H

// ///////////////////////////////////////////////////////////////////////////
// EOF
