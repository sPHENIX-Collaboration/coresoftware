#ifndef CALORECO_BEMCREC_H
#define CALORECO_BEMCREC_H

// Name: BEmcRec.h
// Author: A. Bazilevsky, Apr 2012
// Modified from EmcSectorRec.h and EmcScSectorRec.h

#include <map>
#include <vector>

class EmcCluster;
class EmcModule;
//class BEmcProfile;

typedef struct TowerGeom
{
  float Xcenter;  // Tower center position
  float Ycenter;
  float Zcenter;
  float dX[2]; // Tower i-th trans. dimension spread in global coord X
  float dY[2];
  float dZ[2];

} TowerGeom;

// ///////////////////////////////////////////////////////////////////////////

class BEmcRec
{
 public:
  BEmcRec();
  virtual ~BEmcRec();

  void SetVertex(float *vv)
  {
    fVx = vv[0];
    fVy = vv[1];
    fVz = vv[2];
  }
  void SetDim(int nx, int ny) { fNx = nx; fNy = ny; }

  bool SetTowerGeometry(int ix, int iy, float xx, float yy, float zz);
  bool GetTowerGeometry(int ix, int iy, TowerGeom &geom);
  bool CompleteTowerGeometry();
  void PrintTowerGeometry(const std::string &fname);

  void SetPlanarGeometry() { bCYL = false; }
  void SetCylindricalGeometry() { bCYL = true; }
  bool isCylindrical() const { return bCYL; }

  int GetNx() const { return fNx; }
  int GetNy() const { return fNy; }
  float GetVx() const { return fVx; }
  float GetVy() const { return fVy; }
  float GetVz() const { return fVz; }
  void SetPeakThreshold(float Thresh) { fgMinPeakEnergy = Thresh; }
  float GetPeakThreshold() { return fgMinPeakEnergy; }
  void SetTowerThreshold(float Thresh)  { fgTowerThresh = Thresh; }
  float GetTowerThreshold() { return fgTowerThresh; }

  void SetModules(std::vector<EmcModule> const *modules) { *fModules = *modules; }
  std::vector<EmcModule> *GetModules() { return fModules; }
  std::vector<EmcCluster> *GetClusters() { return fClusters; }

  int iTowerDist(int ix1, int ix2);
  float fTowerDist(float x1, float x2);

  int FindClusters();

  void Momenta(std::vector<EmcModule> *, float&, float&, float&, float&, float&,
               float&);

  void Tower2Global(float E, float xC, float yC, float &xA, float &yA, float &zA);

  virtual float PredictEnergy(float, float, float);

  // Calorimeter specific functions to be specified in respective inherited object
  virtual void CorrectEnergy(float energy, float x, float y, float* ecorr) {*ecorr=energy;}
  virtual void CorrectECore(float ecore, float x, float y, float* ecorecorr) {*ecorecorr=ecore;}
  virtual void CorrectPosition(float energy, float x, float y, float& xcorr, float& ycorr) {xcorr=x; ycorr=y;}
  virtual void CorrectShowerDepth(float energy, float x, float y, float z, float& xc, float& yc, float& zc ) {xc=x; yc=y; zc=z; }
  virtual void LoadProfile(const std::string &fname);
  virtual float GetProb(std::vector<EmcModule> HitList, float e, float xg, float yg, float zg, float &chi2, int &ndf);
  virtual std::string Name() const {return m_ThisName;}
  virtual void Name(const std::string &name) {m_ThisName = name;}

  // Auxiliary static functions
  static int HitNCompare(const void *, const void *);
  static int HitACompare(const void *, const void *);
  static void CopyVector(int *, int *, int);
  static void CopyVector(EmcModule *, EmcModule *, int);
  static void ZeroVector(int *, int);
  static void ZeroVector(float *, int);
  static void ZeroVector(EmcModule *, int);

 protected:
  // Geometry
  bool bCYL;  // Cylindrical?
  int fNx;    // length in X direction
  int fNy;    // length in Y direction
  std::map<int, TowerGeom> fTowerGeom;
  float fVx;        // vertex position (cm)
  float fVy;
  float fVz;

  std::vector<EmcModule> *fModules;
  std::vector<EmcCluster> *fClusters;

  float fgTowerThresh;
  float fgMinPeakEnergy;
  //  static float const fgMinShowerEnergy;
  static int const fgMaxLen;

  //  BEmcProfile *_emcprof;

 private:
  std::string m_ThisName;
  // the default copy ctor will not work
  // we do not use a copy ctor, so just delete it
  BEmcRec(const BEmcRec &) = delete;
};

#endif
