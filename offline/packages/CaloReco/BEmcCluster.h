#ifndef CALORECO_EMCCLUSTER_H
#define CALORECO_EMCCLUSTER_H

// Name: EmcCluster.h
// Author: A. Bazilevsky (RIKEN-BNL)
// Major modifications by M. Volkov (RRC KI) Jan 27 2000

#include <TObject.h>

#include <cmath>
#include <cstdlib>
#include <vector>

// Forward declarations
class BEmcRec;

/** One tower information for internal clustering use.
@ingroup clustering
*/

class EmcModule
{
 public:
  EmcModule();
  EmcModule(int ich_, float amp_, float tof_);

  virtual ~EmcModule() {}

  int ich;    // module id (linear)
  float amp;  // module signal
  float tof;  // module time-of-flight
};

// ///////////////////////////////////////////////////////////////////////////

/** The 1-st level of the EMCal clustering: cluster is a set of contiguous
    towers. 

    Only used internally by clustering routines.
    @ingroup clustering
*/

class EmcCluster : public TObject
{
 public:
  /// Constructor (zero Hit List)
  EmcCluster()
    : fOwner(nullptr)
  {
  }

  explicit EmcCluster(BEmcRec* sector)
    : fOwner(sector)
  {
  }

  /// Constructor (inputs Hit List)

  EmcCluster(const std::vector<EmcModule>& hlist,
             BEmcRec* sector)
    : fHitList(hlist)
    , fOwner(sector)
  {
  }

  ///
  ~EmcCluster() override
  {
  }

  /// Reinitializes EmcCluster supplying new Hit List.

  void ReInitialize(const std::vector<EmcModule>& hlist)
  {
    fHitList = hlist;
  }
  /// Returns number of EmcModules in EmcCluster
  int GetNofHits() { return fHitList.size(); }
  /// Returns EmcCluster fHitList
  std::vector<EmcModule> GetHitList() { return fHitList; };
  /// Returns the EmcModule with the maximum energy
  EmcModule GetMaxTower();
  /// Returns the EmcModule corresponding to the reconstructed impact tower
  //  EmcModule GetImpactTower();
  /// Returns the energy of the ich-tower
  float GetTowerEnergy(int ich);
  /// Returns the energy of the tower ix,iy
  float GetTowerEnergy(int ix, int iy);
  /// Returns the ToF of the ich-tower
  float GetTowerToF(int ich);
  /// Returns the energy in 2x2 towers around the cluster Center of Gravity
  float GetE4();
  /// Returns the energy in 3x3 towers around the cluster Center of Gravity
  float GetE9();
  /// Returns the energy in 3x3 towers around the tower ich
  float GetE9(int ich);
  /// Returns the cluster energy taking into account towers with E>Ethresh
  float GetECore();
  /// Ecore corrected for energy leak sidewise core towers
  float GetECoreCorrected();
  /// Returns the EmcCluster total energy
  float GetTotalEnergy();
  /// Returns EmcCluster 1-st (pxcg,pycg) and 2-d momenta (pxx,pxy,pyy)
  void GetMoments(float& pxcg, float& pycg,
                  float& pxx, float& pxy, float& pyy);
  /// Returns the EmcCluster corrected position in Sector (SM) frame
  void GetCorrPos(float& xc, float& yc);
  /// Returns the EmcCluster position in PHENIX global coord system
  void GetGlobalPos(float& xg, float& yg, float& zg);
  /// Splits the Cluster onto SubClusters; returns list of clusters and list of peak towers corresponding to subclusters
  int GetSubClusters(std::vector<EmcCluster>& sClList, std::vector<EmcModule>& ppeaks);
  float GetProb(float& chi2, int& ndf);

 protected:
  std::vector<EmcModule> fHitList;

  BEmcRec* fOwner;  // what sector it belongs to

  // static members
  static int const fgMaxNofPeaks;
  static int const fgPeakIter;
  static float const fgEmin;

 public:
  // MV 2002/02/28 moved these functions here from #define's

  static int max(int a, int b)
  {
    return a > b ? a : b;
  }
  static float max(float a, float b)
  {
    return a > b ? a : b;
  }
  static double max(double a, double b)
  {
    return a > b ? a : b;
  }

  static int min(int a, int b)
  {
    return a < b ? a : b;
  }
  static float min(float a, float b)
  {
    return a < b ? a : b;
  }
  static double min(double a, double b)
  {
    return a < b ? a : b;
  }

  static int ABS(int x)
  {
    return abs(x);
  }
  static float ABS(float x)
  {
    return fabsf(x);
  }
  static double ABS(double x)
  {
    return fabs(x);
  }

  static int lowint(float x)
  {
    return x < 0. ? int(x - 1) : int(x);
  }
};

// ///////////////////////////////////////////////////////////////////////////

#endif
