#ifndef EMCCLUSTER_H
#define EMCCLUSTER_H

// Name: EmcCluster.h
// Author: A. Bazilevsky (RIKEN-BNL)
// Major modifications by M. Volkov (RRC KI) Jan 27 2000

#include <TObject.h>

#include <cmath>

// Forward declarations
class BEmcRec;
class EmcCluster;
class EmcPeakarea;
class EmcEmshower;

/** One tower information for internal clustering use.
@ingroup clustering
*/

class EmcModule
{
 public:
  EmcModule();
  EmcModule(int ich_, int softkey_, float amp_, float tof_,
	    int deadmap_, int warnmap_, float adc_, float tac_);

  virtual ~EmcModule() {}

  int ich;   // module id (linear)
  int softKey; /* software key = arm/sector/row/column =
		  100000 * iarm +
		  10000 * iS +
		  100 * iy +
		  iz                       */
  float amp; // module signal
  float tof; // module time-of-flight
  int deadmap; // module dead map: see emcCalibratedDataObject.h
  int warnmap; // MV 2001/12/06
  float adc; // ADC amplitude
  float tac; // TAC amplitude
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
  {}

  EmcCluster(BEmcRec *sector): fOwner(sector)
  {}

  /// Constructor (inputs Hit List)
#ifndef __CINT__

  EmcCluster(const std::vector<EmcModule>& hlist,
	     BEmcRec *sector)
    : fOwner(sector)
  {
    fHitList = hlist;
  }
#endif

  ///
  virtual ~EmcCluster()
  {}

  /// Reinitializes EmcCluster supplying new Hit List.
#ifndef __CINT__

  void ReInitialize( const std::vector<EmcModule>& hlist )
  {
    fHitList = hlist;
  }
#endif
  /// Returns number of EmcModules in EmcCluster
  int GetNofHits()
  {
    return fHitList.size();
  }
  /// Returns n EmcModules (sorted) with the maximum energy
  void GetHits(EmcModule* phit, int n);
  /// Returns EmcCluster fHitList
#ifndef __CINT__

  std::vector<EmcModule> GetHitList()
  {
    return fHitList;
  };
#endif
  /// Returns the EmcModule with the maximum energy
  EmcModule GetMaxTower();
  /// Returns the EmcModule corresponding to the reconstructed impact tower
  EmcModule GetImpactTower();
  /// Returns the energy of the ich-tower (numbering from 0)
  float GetTowerEnergy( int ich );
  /// Returns the energy of the tower ix,iy (numbering from 0)
  float GetTowerEnergy( int ix, int iy );
  /// Returns the ToF of the ich-tower (numbering from 0)
  float GetTowerToF( int ich );
  /// Returns the Dead Map of the ich-tower (numbering from 0)
  int GetTowerDeadMap( int ich );
  /// Returns the Warning Map of the ich-tower (numbering from 0)
  int GetTowerWarnMap( int ich ); // MV 2002/02/18 bugfix
  float GetTowerADC( int ich ); // MV 2002/03/12 bugfix
  float GetTowerTAC( int ich ); // MV 2002/03/12 bugfix
  /// Returns the number of dead channels around MaxTower
  int GetNDead();
  int GetDeadMap(); // MV 2001/12/06
  int GetWarnMap(); // MV 2001/12/06
  /// Returns the energy in 2x2 towers around the cluster Center of Gravity
  float GetE4();
  /// Returns the energy in 3x3 towers around the cluster Center of Gravity
  float GetE9();
  /// Returns the energy in 3x3 towers around the tower ich
  float GetE9( int ich );
  /// Returns the cluster energy taking into account towers with E>Ethresh
  float GetECore();
  /// Ecore corrected for energy leak sidewise core towers
  float GetECoreCorrected();
  /// Returns the EmcCluster total energy
  float GetTotalEnergy();
  /// Returns EmcCluster 1-st (pxcg,pycg) and 2-d momenta (pxx,pxy,pyy)
  void GetMoments( float* pxcg, float* pycg,
		   float* pxx, float* pxy, float* pyy );
  /// Returns the EmcCluster corrected position in Sector (SM) frame
  void GetCorrPos( float* pxc, float* pyc );
  /// Returns the EmcCluster position in PHENIX global coord system
  void GetGlobalPos( float* pxg, float* pyg, float* pzg );
  /// Returns the errors for the reconstructed energy and position
  void GetErrors( float* pde, float* pdx, float* pdy, float* pdz);
  /// Substitutes a number of functions above (to save CPU time)
  void GetChar( float* pe,
		float* pxcg, float* pysg,
		float* pxc, float* pyc,
		float* pxg, float* pyg, float* pzg,
		float* pxx, float* pxy, float* pyy,
		float* pde, float* pdx, float* pdy, float* pdz );
  /// Splits the EmcCluster onto peakarea's; also returns peak tower array corresponding to peakarea array
  int GetPeaks( std::vector<EmcPeakarea> *PkList, std::vector<EmcModule> *ppeaks );
  float GetProb(float &chi2, int &ndf);

protected:

#ifndef __CINT__

  std::vector<EmcModule> fHitList;
#endif

  BEmcRec *fOwner; // what sector it belongs to

  // static members
  static int const fgMaxNofPeaks;
  static int const fgPeakIter;
  static float const fgEmin;
  static float const fgChisq;
  static float const fgXABSURD;
  static float const fgYABSURD;

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

/** The 2-d level of the EMCal clustering. 
    Every local maximum in cluster gives peakarea.
    @ingroup clustering
*/

class EmcPeakarea: public EmcCluster
{

public:

  /// Constructor (zero Hit List)
  EmcPeakarea(): fNdf(0), fCL(1.)
  {}

  EmcPeakarea(BEmcRec *sector):
    EmcCluster(sector), fNdf(0), fCL(1.)
  {}

  /// Constructor (inputs Hit List)
#ifndef __CINT__

  EmcPeakarea(const std::vector<EmcModule>& hlist, BEmcRec *sector):
    EmcCluster(hlist, sector), fNdf(0), fCL(1.)
  {}
#endif

  virtual ~EmcPeakarea()
  {}

  /// Returns Chi2
  float GetChi2();
  float GetChi2New();
  int GetNdf() const
  {
    return fNdf;
  } 
  // fetch number of degrees of freedom
  float GetCL() const
  {
    return fCL;
  } 
  // get CL (call only after GetChar())
  float GetCLNew();

  /// Returns peakarea's 1st momentum (COG) after Shower Shape fit
  void GetCGmin( float* pxcgmin, float* pycgmin );

  /// Substitutes a number of functions above (to save CPU time)
  void GetChar( float* pe, float* pec, float* pecore, float* pecorec,
		float* pxcg, float* pysg,
		float* pxcgmin, float* pysgmin,
		float* pxc, float* pyc,
		float* pxg, float* pyg, float* pzg,
		float* pxx, float* pxy, float* pyy,
		float* pchi,
		float* pde, float* pdx, float* pdy, float* pdz );

  /// Splits the peakarea onto 1 or 2 EmcEmshower's
  int GetGammas( EmcEmshower* );

protected:

  int fNdf; // Number of degrees of freedom
  float fCL; // Confidence level

};

// ///////////////////////////////////////////////////////////////////////////

/** The 3-d level of the EMCal clustering: peakarea with bad Chi2 is splitted
    onto two EmcEmshowers. 
    @ingroup clustering
*/
class EmcEmshower
{

public:

  /// Constructor (with energy and position zeroed)
  EmcEmshower();
  EmcEmshower(BEmcRec *sector);

  /// Constructor (inputs energy, position and Chi2)
  EmcEmshower(float e, float x, float y, float chi, int ndf,
	      BEmcRec *sector);

  /// Reinitializes EmcEmshower
  void ReInitialize(float e, float x, float y, float chi, int ndf)
  {
    fEnergy = e;
    fXcg = x;
    fYcg = y;
    fChisq = chi;
    fNdf = ndf;
  }

  /// Returns EmcEmshower total energy
  float GetTotalEnergy()
  {
    return fEnergy;
  }

  /// Returns EmcEmshower 1-at momentum (Center of Gravity)
  void GetCG(float* px, float* py)
  {
    *px = fXcg;
    *py = fYcg;
  }

  /// Returns EmcEmshower corrected position in Sector (SM) frame
  void GetCorrPos(float*, float* );

  /// Returns EmcEmshower position in PHENIX global coord system
  void GetGlobalPos(float*, float*, float*);

  /// Returns EmcEmshower Chi2
    float GetChi2() const
      {
        return fChisq;
      }

  /// Returns CL
  float GetCL() const
  {
    return fCL;
  } 
  // get CL (call only after GetChar())

  /// Returns the errors for the reconstructed energy and position
  void GetErrors( float* pde, float* pdx, float* pdy, float* pdz);
  
  /// Substitutes a number of functions above (to save CPU time)
  void GetChar( float* pe,
		float* pxcg, float* pysg,
		float* pxc, float* pyc,
		float* pxg, float* pyg, float* pzg,
		float* pchi,
		float* pde, float* pdx, float* pdy, float* pdz );
  
protected:

    BEmcRec* fOwner; // what sector it belongs to
    int fNdf; // Number of degrees of freedom
    float fCL; // Confidence level
    float fEnergy;
    float fXcg;
    float fYcg;
    float fChisq;

};

#endif // #ifdef EMCCLUSTER_H
