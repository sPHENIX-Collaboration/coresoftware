#ifndef _EpFinder
#define _EpFinder

class TVector3;
class TProfile;
class TProfile2D;

#include "TH2D.h"
#include "TVector3.h"
#include "EpInfo.h"
#include <vector>
#include <utility>

/*************************************
 * \author Mike Lisa
 * \date 23 June 2018
 * 
 * adjusted for sPHENIX by J. Lajoie 
 * 24 August 2019
 *
 * \description:
 * Finds the Event Plane (EP) and EP-related quantities.
 * Also creates correction files and applies them.
 * Also calculates resolution
 *
 * There is a lot of EP-related information.  Raw, phi-weighted, shifted, etc.
 * 1st, 2nd, nth order.  Q-vector and Psi.  Even if the user does not request
 * all of these things, it is convenient and not so wasteful to simply calculate
 * them all at once.  Therefore, the user is presented with a large object of
 * type StEpdEpInfo.  This avoids "calculate-it-on-the-fly" which can be wasteful
 * if the user requests more than one thing, as well as "has-it-already-been-calculated?"
 * ambiguities.
 *
 *   A word about "EventType": the correction factors and other things about the event
 * plane can depend on centrality, vertex position, etc.  This code will apply the
 * corrections separately for different "EventTypes".  It is up to the user to decide
 * what this denotes.  All I care about is that when you send me an event, you tell
 * me the EventTypeId, which is just an integer.  The rest is up to you (as it should be :).
 *   For many (most?) people, this will just be a centrality index.
 *
 * This class will do "phi-weighting" and "shifting" corrections, but it needs
 *  some information to do it.  It will generate such information itself, but
 *  here is what you need to do, if you are starting from scratch:
 * 1) With whatever multiplicity/Vz/whatever event cuts you are going to analyze, and
 *    on whatever dataset, run your code that calls this class.  This will produce
 *    a file called EpFinderCorrectionHistograms_OUTPUT.root in the cwd.
 *    --> mv EpFinderCorrectionHistograms_OUTPUT.root EpFinderCorrectionHistos_INPUT.root
 *    That takes care of the Phi-weighting
 * 2) Repeat the above step (including the rename of the file, overwriting the old one).
 *      That takes care of the shifting weights.
 * 3) You are good to go.
 *
 *
 * ------------------------------------------
 * This class creates some histograms and uses some histograms.  Since I use GetBinContent()
 *  to extract values from the histograms, the binning is important.  I try to keep things
 *  consistent, but let me explain.
 *
 * 1) Phi Weighting.  Used by EpFinder and created by EpFinder.
 *    This code creates a histogram with the root name 
 *    "PhiWeight"
 *    x-axis and y-axis are USER DEFINED IN THE EpFinder constructor. 
 *    The user has no direct interaction with this histogram.
 *
 * 2) Shifting correction.  Used by EpFinder and created by EpFinder.
 *    This implements equation (6) of arxiv:nucl-ex/9805001
 *    The histogram names are
 *     - Form("EpdShiftdPsi%d_sin",ew,order)   
 *     - Form("EpdShiftdPsi%d_cos",ew,order)   
 *     - Form("EpdShiftFullEventPsi%d_sin",order)
 *     - Form("EpdShiftFullEventPsi%d_cos",order)
 *    In these histograms, order is "n" (as in n=2 for second-order EP)
 *    x-axis is "i" from equation (6) above.  As in <cos(n*i*Psi_n)>
 *       There are _EpTermsMax bins running from 0.5 to _EPtermsMax+0.5, so there should be no confusion with this axis.
 *    y-axis is EventTypeId, the *user-defined* EventType bin number.
 *--------->>>>>>>>>>>>>>>>>> And at this point I must make a demand of the user <<<<<<<<<<<<<<<<<<<
 *    When the user instantiates an EpFinder object, he specifies nEventTypeBins, the number of EventType bins he will use.
 *       >>>> The user MUST number these bins 0,1,2,3,4,...(nEventTypeBins-1) when he interacts with this class <<<<
 *       (If he wants to use a different convention in his code, that's fine, but when talking to EpFinder, use 0..(nEventTypeBins-1)
 *    The y-axis then has nEventTypeBins bins, going from -0.5 to nEventTypeBins-0.5
 *
 *************************************/

#define _EpTermsMax 6

// This is the structure for passing hit information into EpFinder:
// nMip : hit weight (energy of number of MIPs, etc.)
// ix   : detector index in x (user defined)
// iy   : detector index in y (user defined)
// samePhi : pointer is vector of index (ix,iy) pairs for decetor elements in the same phi bin as this hit (NULL is OK)
//
// ix, iy are USER DEFINED and their range is set in the EpFinder constructor. ix,iy and samePhi are used in the phi weighting correction determination. 

typedef struct{

  float nMip; 
  double phi;   
  int ix; 
  int iy; 
  std::vector<std::pair<int,int>> *samePhi; 
  
} EpHit; 


class EpFinder {
 public:

  /// Constructor.  Initializes values and reads correction file, if it exists.
  /// This file is actually PRODUCED by the code in an earlier run.  The user must rename
  /// the file EpFinderCorrectionHistograms_OUTPUT.root if he wants to use it.
  /// \param CorrectionFileName    Full name of the .root file with correction histograms.
  /// \param nEventTypeBins        Number of EventType bins that the user is using.  Up to the user to have a consistent usage, here and in analysis.
  /// \param pbinsx, pbinsy        Number of detector index bins for phi correction
  EpFinder(int nEventTypeBins=10, char const* OutFileName="EpFinderCorrectionHistograms_OUTPUT.root", char const* CorrectionFileName="EpFinderCorrectionHistograms_INPUT.root", 
	   int pbinsx=1, int pbinsy=1);
  ~EpFinder(){/* no-op */};

  /// sets the threshold, in units of nMIP, for determining tile weights
  ///   TileWeight = (EpdHit->nMIP()>thresh)?((EpdHit->nMIP()>MAX)?MAX:EpdHit->nMIP()):0;
  /// \param thresh       threshold.  If epdHit->nMIP() is less than thresh, then weight is zero
  void SetnMipThreshold(double thresh){mThresh=thresh;};

  /// sets the maximum weight, in units of nMIP, for determining tile weights
  ///   TileWeight = (EpdHit->nMIP()>thresh)?((EpdHit->nMIP()>MAX)?MAX:EpdHit->nMIP()):0;
  /// \param MAX          maximum tile weight.  If epdHit->nMIP()>MAX then weight=MAX
  void SetMaxTileWeight(double MAX){mMax=MAX;};

  /// call this method at the end of your run to output correction histograms to a file (you can choose to use these or not)
  ///   and to calculate EP resolutions
  void Finish();

  /// returns all information about the EP.  A large object of type StEpdEpInfo is returned, so you don't
  ///  have to call the EpFinder over and over again for various information
  /// \param EpdHits      Epd Hits in a TClones array.  Will be decoded as StEpdHit, StMuEpdHit, or StPicoEpdHit as dictated by mFormatUsed
  /// \param primVertex   primary vertex position for this event
  /// \param EventTypeID user-defined integer specifying EventType of the event.  User must use same convention in correction histograms and weights
  EpInfo Results(std::vector<EpHit> *EpdHits, int EventTypeID);

  /// Returns a big string that tells in text what the settings were.
  /// This is for your convenience and is of course optional.  I like
  /// to put a concatenation of such Reports into a text file, so I
  /// "autodocument" what were the settings for a given run
  TString Report();

 private:

  bool OrderOutsideRange(int order);         // just makes sure order is between 1 and _EpOrderMax

  double GetPsiInRange(double Qx, double Qy, int order);

  int mNumberOfEventTypeBins;                // user-defined.  Default is 10.  Used for correction histograms

  // tile weight = (0 if ADC< thresh), (MAX if ADC>MAX); (ADC otherwise).
  double mThresh;                            // default is 0.3
  double mMax;                               // default is 2.0

  TProfile* mAveCosDeltaPsi[_EpOrderMax];        // average of cos(Psi_{East,n}-Psi_{West,n}) using phi-weighted and shifted EPs

  std::string OutFileNameString; 

  //  these are shift correction factors that we MAKE now and write out
  TProfile2D* mEpShiftOutput_sin[_EpOrderMax];   
  TProfile2D* mEpShiftOutput_cos[_EpOrderMax];    
  //  these are shift correction factors that we made before, and USE now
  TProfile2D* mEpShiftInput_sin[_EpOrderMax];     
  TProfile2D* mEpShiftInput_cos[_EpOrderMax];     
  //   these are the phi weights
  TH2D* mPhiWeightInput;
//  TH2D* mPhiWeightOutput;
//  TH2D* mPhiAveraged;

};

#endif
