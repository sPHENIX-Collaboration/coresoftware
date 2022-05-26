#ifndef EpInfo_H
#define EpInfo_H

#define _EpOrderMax 3   // maximum order of EP to worry about.

#include "TVector2.h"

#include <phool/PHObject.h>

/// the class EpInfo has only public members.
/// No need to hide information.  Indeed, it's really only a struct,
/// with the possibility to create some methods.
class EpInfo : public PHObject {

  // making EpFinder a "friend class" just gives it direct access to the arrays.
  // But the general User is required to use accessors, since the numbering/index convention can be confusing
  friend class EpFinder;

 public:
  EpInfo();
  ~EpInfo(){/* no op */};

  void Reset(); 

  // in the below, when it says "of order," then "order" begins at 1.  E.g. order=2 means second-order q vector

  //-----------------------------------------------------------------------------------------
  /// Raw (no phi-weighting) Q vector 
  /// \parameter order     order of the Q-vector.  Begins at unity (order=1 means first-order Q)
  TVector2 RawQ(int order);
 
  //-----------------------------------------------------------------------------------------
  /// Phi weighted Q vector 
  /// \parameter order     order of the Q-vector.  Begins at unity (order=1 means first-order Q)
  TVector2 PhiWeightedQ(int order);

  //-----------------------------------------------------------------------------------------
  /// Raw (no phi-weighting, no shifting) Event plane angle 
  /// \parameter order     order of the Q-vector.  Begins at unity (order=1 means first-order Q)
  double RawPsi(int order);

  //-----------------------------------------------------------------------------------------
  /// Phi-weighted (but not shift flattened) Event plane angle 
  /// \parameter order     order of the Q-vector.  Begins at unity (order=1 means first-order Q)
  double PhiWeightedPsi(int order);

  //-----------------------------------------------------------------------------------------
  /// Phi-weighted and shift-corrected Event plane angle 
  /// \parameter order     order of the Q-vector.  Begins at unity (order=1 means first-order Q)
  double PhiWeightedAndShiftedPsi(int order);

  //-----------------------------------------------------------------------------------------
  /// The sum of weights used to calculate Q-vector, for entire detector.  This is RAW, not phi-weighted
  /// This is useful if one wants to un-normalize the Q-vector
  /// ** note that this depends on "order," because the eta-weighting or ring-weighting can depend on order (user sets it)
  /// ** this stands in contrast to the ring-by-ring Qvectors, which do not depend on order.
  /// \parameter order     order of the EP.  Begins at unity (order=1 means first-order EP)
  double SWRaw(int order);

  //-----------------------------------------------------------------------------------------
  /// The sum of weights used to calculate Q-vector, for entire detector.  This is phi-weighted
  /// This is useful if one wants to un-normalize the Q-vector
  /// ** note that this depends on "order," because the eta-weighting or ring-weighting can depend on order (user sets it)
  /// ** this stands in contrast to the ring-by-ring Qvectors, which do not depend on order.
  /// \parameter order     order of the EP.  Begins at unity (order=1 means first-order EP)
  double SWPhiWeighted(int order);

  double PsiRaw[_EpOrderMax];                   /// indices: [order]

 private:

  bool ArgumentOutOfBounds(int order);              /// protection against user selecting "order=0" or order that has not been defined

  double Range(double psi, int order);                     /// puts angle psi into range (0,2pi/n)

  double QrawOneSide[_EpOrderMax][2];           /// indices: [order][x,y]
  double QphiWeightedOneSide[_EpOrderMax][2];   /// indices: [order][x,y]
  double PsiPhiWeighted[_EpOrderMax];           /// indices: [order]
  double PsiPhiWeightedAndShifted[_EpOrderMax]; /// indices: [order]
  double WheelSumWeightsRaw[_EpOrderMax];       /// indices: [order]
  double WheelSumWeightsPhiWeighted[_EpOrderMax]; /// indices: [order]

};

#endif
