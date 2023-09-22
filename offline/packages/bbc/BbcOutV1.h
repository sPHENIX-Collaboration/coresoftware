// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCOUTV1_H
#define BBC_BBCOUTV1_H

#include "BbcOut.h"

#include <iostream>

class TClonesArray;

///
class BbcOutV1 : public BbcOut
{
 public:
  ///
  BbcOutV1();
  ///
  ~BbcOutV1() override;

  /// Clear Event from memory
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  /// get ZVertex determined by Bbc
  float get_zvtx() const override { return bz; }

  /// get Error on ZVertex determined by Bbc
  float get_zvtxerr() const override { return bzerr; }

  /// get T0 determined by Bbc
  float get_t0() const override { return bt0; }

  /// get Error on T0 determined by Bbc
  float get_t0err() const override { return bt0err; }

  /** set T0 for Bbc
      @param t0 Bbc T0
      @param t0err Bbc T0 error
   */
  void set_t0(const float t0, const float t0err = 0) override;

  //! set vertex
  void set_zvtx(const float vtx, const float vtxerr = 0) override;

  /** set Vtx Error for Bbc
      @param vtxerr Bbc Vtx Error
  */
  void set_zvtxerr(const float vtxerr) override;

  /** Add Bbc North/South data containing Number of pmt's, Energy and Timing
      @param npmt Number of PMT's fired
      @param energy Energy in North/South
      @param timing Timing of North/South
      @param iarm  Arm, use Bbc::North and Bbc::South
   */
  void set_arm(const int iarm, const short npmt, const float chargesum, const float timing) override;

  /** get Number of PMT's fired in North/South Bbc
      @param iarm  Arm, use Bbc::North and Bbc::South
   */
  short get_npmt(const int iarm) const override;

  /** get Number of Charged Particles into North/South Bbc
      @param iarm  Arm, use Bbc::North and Bbc::South
   */
  float get_q(const int iarm) const override;

  /** get Timing of North/South Bbc
      @param iarm  Arm, use Bbc::North and Bbc::South
   */
  float get_time(const int iarm) const override;


 private:

  void Init();

  float bz{};
  float bzerr{};
  float bt0{};
  float bt0err{};
  float bns{};
  float bnn{};
  float bqs{};
  float bqn{};
  float bts{};
  float btn{};

  ClassDefOverride(BbcOutV1, 1)
};

#endif
