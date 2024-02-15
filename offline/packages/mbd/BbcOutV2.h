// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCOUTV2_H
#define BBC_BBCOUTV2_H

#include "BbcOut.h"

#include <iostream>
#include <limits>

class TClonesArray;

///
class BbcOutV2 : public BbcOut
{
 public:
  ///
  BbcOutV2();
  ///
  ~BbcOutV2() override;

  /// Clear Event from memory
  void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  /// get ZVertex determined by Bbc
  Float_t get_zvtx() const override { return bz; }

  /// get Error on ZVertex determined by Bbc
  Float_t get_zvtxerr() const override { return bzerr; }

  /// get T0 determined by Bbc
  Float_t get_t0() const override { return bt0; }

  /// get Error on T0 determined by Bbc
  Float_t get_t0err() const override { return bt0err; }

  /** set T0 for Bbc
      @param t0 Bbc T0
      @param t0err Bbc T0 error
   */
  void set_t0(const Float_t t0, const Float_t t0err = 0) override;

  //! set vertex
  void set_zvtx(const Float_t vtx, const Float_t vtxerr = 0) override;

  /** set Vtx Error for Bbc
      @param vtxerr Bbc Vtx Error
  */
  void set_zvtxerr(const Float_t vtxerr) override;

  /** Add Bbc North/South data containing Number of pmt's, Energy and Timing
      @param npmt Number of PMT's fired
      @param energy Energy in North/South
      @param timing Timing of North/South
      @param iarm  Arm, use Bbc::North and Bbc::South
   */
  void set_arm(const int iarm, const Short_t npmt, const Float_t chargesum, const Float_t timing) override;

  /** Add Bbc data containing evt, clk, and femclk
      @param ievt   Event number
      @param iclk    XMIT clock
      @param ifemclk FEM clock
   */
  virtual void set_clocks(const Int_t ievt, const UShort_t iclk, const UShort_t ifemclk) override;

  /** get Number of PMT's fired in North/South Bbc
      @param iarm  Arm, use Bbc::North and Bbc::South
   */
  Short_t get_npmt(const int iarm) const override;

  /** get Number of Charged Particles into North/South Bbc
      @param iarm  Arm, use Bbc::North and Bbc::South
   */
  Float_t get_q(const int iarm) const override;

  /** get Timing of North/South Bbc
      @param iarm  Arm, use Bbc::North and Bbc::South
   */
  Float_t get_time(const int iarm) const override;

  /** get Event Number
   */
  virtual Int_t get_evt() const override;

  /** get XMIT Clock Counter
   */
  virtual UShort_t get_clock() const override;

  /** get FEM Clock Counter
   */
  virtual UShort_t get_femclock() const override;

 private:
  Float_t bz{std::numeric_limits<Float_t>::quiet_NaN()};
  Float_t bzerr{std::numeric_limits<Float_t>::quiet_NaN()};
  Float_t bt0{std::numeric_limits<Float_t>::quiet_NaN()};
  Float_t bt0err{std::numeric_limits<Float_t>::quiet_NaN()};
  Float_t bqs{0};
  Float_t bqn{0};
  Float_t bts{std::numeric_limits<Float_t>::quiet_NaN()};
  Float_t btn{std::numeric_limits<Float_t>::quiet_NaN()};
  Short_t bns{0};
  Short_t bnn{0};
  Int_t evt{-1};
  UShort_t clk{0};
  UShort_t femclk{0};

  ClassDefOverride(BbcOutV2, 1)
};

#endif
