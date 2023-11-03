// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MBD_MBDOUTV1_H
#define MBD_MBDOUTV1_H

#include "MbdOut.h"

#include <iostream>
#include <limits>

class TClonesArray;

///
class MbdOutV1 : public MbdOut
{
public:
  ///
  MbdOutV1();
  ///
  ~MbdOutV1() override;

  /// Clear Event from memory
  virtual void Reset() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream &os = std::cout) const override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  /// get ZVertex determined by Mbd
  Float_t get_zvtx() const override { return bz; }

  /// get Error on ZVertex determined by Mbd
  Float_t get_zvtxerr() const override { return bzerr; }

  /// get T0 determined by Mbd
  Float_t get_t0() const override { return bt0; }

  /// get Error on T0 determined by Mbd
  Float_t get_t0err() const override { return bt0err; }

  /** set T0 for Mbd
      @param t0 Mbd T0
      @param t0err Mbd T0 error
   */
  void set_t0(const Float_t t0, const Float_t t0err = 0) override;

  //! set vertex
  void set_zvtx(const Float_t vtx, const Float_t vtxerr = 0) override;

  /** set Vtx Error for Mbd
      @param vtxerr Mbd Vtx Error
  */
  void set_zvtxerr(const Float_t vtxerr) override;

  /** Add Mbd North/South data containing Number of pmt's, Energy and Timing
      @param npmt Number of PMT's fired
      @param energy Energy in North/South
      @param timing Timing of North/South
      @param iarm  Arm, use Mbd::North and Mbd::South
   */
  void set_arm(const int iarm, const Short_t npmt, const Float_t chargesum, const Float_t timing) override;

  /** get Number of PMT's fired in North/South Mbd
      @param iarm  Arm, use Mbd::North and Mbd::South
   */
  Short_t get_npmt(const int iarm) const override;

  /** get Number of Charged Particles into North/South Mbd
      @param iarm  Arm, use Mbd::North and Mbd::South
   */
  Float_t get_q(const int iarm) const override;

  /** get Timing of North/South Mbd
      @param iarm  Arm, use Mbd::North and Mbd::South
   */
  Float_t get_time(const int iarm) const override;


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

  ClassDefOverride(MbdOutV1, 1)
};

#endif
