// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MBD_MBDOUT_H
#define MBD_MBDOUT_H

#include <phool/PHObject.h>

#include <string>

///
class MbdOut : public PHObject
{
 public:
  ///
  ~MbdOut() override {}

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  /// Clear Event
  virtual void Reset() override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  /// get ZVertex determined by Mbd
  virtual Float_t get_zvtx() const;

  /// get Error on ZVertex determined by Mbd
  virtual Float_t get_zvtxerr() const;

  /// get T0 determined by Mbd
  virtual Float_t get_t0() const;

  /// get Error on T0 determined by Mbd
  virtual Float_t get_t0err() const;

  /** set T0, Error on T0, ZVertex and Error on ZVertex
      @param t0 Mbd T0
      @param t0err Mbd Error on T0
      @param vtx Mbd ZVertex
      @param vtxerr Mbd Error on ZVertex
   */
  virtual void set_t0zvtx(const Float_t t0, const Float_t t0err, const Float_t vtx, const Float_t vtxerr)
  {
    set_t0(t0, t0err);
    set_zvtx(vtx, vtxerr);
  }

  /** set T0 for Mbd
      @param t0 Mbd T0
      @param t0err Mbd T0 error
   */
  virtual void set_t0(const Float_t t0, const Float_t t0err = 0);

  //! set vertex
  virtual void set_zvtx(const Float_t vtx, const Float_t vtxerr);

  /** set Vtx Error for Mbd
      @param vtxerr Mbd Vtx Error
   */
  virtual void set_zvtxerr(const Float_t vtxerr);

  /** Add Mbd North/South data containing Number of pmt's, Energy and Timing
      @param npmt Number of PMT's fired
      @param ncharge Number of Charged Particles into North/South
      @param timing Timing of North/South
      @param iarm  Arm, use Mbd::North and Mbd::South
   */
  virtual void set_arm(const int iarm, const Short_t npmt, const Float_t ncharge, const Float_t timing);

  /** Add Mbd data containing evt, clk, and femclk
      @param ievt   Event number
      @param iclk    XMIT clock
      @param ifemclk FEM clock
   */
  virtual void set_clocks(const Int_t ievt, const UShort_t iclk, const UShort_t ifemclk);

  /** get Number of PMT's fired in North/South Mbd
      @param iarm  Arm, use Mbd::North and Mbd::South
   */
  virtual Short_t get_npmt(const int iarm) const;

  /** get Number of Charged Particles into North/South Mbd
      @param iarm  Arm, use Mbd::North and Mbd::South
   */
  virtual Float_t get_q(const int iarm) const;

  /** get Timing of North/South Mbd
      @param iarm  Arm, use Mbd::North and Mbd::South
   */
  virtual Float_t get_time(const int iarm) const;

  /** get Event Number
   */
  virtual Int_t get_evt() const;

  /** get XMIT Clock Counter
   */
  virtual UShort_t get_clock() const;

  /** get FEM Clock Counter
   */
  virtual UShort_t get_femclock() const;

  virtual void FillFromClass(const MbdOut& old);

 private:
  void virtual_warning(const std::string& funcname) const;

  ClassDefOverride(MbdOut, 1)
};

#endif
