// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef BBC_BBCOUT_H
#define BBC_BBCOUT_H

#include <phool/PHObject.h>

#include <string>

///
class BbcOut : public PHObject
{
 public:
  ///
  ~BbcOut() override {}

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  /// Clear Event
  virtual void Reset() override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  /// get ZVertex determined by Bbc
  virtual Float_t get_zvtx() const;

  /// get Error on ZVertex determined by Bbc
  virtual Float_t get_zvtxerr() const;

  /// get T0 determined by Bbc
  virtual Float_t get_t0() const;

  /// get Error on T0 determined by Bbc
  virtual Float_t get_t0err() const;

  /** set T0, Error on T0, ZVertex and Error on ZVertex
      @param t0 Bbc T0
      @param t0err Bbc Error on T0
      @param vtx Bbc ZVertex
      @param vtxerr Bbc Error on ZVertex
   */
  virtual void set_t0zvtx(const Float_t t0, const Float_t t0err, const Float_t vtx, const Float_t vtxerr)
  {
    set_t0(t0, t0err);
    set_zvtx(vtx, vtxerr);
  }

  /** set T0 for Bbc
      @param t0 Bbc T0
      @param t0err Bbc T0 error
   */
  virtual void set_t0(const Float_t t0, const Float_t t0err = 0);

  //! set vertex
  virtual void set_zvtx(const Float_t vtx, const Float_t vtxerr);

  /** set Vtx Error for Bbc
      @param vtxerr Bbc Vtx Error
   */
  virtual void set_zvtxerr(const Float_t vtxerr);

  /** Add Bbc North/South data containing Number of pmt's, Energy and Timing
      @param npmt Number of PMT's fired
      @param ncharge Number of Charged Particles into North/South
      @param timing Timing of North/South
      @param iarm  Arm, use Bbc::North and Bbc::South
   */
  virtual void set_arm(const int iarm, const Short_t npmt, const Float_t ncharge, const Float_t timing);

  /** Add Bbc data containing evt, clk, and femclk
      @param ievt   Event number
      @param iclk    XMIT clock
      @param ifemclk FEM clock
   */
  virtual void set_clocks(const Int_t ievt, const UShort_t iclk, const UShort_t ifemclk);

  /** get Number of PMT's fired in North/South Bbc
      @param iarm  Arm, use Bbc::North and Bbc::South
   */
  virtual Short_t get_npmt(const int iarm) const;

  /** get Number of Charged Particles into North/South Bbc
      @param iarm  Arm, use Bbc::North and Bbc::South
   */
  virtual Float_t get_q(const int iarm) const;

  /** get Timing of North/South Bbc
      @param iarm  Arm, use Bbc::North and Bbc::South
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

  virtual void FillFromClass(const BbcOut& old);

 private:
  void virtual_warning(const std::string& funcname) const;

  ClassDefOverride(BbcOut, 1)
};

#endif
