// virtual Mbd RAW Container class

#ifndef MBD_MBDRAWCONTAINER_H__
#define MBD_MBDRAWCONTAINER_H__

#include <phool/PHObject.h>

#include <iostream>
#include <string>

class MbdRawHit;

///
class MbdRawContainer : public PHObject
{
 public:
  /// dtor
  virtual ~MbdRawContainer() {}

  /** identify Function from PHObject
      @param os Output Stream
   */
  virtual void identify(std::ostream& os = std::cout) const override;

  /// Clear Event
  virtual void Reset() override;

  /// isValid returns non zero if object contains valid data
  virtual int isValid() const override;

  /** Add Mbd data containing evt, clk, and femclk
      @param ievt   Event number
      @param iclk    XMIT clock
      @param ifemclk FEM clock
   */
  virtual void set_clocks(const Int_t ievt, const UShort_t iclk, const UShort_t ifemclk);

  //! get Event Number
  virtual Int_t get_evt() const;

  //! get XMIT Clock Counter
  virtual UShort_t get_clock() const;

  //! get FEM Clock Counter
  virtual UShort_t get_femclock() const;

  /** set number of PMT's for Mbd
      @param ival Number of Mbd Pmt's
   */
  virtual void set_npmt(const short ival);

  //! get Number of Mbd Pmt's
  virtual Short_t get_npmt() const;

  //! get MbdRawHit Object
  virtual MbdRawHit *get_pmt(const int ipmt) const;

  virtual void Print(Option_t *option="") const override;

 private:
  static void virtual_warning(const std::string& funcsname) ;

  ClassDefOverride(MbdRawContainer, 1)
};

#endif  // __MBD_MBDRAWCONTAINER_H__
