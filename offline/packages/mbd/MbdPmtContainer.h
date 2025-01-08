// virtual Mbd PMT Container class

#ifndef __MBD_MBDPMTCONTAINER_H__
#define __MBD_MBDPMTCONTAINER_H__

#include <phool/PHObject.h>

#include <iostream>
#include <string>

class MbdPmtHit;

///
class MbdPmtContainer : public PHObject
{
 public:
  /// dtor
  virtual ~MbdPmtContainer() {}

  /** identify Function from PHObject
      @param os Output Stream
   */
  virtual void identify(std::ostream& os = std::cout) const override;

  /// Clear Event
  virtual void Reset() override;

  /// isValid returns non zero if object contains valid data
  virtual int isValid() const override;

  /** set number of PMTs for Mbd
      @param ival Number of Mbd Pmt's
   */
  virtual void set_npmt(const short ival);

  //! get Number of Mbd Pmt's
  virtual Short_t get_npmt() const;

  //! get MbdPmtHit Object
  virtual MbdPmtHit *get_pmt(const int ipmt) const;

 private:
  void virtual_warning(const std::string& funcname) const;

  ClassDefOverride(MbdPmtContainer, 1)
};

#endif  // __MBD_MBDPMTCONTAINER_H__
