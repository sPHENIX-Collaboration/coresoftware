// virtual Bbc PMT Container class

#ifndef __BBC_BBCPMTCONTAINER_H__
#define __BBC_BBCPMTCONTAINER_H__

#include <phool/PHObject.h>

#include <iostream>
#include <string>

class BbcPmtHit;

///
class BbcPmtContainer : public PHObject
{
 public:
  /// dtor
  virtual ~BbcPmtContainer() {}

  /** identify Function from PHObject
      @param os Output Stream
   */
  virtual void identify(std::ostream& os = std::cout) const override;

  /// Clear Event
  virtual void Reset() override;

  /// isValid returns non zero if object contains valid data
  virtual int isValid() const override;

  /** set number of PMTs for Bbc
      @param ival Number of Bbc Pmt's
   */
  virtual void set_npmt(const short ival);

  //! get Number of Bbc Pmt's
  virtual Short_t get_npmt() const;

  //! get BbcPmtHit Object
  virtual BbcPmtHit *get_pmt(const int ipmt) const;

 private:
  void virtual_warning(const std::string& funcname) const;

  ClassDefOverride(BbcPmtContainer, 1)
};

#endif  // __BBC_BBCPMTCONTAINER_H__
