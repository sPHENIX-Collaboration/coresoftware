// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_FLAGSAVE_H
#define FFAOBJECTS_FLAGSAVE_H

#include <phool/PHObject.h>
#include <phool/phool.h>

class PHFlag;

///
class FlagSave : public PHObject
{
 public:
  /// dtor
  ~FlagSave() override {}

  /// Clear Flag
  void Reset() override
  {
    std::cout << PHWHERE << "ERROR Reset() not implemented by daughter class" << std::endl;
    return;
  }

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream& os = std::cout) const override
  {
    os << "identify yourself: virtual FlagSave Object" << std::endl;
    return;
  }

  /// isValid returns non zero if object contains valid data
  int isValid() const override
  {
    std::cout << PHWHERE << "isValid not implemented by daughter class" << std::endl;
    return 0;
  }

  virtual int FillFromPHFlag(const PHFlag* /*flags*/) { return -1; }
  virtual int PutFlagsBack(PHFlag* /*flags*/) { return -1; }

 private:
  ClassDefOverride(FlagSave, 1)
};

#endif
