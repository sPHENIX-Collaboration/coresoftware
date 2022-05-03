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

  /// Flags are read during InitRun() and written during End()
  /// Fills DST object with flags, if clearold is set, old flags from previous files
  /// which were deleted will not be saved
  virtual int FillFromPHFlag(const PHFlag* /*flags*/, const bool /* clearold */) { return -1; }
  /// Read back flags from the DST, if overwrite is set: flags from DST object will overwrite
  /// flag values set in the macro
  virtual int PutFlagsBack(PHFlag* /*flags*/, const bool /* overwrite */) { return -1; }

 private:
  ClassDefOverride(FlagSave, 1)
};

#endif
