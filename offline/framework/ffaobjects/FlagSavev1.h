// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_FLAGSAVEV1_H
#define FFAOBJECTS_FLAGSAVEV1_H

#include "FlagSave.h"

#include <iostream>
#include <map>
#include <string>

class PHFlag;
class PHObject;

///
class FlagSavev1 : public FlagSave
{
 public:
  /// ctor
  FlagSavev1() = default;
  /// dtor
  ~FlagSavev1() override = default;

  PHObject *CloneMe() const override;

  ///  Clear Event
  void Reset() override {}
  int isValid() const override;

  /** identify Function from PHObject
      @param os Output Stream 
   */
  void identify(std::ostream &os = std::cout) const override;

  int FillFromPHFlag(const PHFlag *flags) override;
  int PutFlagsBack(PHFlag *flags) override;

 private:
  int FillIntFromPHFlag(const PHFlag *flags);
  int FillDoubleFromPHFlag(const PHFlag *flags);
  int FillFloatFromPHFlag(const PHFlag *flags);
  int FillCharFromPHFlag(const PHFlag *flags);

  int PutIntToPHFlag(PHFlag *flags);
  int PutDoubleToPHFlag(PHFlag *flags);
  int PutFloatToPHFlag(PHFlag *flags);
  int PutCharToPHFlag(PHFlag *flags);

  void PrintIntFlag(std::ostream &os) const;
  void PrintDoubleFlag(std::ostream &os) const;
  void PrintFloatFlag(std::ostream &os) const;
  void PrintStringFlag(std::ostream &os) const;

  std::map<std::string, int> intflag;
  std::map<std::string, double> doubleflag;
  std::map<std::string, float> floatflag;
  std::map<std::string, std::string> stringflag;

  ClassDefOverride(FlagSavev1, 1)

};

#endif
