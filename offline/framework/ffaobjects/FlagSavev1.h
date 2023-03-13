// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFAOBJECTS_FLAGSAVEV1_H
#define FFAOBJECTS_FLAGSAVEV1_H

#include "FlagSave.h"

#include <cstdint>  // for uint64_t
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

  int FillFromPHFlag(const PHFlag *flags, const bool clearold) override;
  int PutFlagsBack(PHFlag *flags, const bool overwrite) override;

 private:
  void ClearAll();
  int FillIntFromPHFlag(const PHFlag *flags);
  int Filluint64FromPHFlag(const PHFlag *flags);
  int FillDoubleFromPHFlag(const PHFlag *flags);
  int FillFloatFromPHFlag(const PHFlag *flags);
  int FillStringFromPHFlag(const PHFlag *flags);

  int PutIntToPHFlag(PHFlag *flags, const bool overwrite);
  int Putuint64ToPHFlag(PHFlag *flags, const bool overwrite);
  int PutDoubleToPHFlag(PHFlag *flags, const bool overwrite);
  int PutFloatToPHFlag(PHFlag *flags, const bool overwrite);
  int PutStringToPHFlag(PHFlag *flags, const bool overwrite);

  void PrintIntFlag(std::ostream &os) const;
  void Printuint64Flag(std::ostream &os) const;
  void PrintDoubleFlag(std::ostream &os) const;
  void PrintFloatFlag(std::ostream &os) const;
  void PrintStringFlag(std::ostream &os) const;

  std::map<std::string, int> intflag;
  std::map<std::string, double> doubleflag;
  std::map<std::string, float> floatflag;
  std::map<std::string, std::string> stringflag;
  std::map<std::string, uint64_t> m_uint64flag_map;

  ClassDefOverride(FlagSavev1, 2)
};

#endif
