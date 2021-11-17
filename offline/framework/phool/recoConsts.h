// Do yourself and others a favour, please sort variable/function name
// according to the roman alphabet

#ifndef PHOOL_RECOCONSTS_H
#define PHOOL_RECOCONSTS_H

#include "PHFlag.h"

class recoConsts : public PHFlag
{
 public:
  static recoConsts *instance()
  {
    if (__instance) return __instance;
    __instance = new recoConsts();
    return __instance;
  }

  void Print() const override;

 private:
  recoConsts() {}
  static recoConsts *__instance;
};

#endif
