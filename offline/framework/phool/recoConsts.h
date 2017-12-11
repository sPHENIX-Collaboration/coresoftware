// Do yourself and others a favour, please sort variable/function name
// according to the roman alphabet

#ifndef RECOCONSTS_H__
#define RECOCONSTS_H__

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

  void Print() const;

 protected:
  recoConsts() {}
  static recoConsts *__instance;
};

#endif /* __RECOCONSTS_H__ */
