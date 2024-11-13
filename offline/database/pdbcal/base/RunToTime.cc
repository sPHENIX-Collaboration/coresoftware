#include "RunToTime.h"

#include <iostream>

RunToTime *RunToTime::__instance = nullptr;

RunToTime *RunToTime::instance()
{
  if (!__instance)
  {
    std::cout << __FILE__ << "  " << __LINE__ << " No instance of RunToTime available" << std::endl;
  }

  return __instance;
}
