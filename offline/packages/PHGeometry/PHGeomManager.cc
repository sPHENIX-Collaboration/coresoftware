
#include "PHGeomManager.h"

PHGeomManager *PHGeomManager::__instance = NULL;


PHGeomManager::PHGeomManager()
{
}

PHGeomManager::~PHGeomManager()
{
  __instance = NULL;
}

PHGeomManager *PHGeomManager::instance()
{
  if ( ! __instance )
    {
      std::cout << __FILE__ << "  " << __LINE__ << 
	" No instance of PHGeomManager available" << std::endl;
    }
  
  return __instance;
}

