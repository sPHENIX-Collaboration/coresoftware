#ifndef PHGeomManager_HH__
#define PHGeomManager_HH__


#include <ctime>
#include <map>
#include <set>
#include <string>
#include <iostream>

class TGeoManager;

class PHGeomManager
{

protected:

  PHGeomManager();
  virtual
  ~PHGeomManager();

public:

  static PHGeomManager *
  instance();

protected:



  static PHGeomManager *__instance;

};

#endif /* PHGeomManager_HH__ */
