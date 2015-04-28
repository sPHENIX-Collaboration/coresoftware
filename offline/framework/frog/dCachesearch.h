#ifndef __DCACHESEARCH__H
#define __DCACHESEARCH__H

#include <string>

class dCachesearch
{

public:
  dCachesearch(){}
  virtual ~dCachesearch(){}

  // Make the search method public for external operators...
  void search(const std::string  &fileName, std::string &cp);

};

#endif
