#ifndef __PGSEARCH__H
#define __PGSEARCH__H

#include <string>

class pgsearch
{

public:
  pgsearch(){}
  virtual ~pgsearch(){}

  // Make the search method public for external operators...
  void search(const std::string  &fileName, std::string &cp);

};

#endif
