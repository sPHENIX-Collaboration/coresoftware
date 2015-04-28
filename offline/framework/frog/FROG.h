#ifndef ROOT__FROG
#define ROOT__FROG

#include <string>

class FROG 
{

private:
  std::string pfn;
  //  char pfn[500];

public:
  FROG(){}
  virtual ~FROG(){}

  const char * location(const char * lname);
  const char * localSearch(const std::string &lname);
  void searchPG(const char * lname, std::string &cp);
  void searchDC(const char * lname, std::string &cp);
};

#endif
