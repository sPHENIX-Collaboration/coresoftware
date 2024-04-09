#ifndef INTT_LOADABLE_H
#define INTT_LOADABLE_H

#include <string>

class CDBTTree;

class InttLoadable
{
 public:
  virtual ~InttLoadable() = default;

  int LoadFromFile(std::string const&);
  int LoadFromCDB(std::string const&);

  virtual int LoadFromFile();
  virtual int LoadFromCDB();

 protected:
  InttLoadable() = default;
  virtual int LoadFromCDBTTree(CDBTTree&);
};

#endif  // INTT_LOADABLE_H
