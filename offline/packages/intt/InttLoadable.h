#ifndef INTT_LOADABLE_H
#define INTT_LOADABLE_H

#include <string>

class CDBTTree;

class InttLoadable
{
 public:
  virtual ~InttLoadable() = default;

  int LoadFromFile(std::string = "");
  int LoadFromCDB(std::string = "");
  int Loaded() const;

  virtual std::string DefaultFileName() const {return "";}
  virtual std::string DefaultCDBName() const {return "";}

 protected:
  InttLoadable() = default;

  virtual int LoadFromCDBTTree(CDBTTree&) = 0; // pure virtual

 private:
  int m_loaded{0};
};

#endif  // INTT_LOADABLE_H
