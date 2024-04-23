#ifndef INTT_LOADABLE_H
#define INTT_LOADABLE_H

#include <string>

class CDBTTree;

class InttLoadable
{
 public:
  InttLoadable() = default;
  virtual ~InttLoadable() = default;

  int Load(std::string = "");
  bool IsLoaded() const;

  int Verbosity() const;
  int Verbosity(int const&);

  virtual std::string DefaultName() const {return "";}

 protected:
  virtual int LoadFromCDBTTree(CDBTTree&);
  int m_verbosity{0};

 private:
  bool m_loaded{0};
};

#endif  // INTT_LOADABLE_H
