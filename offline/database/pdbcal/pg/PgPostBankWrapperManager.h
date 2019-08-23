// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PDBCALPG_PGPOSTBANKWRAPPERMANAGER_H
#define PDBCALPG_PGPOSTBANKWRAPPERMANAGER_H

#include <vector>

class PgPostBankWrapper;

class PgPostBankWrapperManager
{
 public:
  static PgPostBankWrapperManager& instance();

  void clear();

  bool commit();

  void print() const;

  bool registerWrapper(PgPostBankWrapper* wrapper);

  bool unregisterWrapper(PgPostBankWrapper* wrapper);

 private:
  PgPostBankWrapperManager(){}

  typedef std::vector<PgPostBankWrapper*> WVECTOR;

  WVECTOR fWrappers;  //!
};

#endif
