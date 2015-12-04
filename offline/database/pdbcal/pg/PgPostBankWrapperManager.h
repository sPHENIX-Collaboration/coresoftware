#ifndef __PGPOSTBANKWRAPPERMANAGER_HH_
#define __PGPOSTBANKWRAPPERMANAGER_HH__

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

  PgPostBankWrapperManager();

  typedef std::vector<PgPostBankWrapper*> WVECTOR;

  WVECTOR fWrappers; //!

};

#endif
