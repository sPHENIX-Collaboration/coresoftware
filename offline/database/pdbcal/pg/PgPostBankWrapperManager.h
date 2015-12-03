#ifndef __PGPOSTBANKWRAPPERMANAGER_HH_
#define __PGPOSTBANKWRAPPERMANAGER_HH__

#include <vector>

class PgPostBankWrapper;
class PgPostBankWrapper2;

class PgPostBankWrapperManager
{
 public:

  static PgPostBankWrapperManager& instance();

  void clear();

  bool commit();

  void print() const;

  bool registerWrapper(PgPostBankWrapper* wrapper);

  bool unregisterWrapper(PgPostBankWrapper* wrapper);

  bool registerWrapper2(PgPostBankWrapper2* wrapper);

  bool unregisterWrapper2(PgPostBankWrapper2* wrapper);

 private:

  PgPostBankWrapperManager();

  typedef std::vector<PgPostBankWrapper*> WVECTOR;

  WVECTOR fWrappers; //!

  typedef std::vector<PgPostBankWrapper2*> WVECTOR2;

  WVECTOR2 fWrappers2; //!

};

#endif
