#ifndef PDBCAL_BASE_PDBBANKLIST_H
#define PDBCAL_BASE_PDBBANKLIST_H

#include "PdbCalBank.h"

#include <phool/PHPointerList.h>

class PdbBankList : public PHPointerList<PdbCalBank>
{
public:
  PdbBankList() {}
  virtual ~PdbBankList() {}
};

#endif /* PDBCAL_BASE_PDBBANKLIST_H */
