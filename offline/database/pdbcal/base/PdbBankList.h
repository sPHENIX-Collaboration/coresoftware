#ifndef __PDBBANKLIST_HH__
#define __PDBBANKLIST_HH__

#include "PdbCalBank.h"

#include <phool/PHPointerList.h>

class PdbBankList : public PHPointerList<PdbCalBank>
{
public:
  PdbBankList() {}
  virtual ~PdbBankList() {}
};

#endif /* __PDBBANKLIST_HH__ */
