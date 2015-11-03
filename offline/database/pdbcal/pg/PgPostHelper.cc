#include "PgPostHelper.hh"

#include "PHString.h"
#include "PHPointerList.h"

#include <algorithm>
#include <cctype>

std::string getTableName(const char* bankName)
  {
    PHString phs(bankName);
    PHPointerList<PHString> list;
    phs.split(list, ".");
    std::string tablename;
    for ( size_t i = 0; i < list.length(); ++i )
      {
        std::string part = list[i]->getString();
        std::transform(part.begin(), part.end(), part.begin(),
                       (int(*)(int))std::tolower);
	tablename += part;
	
      }
    list.clearAndDestroy();
    return tablename;
  }
