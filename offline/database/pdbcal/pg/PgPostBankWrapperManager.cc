#include <PgPostBankWrapperManager.h>
#include <PgPostBankWrapper.h>
#include <PgPostBankWrapper2.h>

#include <phool/phool.h>

#include <algorithm>
#include <iomanip>
#include <iostream>

//_____________________________________________________________________________
PgPostBankWrapperManager::PgPostBankWrapperManager()
{
}

//_____________________________________________________________________________
PgPostBankWrapperManager&
PgPostBankWrapperManager::instance()
{
  static PgPostBankWrapperManager instance__;
  return instance__;
}

//_____________________________________________________________________________
void
PgPostBankWrapperManager::clear()
{
  // We simply clear the vector. Does not care about the pointers,
  // as we are *not* owners.
  fWrappers.clear();
}

//_____________________________________________________________________________
bool
PgPostBankWrapperManager::commit()
{
  WVECTOR::const_iterator it;

  for ( it = fWrappers.begin(); it != fWrappers.end(); ++it )
    {
      (*it)->commit();
    }
  // Once committed, forget about the registered wrappers.
  clear();
  return true;
}

//_____________________________________________________________________________
void
PgPostBankWrapperManager::print() const
{
  WVECTOR::const_iterator it;
  size_t i = 0;

  for ( it = fWrappers.begin(); it != fWrappers.end(); ++it )
    {
      PgPostBankWrapper* wrapper = *it;
      if (!wrapper)
	{
	  std::cout << PHWHERE << " invalid wrapper pointer in list, exiting" << std::endl;
	  exit(1);
	}
      std::cout << "Wrapper " << std::setw(4) << i
		<< " " << " is " << std::hex << wrapper
		<< " BankID=" << std::dec 
		<< wrapper->getBankID().getInternalValue()
		<< " Table=" << wrapper->getTableName()
		<< std::endl;
      ++i;
    }
}

//_____________________________________________________________________________
bool
PgPostBankWrapperManager::registerWrapper(PgPostBankWrapper* wrapper)
{
  if (!wrapper)
    {
      std::cerr << PHWHERE << "Cannot register a null wrapper!!!"
		<< std::endl;
      return false;
    }

  WVECTOR::const_iterator it = 
    std::find(fWrappers.begin(),fWrappers.end(),wrapper);

  if ( it != fWrappers.end() )
    {
      std::cerr << PHWHERE << "This wrapper is already registered!!!"
		<< std::endl;
      return false;
    }
  fWrappers.push_back(wrapper);
  return true;
}
bool
PgPostBankWrapperManager::registerWrapper2(PgPostBankWrapper2* wrapper2)
{
  if (!wrapper2)
    {
      std::cerr << PHWHERE << "Cannot register a null wrapper2!!!"
		<< std::endl;
      return false;
    }

  WVECTOR2::const_iterator it = 
    std::find(fWrappers2.begin(),fWrappers2.end(),wrapper2);

  if ( it != fWrappers2.end() )
    {
      std::cerr << PHWHERE << "This wrapper2 is already registered!!!"
		<< std::endl;
      return false;
    }
  fWrappers2.push_back(wrapper2);
  return true;
}
//_____________________________________________________________________________
bool
PgPostBankWrapperManager::unregisterWrapper(PgPostBankWrapper* wrapper)
{ 
  if (!wrapper)
    {
      std::cerr << PHWHERE << "Cannot unregister a null wrapper!!!"
		<< std::endl;
      return false;
    }
  
  WVECTOR::iterator it = 
    std::find(fWrappers.begin(),fWrappers.end(),wrapper);

  if ( it != fWrappers.end() )
    {
      fWrappers.erase(it);
      return true;
    }
  else
    {
      //      std::cerr << PHWHERE << "This wrapper cannot be unregistered, as it's "
      //	<< " not there !!!!"
      //	<< std::endl;
      return false;
    }
}


bool
PgPostBankWrapperManager::unregisterWrapper2(PgPostBankWrapper2* wrapper2)
{ 
  if (!wrapper2)
    {
      std::cerr << PHWHERE << "Cannot unregister a null wrapper2!!!"
		<< std::endl;
      return false;
    }
  
  WVECTOR2::iterator it = 
    std::find(fWrappers2.begin(),fWrappers2.end(),wrapper2);

  if ( it != fWrappers2.end() )
    {
      fWrappers2.erase(it);
      return true;
    }
  else
    {
      //      std::cerr << PHWHERE << "This wrapper2 cannot be unregistered, as it's "
      //	<< " not there !!!!"
      //	<< std::endl;
      return false;
    }
}
