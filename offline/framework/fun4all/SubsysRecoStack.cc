/*
 * a SubsysReco that can stack other modules.
 * (see header file for more info)
 *
 */

#include "SubsysRecoStack.h"
#include "Fun4AllReturnCodes.h"

#include <iomanip>
#include <string>
#include <typeinfo>


SubsysRecoStack::SubsysRecoStack(const char* name, PHCompositeNode * chroot): SubsysReco(name)
{
  this->chroot = chroot;
}


SubsysRecoStack::~SubsysRecoStack()
{
  while( !empty() )
    {
      delete back();
      pop_back();
    }
}




int SubsysRecoStack::Init(PHCompositeNode *topNode)
{
  PHCompositeNode * realroot = getroot( topNode );

  for( std::list<SubsysReco *>::iterator i = begin(); i != end(); ++i)
    {
      if( verbosity )
	std::cout << "calling " << (*i)->Name() << "::Init()" << std::endl;

      int rc = (*i)->Init(realroot);

      if( rc != Fun4AllReturnCodes::EVENT_OK )
	{
	  std::cout << "non EVENT_OK return value from " << (*i)->Name() << "::Init()" << std::endl;
	  return rc;
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}



int SubsysRecoStack::InitRun(PHCompositeNode *topNode)
{
  PHCompositeNode * realroot = getroot( topNode );

  for( std::list<SubsysReco *>::iterator i = begin(); i != end(); ++i)
    {
      if( verbosity )
	std::cout << "calling " << (*i)->Name() << "::InitRun()" << std::endl;

      int rc = (*i)->InitRun(realroot);

      if( rc != Fun4AllReturnCodes::EVENT_OK )
	{
	  std::cout << "non EVENT_OK return value from " << (*i)->Name() << "::InitRun()" << std::endl;
	  return rc;
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}



int SubsysRecoStack::process_event(PHCompositeNode *topNode)
{
  PHCompositeNode * realroot = getroot( topNode );

  for( std::list<SubsysReco *>::iterator i = begin(); i != end(); ++i)
    {
      if( verbosity )
	std::cout << "calling " << (*i)->Name() << "::process_event()" << std::endl;

      int rc = (*i)->process_event(realroot);

      if( rc != Fun4AllReturnCodes::EVENT_OK )
	{
	  std::cout << "non EVENT_OK return value from " << (*i)->Name() << "::process_event()" << std::endl;
	  return rc;
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}



int SubsysRecoStack::Reset(PHCompositeNode *topNode)
{
  PHCompositeNode * realroot = getroot( topNode );

  for( std::list<SubsysReco *>::iterator i = begin(); i != end(); ++i)
    {
      if( verbosity )
	std::cout << "calling " << (*i)->Name() << "::Reset()" << std::endl;

      int rc = (*i)->Reset(realroot);

      if( rc != Fun4AllReturnCodes::EVENT_OK )
	{
	  std::cout << "non EVENT_OK return value from " << (*i)->Name() << "::Reset()" << std::endl;
	  return rc;
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}



int SubsysRecoStack::ResetEvent(PHCompositeNode *topNode)
{
  PHCompositeNode * realroot = getroot( topNode );

  for( std::list<SubsysReco *>::iterator i = begin(); i != end(); ++i)
    {
      if( verbosity )
	std::cout << "calling " << (*i)->Name() << "::ResetEvent()" << std::endl;

      int rc = (*i)->ResetEvent(realroot);

      if( rc != Fun4AllReturnCodes::EVENT_OK )
	{
	  std::cout << "non EVENT_OK return value from " << (*i)->Name() << "::ResetEvent()" << std::endl;
	  return rc;
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}



int SubsysRecoStack::End(PHCompositeNode *topNode)
{
  PHCompositeNode * realroot = getroot( topNode );

  for( std::list<SubsysReco *>::iterator i = begin(); i != end(); ++i)
    {
      if( verbosity )
	std::cout << "calling " << (*i)->Name() << "::End()" << std::endl;

      int rc = (*i)->End(realroot);

      if( rc != Fun4AllReturnCodes::EVENT_OK )
	{
	  std::cout << "non EVENT_OK return value from " << (*i)->Name() << "::End()" << std::endl;
	  return rc;
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}



int SubsysRecoStack::EndRun(const int runnumber)
{

  for( std::list<SubsysReco *>::iterator i = begin(); i != end(); ++i)
    {
      if( verbosity )
	std::cout << "calling " << (*i)->Name() << "::EndRun()" << std::endl;

      int rc = (*i)->EndRun(runnumber);

      if( rc != Fun4AllReturnCodes::EVENT_OK )
	{
	  std::cout << "non EVENT_OK return value from " << (*i)->Name() << "::EndRun()" << std::endl;
	  return rc;
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}



void SubsysRecoStack::Print(const std::string &what) const
{
  Print2(what);
}



#define NAMEW(x) std::left << std::setw(24) << (x)->Name()
#define TYPEADDRW(x) "(" << std::left << std::setw(24) << typeid(*(x)).name() << " @ " << (x) << ")"

void SubsysRecoStack::Print2(const std::string &what, std::string prefix) const
{
    std::cout << prefix << NAMEW(this) << " " << TYPEADDRW(this) << ": ";
    std::cout << "chroot = " << chroot << std::endl;
    prefix += "    ";


    for( std::list<SubsysReco *>::const_iterator i = begin(); i != end(); ++i)
        {
            SubsysRecoStack * p = dynamic_cast< SubsysRecoStack * >( *i );

            if( p )
                {
                    p->Print2(what, prefix);
                }
            else
                {
                    std::cout << prefix << NAMEW(*i) << " " << TYPEADDRW(*i) << ": ";
                    (*i)->Print(what);
                }
        }

}



void SubsysRecoStack::Verbosity(const int ival)
{
    verbosity = ival;
    for( std::list<SubsysReco *>::iterator i = begin(); i != end(); ++i)
        (*i)->Verbosity(ival);
}


