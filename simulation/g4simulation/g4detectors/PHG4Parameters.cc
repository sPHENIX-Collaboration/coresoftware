#include "PHG4Parameters.h"

#include <pdbcalbase/PdbBankManager.h>
#include <pdbcalbase/PdbApplication.h>
#include <pdbcalbase/PdbBankList.h>
#include <pdbcalbase/PdbCalBank.h>
#include <pdbcalbase/PdbParameterMap.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHTimeStamp.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

void
PHG4Parameters::set_int_param(const std::string &name, const int ival)
{
  intparams[name] = ival;
}

int
PHG4Parameters::get_int_param(const std::string &name) const
{
  if (intparams.find(name) != intparams.end())
    {
      return intparams.find(name)->second;
    }
  cout << PHWHERE << " integer parameter " << name << " does not exist (forgot to set?)" << endl;

  exit(1);
}

void
PHG4Parameters::printint() const
{
  cout << "int parameters: " << endl;
  for (map<const string, int>::const_iterator iter = intparams.begin(); iter != intparams.end(); ++iter)
    {
      cout << iter->first << ": " << iter->second << endl;
    }
  return;
}

void
PHG4Parameters::set_double_param(const std::string &name, const double dval)
{
  doubleparams[name] = dval;
}

double
PHG4Parameters::get_double_param(const std::string &name) const
{
  if (doubleparams.find(name) != doubleparams.end())
    {
      return doubleparams.find(name)->second;
    }
  cout << PHWHERE << " double parameter " << name << " does not exist (forgot to set?)" << endl;

  exit(1);
}

void
PHG4Parameters::print() const
{
  cout << "Parameters for " << detname << endl;
  printint();
  printdouble();
  printstring();
  return;
}

void
PHG4Parameters::printdouble() const
{
  cout << "double parameters: " << endl;
  for (map<const string, double>::const_iterator iter = doubleparams.begin(); iter != doubleparams.end(); ++iter)
    {
      cout << iter->first << ": " << iter->second << endl;
    }
  return;
}

void
PHG4Parameters::set_string_param(const std::string &name, const string &str)
{
  stringparams[name] = str;
}

string
PHG4Parameters::get_string_param(const std::string &name) const
{
  if (stringparams.find(name) != stringparams.end())
    {
      return stringparams.find(name)->second;
    }
  cout << PHWHERE << " string parameter " << name << " does not exist (forgot to set?)" << endl;

  exit(1);
}

void
PHG4Parameters::printstring() const
{
  cout << "string parameters: " << endl;
  for (map<const string, string>::const_iterator iter = stringparams.begin(); iter != stringparams.end(); ++iter)
    {
      cout << iter->first << ": " << iter->second << endl;
    }
  return;
}


void
PHG4Parameters::FillFrom(const PdbParameterMap *saveparams)
{
  pair<std::map<const std::string, double>::const_iterator, std::map<const std::string, double>::const_iterator> begin_end_d = saveparams->get_dparam_iters();
  for (map<const std::string, double>::const_iterator iter = begin_end_d.first; iter != begin_end_d.second;++iter)
    {
      doubleparams[iter->first] = iter->second;
    }
  pair<std::map<const std::string, int>::const_iterator, std::map<const std::string, int>::const_iterator> begin_end_i = saveparams->get_iparam_iters();
  for (map<const std::string, int>::const_iterator iter = begin_end_i.first; iter != begin_end_i.second;++iter)
    {
      intparams[iter->first] = iter->second;
    }
  pair<std::map<const std::string, string>::const_iterator, std::map<const std::string, string>::const_iterator> begin_end_s = saveparams->get_cparam_iters();
  for (map<const std::string, string>::const_iterator iter = begin_end_s.first; iter != begin_end_s.second;++iter)
    {
      stringparams[iter->first] = iter->second;
    }

  return;
}


void
PHG4Parameters::SaveToNodeTree(PHCompositeNode *topNode, const string &nodename)
{
  // write itself since this class is fine with saving by root
  PdbParameterMap *nodeparams = findNode::getClass<PdbParameterMap>(topNode,nodename);
  if (!nodeparams)
    {
      nodeparams = new PdbParameterMap();
      PHIODataNode<PdbParameterMap> *newnode =  new PHIODataNode<PdbParameterMap>(nodeparams,nodename);
      topNode->addNode(newnode);
    }
  for (map<const string,double>::const_iterator iter = doubleparams.begin(); iter != doubleparams.end(); ++iter)
    {
      nodeparams->set_double_param(iter->first,iter->second);
    }
  for (map<const string,int>::const_iterator iter = intparams.begin(); iter != intparams.end(); ++iter)
    {
      nodeparams->set_int_param(iter->first,iter->second);
    }
  for (map<const string,string>::const_iterator iter = stringparams.begin(); iter != stringparams.end(); ++iter)
    {
      nodeparams->set_string_param(iter->first,iter->second);
    }
  return;
}

void
PHG4Parameters::WriteToDB()
{
  PdbBankManager* bankManager = PdbBankManager::instance();
  PdbApplication *application = bankManager->getApplication();
  if (!application->startUpdate())
    {
      PHMessage("BunchCrossCal::", PHError, "Aborting ... Database not writable");
      application->abort();
    }

  //  Make a bank ID...
  PdbBankID bankID(0); // lets start at zero
  PHTimeStamp TStart(0);
  PHTimeStamp TStop(0xffffffff);
  PdbCalBank *NewBank = bankManager->createBank("PdbParameterMapBank",
						bankID,
						"hcaltest",
						TStart, TStop,
						"geo");
  if (NewBank)
    {
      NewBank->setLength(1);
	  PdbParameterMap *myparm = (PdbParameterMap*)&NewBank->getEntry(0);
	  for (map<const string, double>::const_iterator iter = doubleparams.begin(); iter != doubleparams.end(); ++iter)
	    {
	      myparm->set_double_param(iter->first, iter->second);
	    }
	  for (map<const string, int>::const_iterator iter = intparams.begin(); iter != intparams.end(); ++iter)
	    {
	      myparm->set_int_param(iter->first, iter->second);
	    }
	  for (map<const string, string>::const_iterator iter = stringparams.begin(); iter != stringparams.end(); ++iter)
	    {
	      myparm->set_string_param(iter->first, iter->second);
	    }
	//   myparm->setName(iter->first);
	//   myparm->setParameter(iter->second);
	// }
      application->commit(NewBank);
      delete NewBank;
    }
  return;
}
