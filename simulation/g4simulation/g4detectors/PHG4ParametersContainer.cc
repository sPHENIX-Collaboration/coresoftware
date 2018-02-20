#include "PHG4ParametersContainer.h"
#include "PHG4Parameters.h"

#include <pdbcalbase/PdbBankManager.h>
#include <pdbcalbase/PdbApplication.h>
#include <pdbcalbase/PdbBankList.h>
#include <pdbcalbase/PdbCalBank.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <pdbcalbase/PdbParameterMapContainer.h>

#include <phool/phool.h>
#include <phool/getClass.h>


#include <TBufferXML.h>
#include <TFile.h>
#include <TSystem.h>

#include <algorithm>
#include <iostream>
#include <sstream>

using namespace std;

PHG4ParametersContainer::PHG4ParametersContainer(const string &name):
  superdetectorname(name)
{}

PHG4ParametersContainer::~PHG4ParametersContainer()
{
  while(parametermap.begin() != parametermap.end())
    {
      delete parametermap.begin()->second;
      parametermap.erase(parametermap.begin());
    }

}

void
PHG4ParametersContainer::FillFrom(const PdbParameterMapContainer *saveparamcontainer)
{
// this fill only existing detids - no new ones are created (if the PdbParameterMapContainer contains
// entries from another detector)
  PdbParameterMapContainer::parConstRange begin_end =  saveparamcontainer->get_ParameterMaps();
  for (PdbParameterMapContainer::parIter iter = begin_end.first; iter != begin_end.second; ++iter)
  {
    Iterator pariter = parametermap.find(iter->first);
    if (pariter != parametermap.end())
    {
      PHG4Parameters *params = pariter->second;
      params->FillFrom(iter->second);
    }
  }
  return;
}

void
PHG4ParametersContainer::CreateAndFillFrom(const PdbParameterMapContainer *saveparamcontainer, const string &name)
{
  PdbParameterMapContainer::parConstRange begin_end =  saveparamcontainer->get_ParameterMaps();
  for (PdbParameterMapContainer::parIter iter = begin_end.first; iter != begin_end.second; ++iter)
  {
    Iterator pariter = parametermap.find(iter->first);
    if (pariter != parametermap.end())
    {
      PHG4Parameters *params = pariter->second;
      params->FillFrom(iter->second);
    }
    else
    {
      PHG4Parameters *params = new PHG4Parameters(name);
      params->FillFrom(iter->second);
      AddPHG4Parameters(iter->first,params);
    }
  }
  return;
}

void
PHG4ParametersContainer::AddPHG4Parameters(const int layer, PHG4Parameters *params)
{
  if (parametermap.find(layer) != parametermap.end())
    {
      cout << PHWHERE << " layer " << layer << " already exists for " 
	   << (parametermap.find(layer))->second->Name() << endl;
      gSystem->Exit(1);
    }
  parametermap[layer] = params;
}

const PHG4Parameters *
PHG4ParametersContainer::GetParameters(const int layer) const
{
  map<int, PHG4Parameters *>::const_iterator iter = parametermap.find(layer);
if (iter == parametermap.end())
  {
    cout << "could not find parameters for layer " << layer
	 << endl;
    return NULL;
  }
 return iter->second;
}

PHG4Parameters *
PHG4ParametersContainer::GetParametersToModify(const int layer)
{
  map<int, PHG4Parameters *>::iterator iter = parametermap.find(layer);
  if (iter == parametermap.end())
    {
      return NULL;
    }
  return iter->second;
}

int
PHG4ParametersContainer::WriteToFile(const string &extension, const string &dir)
{
  ostringstream fullpath;
  ostringstream fnamestream;
  PdbBankID bankID(0); // lets start at zero
  PHTimeStamp TStart(0);
  PHTimeStamp TStop(0xffffffff);
  fullpath << dir;
  // add / if directory lacks ending /
  if (*(dir.rbegin()) != '/')
    {
      fullpath << "/";
    }
  fnamestream << superdetectorname << "_geoparams" << "-" << bankID.getInternalValue()
      << "-" << TStart.getTics() << "-" << TStop.getTics() << "-" << time(0)
      << "." << extension;
  string fname = fnamestream.str();
  std::transform(fname.begin(), fname.end(), fname.begin(), ::tolower);
  fullpath << fname;

  cout <<"PHG4Parameters::WriteToFile - save to "<<fullpath.str()<<endl;

  PdbParameterMapContainer *myparm = new PdbParameterMapContainer();
  CopyToPdbParameterMapContainer(myparm);
  TFile *f = TFile::Open(fullpath.str().c_str(), "recreate");
  // force xml file writing to use extended precision shown experimentally
  // to not modify input parameters (.15e)
  string floatformat = TBufferXML::GetFloatFormat();
  TBufferXML::SetFloatFormat("%.15e");
  myparm->Write();
  delete f;
  // restore previous xml float format
  TBufferXML::SetFloatFormat(floatformat.c_str());
  cout << "sleeping 1 second to prevent duplicate inserttimes" << endl;
  sleep(1);
  return 0;
}

int
PHG4ParametersContainer::WriteToDB()
{
  PdbBankManager* bankManager = PdbBankManager::instance();
  PdbApplication *application = bankManager->getApplication();
  if (!application->startUpdate())
    {
      cout << PHWHERE << " Aborting, Database not writable" << endl;
      application->abort();
      exit(1);
    }

  //  Make a bank ID...
  PdbBankID bankID(0); // lets start at zero
  PHTimeStamp TStart(0);
  PHTimeStamp TStop(0xffffffff);

  string tablename = superdetectorname + "_geoparams";
  std::transform(tablename.begin(), tablename.end(), tablename.begin(),
      ::tolower);
  PdbCalBank *NewBank = bankManager->createBank("PdbParameterMapContainerBank", bankID,
      "Geometry Parameters", TStart, TStop, tablename);
  if (NewBank)
    {
      NewBank->setLength(1);
      PdbParameterMapContainer *myparm = (PdbParameterMapContainer*) &NewBank->getEntry(0);
      CopyToPdbParameterMapContainer(myparm);
      application->commit(NewBank);
      delete NewBank;
    }
  else
    {
      cout << PHWHERE " Committing to DB failed" << endl;
      return -1;
    }
  return 0;
}

void
PHG4ParametersContainer::CopyToPdbParameterMapContainer(PdbParameterMapContainer *myparmap)
{
  std::map<int, PHG4Parameters *>::const_iterator iter;
  for (iter = parametermap.begin(); iter != parametermap.end(); ++iter)
    {
      PdbParameterMap *myparm = new PdbParameterMap();
      iter->second->CopyToPdbParameterMap(myparm);
      myparmap->AddPdbParameterMap(iter->first,myparm);
    }
  return;
}

void
PHG4ParametersContainer::Print() const
{
  cout << "Name: " << Name() << endl;
  map<int, PHG4Parameters *>::const_iterator iter;
  for (iter = parametermap.begin(); iter != parametermap.end(); ++iter)
    {
      cout << "parameter detid: " << iter->first << endl;
      iter->second->Print();
    }
  return;
}

void
PHG4ParametersContainer::SaveToNodeTree(PHCompositeNode *topNode, const string &nodename)
{
  PdbParameterMapContainer *myparmap = findNode::getClass<PdbParameterMapContainer>(topNode, nodename);
  if (! myparmap)
    {
      myparmap = new PdbParameterMapContainer();
      PHIODataNode<PdbParameterMapContainer> *newnode =
          new PHIODataNode<PdbParameterMapContainer>(myparmap, nodename);
      topNode->addNode(newnode);
    }
  else
    {
       myparmap->Reset();
    }
  CopyToPdbParameterMapContainer(myparmap);
  return;
}

int
PHG4ParametersContainer::ExistDetid(const int detid) const
{
  if (parametermap.find(detid) != parametermap.end())
    {
      return 1;
    }
  return 0;
}
