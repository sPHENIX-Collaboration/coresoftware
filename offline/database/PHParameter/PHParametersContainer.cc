#include "PHParametersContainer.h"
#include "PHParameters.h"

#include <pdbcalbase/PdbApplication.h>
#include <pdbcalbase/PdbBankID.h>
#include <pdbcalbase/PdbBankManager.h>
#include <pdbcalbase/PdbCalBank.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <pdbcalbase/PdbParameterMapContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHTimeStamp.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TBufferXML.h>
#include <TFile.h>
#include <TSystem.h>

// stacktrace gives a shadow warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/stacktrace.hpp>
#pragma GCC diagnostic pop

#include <unistd.h>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>

PHParametersContainer::PHParametersContainer(const std::string &name)
  : superdetectorname(name)
{
}

PHParametersContainer::~PHParametersContainer()
{
  while (parametermap.begin() != parametermap.end())
  {
    delete parametermap.begin()->second;
    parametermap.erase(parametermap.begin());
  }
}

void PHParametersContainer::FillFrom(const PdbParameterMapContainer *saveparamcontainer)
{
  // this fill only existing detids - no new ones are created (if the PdbParameterMapContainer contains
  // entries from another detector)
  PdbParameterMapContainer::parConstRange begin_end = saveparamcontainer->get_ParameterMaps();
  for (PdbParameterMapContainer::parIter iter = begin_end.first; iter != begin_end.second; ++iter)
  {
    Iterator pariter = parametermap.find(iter->first);
    if (pariter != parametermap.end())
    {
      PHParameters *params = pariter->second;
      params->FillFrom(iter->second);
    }
  }
  return;
}

void PHParametersContainer::CreateAndFillFrom(const PdbParameterMapContainer *saveparamcontainer, const std::string &name)
{
  PdbParameterMapContainer::parConstRange begin_end = saveparamcontainer->get_ParameterMaps();
  for (PdbParameterMapContainer::parIter iter = begin_end.first; iter != begin_end.second; ++iter)
  {
    Iterator pariter = parametermap.find(iter->first);
    if (pariter != parametermap.end())
    {
      PHParameters *params = pariter->second;
      params->FillFrom(iter->second);
    }
    else
    {
      PHParameters *params = new PHParameters(name);
      params->FillFrom(iter->second);
      AddPHParameters(iter->first, params);
    }
  }
  return;
}

void PHParametersContainer::AddPHParameters(const int detid, PHParameters *params)
{
  if (parametermap.find(detid) != parametermap.end())
  {
    std::cout << PHWHERE << " detector id " << detid << " already exists for "
              << (parametermap.find(detid))->second->Name() << std::endl;
    gSystem->Exit(1);
  }
  parametermap[detid] = params;
}

const PHParameters *
PHParametersContainer::GetParameters(const int detid) const
{
  std::map<int, PHParameters *>::const_iterator iter = parametermap.find(detid);
  if (iter == parametermap.end())
  {
    std::cout << "could not find parameters for detector id " << detid
              << std::endl;
    std::cout << "Here is the stacktrace: " << std::endl;
    std::cout << boost::stacktrace::stacktrace();
    std::cout << std::endl
              << "DO NOT PANIC - this is not a segfault" << std::endl;
    std::cout << "Check the stacktrace for the guilty party (typically #2)" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  return iter->second;
}

PHParameters *
PHParametersContainer::GetParametersToModify(const int detid)
{
  std::map<int, PHParameters *>::iterator iter = parametermap.find(detid);
  if (iter == parametermap.end())
  {
    return nullptr;
  }
  return iter->second;
}

int PHParametersContainer::WriteToFile(const std::string &extension, const std::string &dir)
{
  std::ostringstream fullpath;
  std::ostringstream fnamestream;
  PdbBankID bankID(0);  // lets start at zero
  PHTimeStamp TStart(0);
  PHTimeStamp TStop(0xffffffff);
  fullpath << dir;
  // add / if directory lacks ending /
  if (*(dir.rbegin()) != '/')
  {
    fullpath << "/";
  }
  fnamestream << superdetectorname << "_geoparams"
              << "-" << bankID.getInternalValue()
              << "-" << TStart.getTics() << "-" << TStop.getTics() << "-" << time(nullptr)
              << "." << extension;
  std::string fname = fnamestream.str();
  std::transform(fname.begin(), fname.end(), fname.begin(), ::tolower);
  fullpath << fname;

  std::cout << "PHParameters::WriteToFile - save to " << fullpath.str() << std::endl;

  PdbParameterMapContainer *myparm = new PdbParameterMapContainer();
  CopyToPdbParameterMapContainer(myparm);
  TFile *f = TFile::Open(fullpath.str().c_str(), "recreate");
  // force xml file writing to use extended precision shown experimentally
  // to not modify input parameters (.15e)
  std::string floatformat = TBufferXML::GetFloatFormat();
  TBufferXML::SetFloatFormat("%.15e");
  myparm->Write();
  delete f;
  // restore previous xml float format
  TBufferXML::SetFloatFormat(floatformat.c_str());
  std::cout << "sleeping 1 second to prevent duplicate inserttimes" << std::endl;
  sleep(1);
  return 0;
}

int PHParametersContainer::WriteToDB()
{
  PdbBankManager *bankManager = PdbBankManager::instance();
  PdbApplication *application = bankManager->getApplication();
  if (!application->startUpdate())
  {
    std::cout << PHWHERE << " Aborting, Database not writable" << std::endl;
    application->abort();
    exit(1);
  }

  //  Make a bank ID...
  PdbBankID bankID(0);  // lets start at zero
  PHTimeStamp TStart(0);
  PHTimeStamp TStop(0xffffffff);

  std::string tablename = superdetectorname + "_geoparams";
  std::transform(tablename.begin(), tablename.end(), tablename.begin(),
                 ::tolower);
  PdbCalBank *NewBank = bankManager->createBank("PdbParameterMapContainerBank", bankID,
                                                "Geometry Parameters", TStart, TStop, tablename);
  if (NewBank)
  {
    NewBank->setLength(1);
    PdbParameterMapContainer *myparm = (PdbParameterMapContainer *) &NewBank->getEntry(0);
    CopyToPdbParameterMapContainer(myparm);
    application->commit(NewBank);
    delete NewBank;
  }
  else
  {
    std::cout << PHWHERE " Committing to DB failed" << std::endl;
    return -1;
  }
  return 0;
}

void PHParametersContainer::CopyToPdbParameterMapContainer(PdbParameterMapContainer *myparmap)
{
  std::map<int, PHParameters *>::const_iterator iter;
  for (iter = parametermap.begin(); iter != parametermap.end(); ++iter)
  {
    PdbParameterMap *myparm = new PdbParameterMap();
    iter->second->CopyToPdbParameterMap(myparm);
    myparmap->AddPdbParameterMap(iter->first, myparm);
  }
  return;
}

void PHParametersContainer::UpdatePdbParameterMapContainer(PdbParameterMapContainer *myparmap)
{
  std::map<int, PHParameters *>::const_iterator iter;
  for (iter = parametermap.begin(); iter != parametermap.end(); ++iter)
  {
    PdbParameterMap *nodeparams = myparmap->GetParametersToModify(iter->first);
    iter->second->CopyToPdbParameterMap(nodeparams);
  }
  return;
}

void PHParametersContainer::Print(Option_t * /*option*/) const
{
  std::cout << "Name: " << Name() << std::endl;
  std::map<int, PHParameters *>::const_iterator iter;
  for (iter = parametermap.begin(); iter != parametermap.end(); ++iter)
  {
    std::cout << "parameter detid: " << iter->first << std::endl;
    iter->second->Print();
  }
  return;
}

void PHParametersContainer::SaveToNodeTree(PHCompositeNode *topNode, const std::string &nodename)
{
  PdbParameterMapContainer *myparmap = findNode::getClass<PdbParameterMapContainer>(topNode, nodename);
  if (!myparmap)
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

void PHParametersContainer::UpdateNodeTree(PHCompositeNode *topNode, const std::string &nodename)
{
  PdbParameterMapContainer *myparmap = findNode::getClass<PdbParameterMapContainer>(topNode, nodename);
  if (!myparmap)
  {
    std::cout << PHWHERE << " could not find PdbParameterMapContainer " << nodename
              << " which must exist" << std::endl;
    gSystem->Exit(1);
  }
  UpdatePdbParameterMapContainer(myparmap);
  return;
}

int PHParametersContainer::ExistDetid(const int detid) const
{
  if (parametermap.find(detid) != parametermap.end())
  {
    return 1;
  }
  return 0;
}
