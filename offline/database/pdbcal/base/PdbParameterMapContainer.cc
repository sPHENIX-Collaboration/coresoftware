#include "PdbParameterMapContainer.h"
#include "PdbBankID.h"
#include "PdbParameterMap.h"

#include <phool/PHTimeStamp.h>
#include <phool/phool.h>

#include <TBufferXML.h>
#include <TFile.h>
#include <TSystem.h>

#include <boost/stacktrace.hpp>

#include <unistd.h>
#include <algorithm>
#include <cctype>
#include <ctime>
#include <iostream>
#include <sstream>

PdbParameterMapContainer::~PdbParameterMapContainer()
{
  while (parametermap.begin() != parametermap.end())
  {
    delete parametermap.begin()->second;
    parametermap.erase(parametermap.begin());
  }
  return;
}

void PdbParameterMapContainer::print() const
{
  for (auto iter : parametermap)
  {
    std::cout << "layer " << iter.first << std::endl;
    iter.second->print();
  }
  return;
}
void PdbParameterMapContainer::Reset()
{
  while (parametermap.begin() != parametermap.end())
  {
    delete parametermap.begin()->second;
    parametermap.erase(parametermap.begin());
  }
  return;
}

void PdbParameterMapContainer::AddPdbParameterMap(const int layer, PdbParameterMap *params)
{
  if (parametermap.contains(layer))
  {
    std::cout << PHWHERE << " layer " << layer << " already exists" << std::endl;
    std::cout << "Here is the stacktrace: " << std::endl;
    std::cout << boost::stacktrace::stacktrace();
    std::cout << std::endl
              << "DO NOT PANIC - this is not a segfault" << std::endl;
    std::cout << "Check the stacktrace for the guilty party (typically #2)" << std::endl;
    gSystem->Exit(1);
  }
  parametermap[layer] = params;
}

const PdbParameterMap *
PdbParameterMapContainer::GetParameters(const int layer) const
{
  std::map<int, PdbParameterMap *>::const_iterator iter = parametermap.find(layer);
  if (iter == parametermap.end())
  {
    return nullptr;
  }
  return iter->second;
}

PdbParameterMap *
PdbParameterMapContainer::GetParametersToModify(const int layer)
{
  std::map<int, PdbParameterMap *>::iterator iter = parametermap.find(layer);
  if (iter == parametermap.end())
  {
    return nullptr;
  }
  return iter->second;
}

int PdbParameterMapContainer::WriteToFile(const std::string &detector_name,
                                          const std::string &extension, const std::string &dir)
{
  // Note the naming convention should be consistent with PHParameters::WriteToFile

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
  fnamestream << detector_name << "_geoparams"
              << "-"
              << bankID.getInternalValue() << "-" << TStart.getTics() << "-"
              << TStop.getTics() << "-" << time(nullptr) << "." << extension;
  std::string fname = fnamestream.str();
  std::transform(fname.begin(), fname.end(), fname.begin(), ::tolower);
  fullpath << fname;

  std::cout << "PdbParameterMapContainer::WriteToFile - save to " << fullpath.str()
            << std::endl;

  TFile *f = TFile::Open(fullpath.str().c_str(), "recreate");

  PdbParameterMapContainer *container = new PdbParameterMapContainer();
  for (std::map<int, PdbParameterMap *>::const_iterator it =
           parametermap.begin();
       it != parametermap.end(); ++it)
  {
    PdbParameterMap *myparm = dynamic_cast<PdbParameterMap *>(it->second->CloneMe());
    container->AddPdbParameterMap(it->first, myparm);
  }

  // force xml file writing to use extended precision shown experimentally
  // to not modify input parameters (.15e)
  std::string floatformat = TBufferXML::GetFloatFormat();
  TBufferXML::SetFloatFormat("%.17g");  // for IEEE 754 double
  container->Write("PdbParameterMapContainer");
  delete f;
  // restore previous xml float format
  TBufferXML::SetFloatFormat(floatformat.c_str());
  std::cout << "sleeping 1 second to prevent duplicate inserttimes" << std::endl;
  sleep(1);
  return 0;
}
