#include "PHParameters.h"

#include <pdbcalbase/PdbApplication.h>
#include <pdbcalbase/PdbBankID.h>
#include <pdbcalbase/PdbBankManager.h>
#include <pdbcalbase/PdbCalBank.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <pdbcalbase/PdbParameterMapContainer.h>

#include <ffamodules/XploadInterface.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHTimeStamp.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TBufferXML.h>
#include <TFile.h>
#include <TSystem.h>

#include <boost/foreach.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/functional/hash.hpp>
#pragma GCC diagnostic pop
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
// stacktrace gives a shadow warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/stacktrace.hpp>
#pragma GCC diagnostic pop

#include <unistd.h>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <iterator>  // for reverse_iterator
#include <sstream>

PHParameters::PHParameters(const PHParameters &params, const std::string &name)
  : m_Detector(name)
{
  FillFrom(&params);
}

PHParameters::~PHParameters()
{
  m_DoubleParMap.clear();
  m_IntParMap.clear();
  m_StringParMap.clear();
}

void PHParameters::Reset()
{
  m_DoubleParMap.clear();
  m_IntParMap.clear();
  m_StringParMap.clear();
}

void PHParameters::set_int_param(const std::string &name, const int ival)
{
  m_IntParMap[name] = ival;
}

int PHParameters::get_int_param(const std::string &name) const
{
  if (m_IntParMap.find(name) != m_IntParMap.end())
  {
    return m_IntParMap.find(name)->second;
  }
  std::cout << PHWHERE << " integer parameter " << name
            << " does not exist (forgot to set?)" << std::endl;
  std::cout << "Here is the stacktrace: " << std::endl;
  std::cout << boost::stacktrace::stacktrace();
  std::cout << std::endl
            << "DO NOT PANIC - this is not a segfault" << std::endl;
  std::cout << "Check the stacktrace for the guilty party (typically #2)" << std::endl;
  gSystem->Exit(1);
  exit(1);
}

bool PHParameters::exist_int_param(const std::string &name) const
{
  if (m_IntParMap.find(name) != m_IntParMap.end())
  {
    return true;
  }
  return false;
}

void PHParameters::printint() const
{
  std::cout << "int parameters: " << std::endl;
  for (const auto &iter : m_IntParMap)
  {
    std::cout << iter.first << ": " << iter.second << std::endl;
  }
  return;
}

void PHParameters::set_double_param(const std::string &name, const double dval)
{
  m_DoubleParMap[name] = dval;
}

double
PHParameters::get_double_param(const std::string &name) const
{
  if (m_DoubleParMap.find(name) != m_DoubleParMap.end())
  {
    return m_DoubleParMap.find(name)->second;
  }
  std::cout << PHWHERE << " double parameter " << name
            << " does not exist (forgot to set?)" << std::endl;
  std::cout << "Here is the stacktrace: " << std::endl;
  std::cout << boost::stacktrace::stacktrace();
  std::cout << std::endl
            << "DO NOT PANIC - this is not a segfault" << std::endl;
  std::cout << "Check the stacktrace for the guilty party (typically #2)" << std::endl;

  gSystem->Exit(1);
  exit(1);
}

bool PHParameters::exist_double_param(const std::string &name) const
{
  if (m_DoubleParMap.find(name) != m_DoubleParMap.end())
  {
    return true;
  }
  return false;
}

void PHParameters::Print(Option_t * /*option*/) const
{
  std::cout << "Parameters for " << m_Detector << std::endl;
  printint();
  printdouble();
  printstring();
  return;
}

size_t
PHParameters::get_hash() const
{
  size_t seed = 0;

  for (const auto &iter : m_DoubleParMap)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, iter.first);
    boost::hash_combine(seed, iter.second);
    //      std::cout << iter->first << ": " << iter->second <<" -> "<<seed<< std::endl;
  }

  for (const auto &iter : m_IntParMap)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, iter.first);
    boost::hash_combine(seed, iter.second);
    //      std::cout << iter->first << ": " << iter->second <<" -> "<<seed<< std::endl;
  }

  for (const auto &iter : m_StringParMap)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, iter.first);
    boost::hash_combine(seed, iter.second);
    //      std::cout << iter->first << ": " << iter->second <<" -> "<<seed<< std::endl;
  }

  return seed;
}

void PHParameters::printdouble() const
{
  std::cout << "double parameters: " << std::endl;
  for (const auto &iter : m_DoubleParMap)
  {
    std::cout << iter.first << ": " << iter.second << std::endl;
  }
  return;
}

void PHParameters::set_string_param(const std::string &name, const std::string &str)
{
  m_StringParMap[name] = str;
}

std::string
PHParameters::get_string_param(const std::string &name) const
{
  if (m_StringParMap.find(name) != m_StringParMap.end())
  {
    return m_StringParMap.find(name)->second;
  }
  std::cout << PHWHERE << " string parameter " << name
            << " does not exist (forgot to set?)" << std::endl;
  std::cout << "Here is the stacktrace: " << std::endl;
  std::cout << boost::stacktrace::stacktrace();
  std::cout << std::endl
            << "DO NOT PANIC - this is not a segfault" << std::endl;
  std::cout << "Check the stacktrace for the guilty party (typically #2)" << std::endl;
  gSystem->Exit(1);
  exit(1);
}

bool PHParameters::exist_string_param(const std::string &name) const
{
  if (m_StringParMap.find(name) != m_StringParMap.end())
  {
    return true;
  }
  return false;
}

void PHParameters::printstring() const
{
  std::cout << "string parameters: " << std::endl;
  for (const auto &iter : m_StringParMap)
  {
    std::cout << iter.first << ": " << iter.second << std::endl;
  }
  return;
}

void PHParameters::FillFrom(const PdbParameterMap *saveparams)
{
  assert(saveparams);

  std::pair<std::map<const std::string, double>::const_iterator,
            std::map<const std::string, double>::const_iterator>
      begin_end_d = saveparams->get_dparam_iters();
  for (std::map<const std::string, double>::const_iterator iter = begin_end_d.first;
       iter != begin_end_d.second; ++iter)
  {
    m_DoubleParMap[iter->first] = iter->second;
  }
  std::pair<std::map<const std::string, int>::const_iterator,
            std::map<const std::string, int>::const_iterator>
      begin_end_i = saveparams->get_iparam_iters();
  for (std::map<const std::string, int>::const_iterator iter = begin_end_i.first;
       iter != begin_end_i.second; ++iter)
  {
    m_IntParMap[iter->first] = iter->second;
  }
  std::pair<std::map<const std::string, std::string>::const_iterator,
            std::map<const std::string, std::string>::const_iterator>
      begin_end_s = saveparams->get_cparam_iters();
  for (std::map<const std::string, std::string>::const_iterator iter = begin_end_s.first;
       iter != begin_end_s.second; ++iter)
  {
    m_StringParMap[iter->first] = iter->second;
  }

  return;
}

void PHParameters::FillFrom(const PdbParameterMapContainer *saveparamcontainer, const int detid)
{
  //  assert(saveparamcontainer != nullptr);

  const PdbParameterMap *saveparams = saveparamcontainer->GetParameters(detid);
  if (!saveparams)
  {
    return;
  }
  std::pair<std::map<const std::string, double>::const_iterator,
            std::map<const std::string, double>::const_iterator>
      begin_end_d = saveparams->get_dparam_iters();
  for (std::map<const std::string, double>::const_iterator iter = begin_end_d.first;
       iter != begin_end_d.second; ++iter)
  {
    m_DoubleParMap[iter->first] = iter->second;
  }
  std::pair<std::map<const std::string, int>::const_iterator,
            std::map<const std::string, int>::const_iterator>
      begin_end_i = saveparams->get_iparam_iters();
  for (std::map<const std::string, int>::const_iterator iter = begin_end_i.first;
       iter != begin_end_i.second; ++iter)
  {
    m_IntParMap[iter->first] = iter->second;
  }
  std::pair<std::map<const std::string, std::string>::const_iterator,
            std::map<const std::string, std::string>::const_iterator>
      begin_end_s = saveparams->get_cparam_iters();
  for (std::map<const std::string, std::string>::const_iterator iter = begin_end_s.first;
       iter != begin_end_s.second; ++iter)
  {
    m_StringParMap[iter->first] = iter->second;
  }

  return;
}

void PHParameters::FillFrom(const PHParameters *saveparams)
{
  assert(saveparams);

  for (const auto &iter : saveparams->m_DoubleParMap)
  {
    m_DoubleParMap[iter.first] = iter.second;
  }

  for (const auto &iter : saveparams->m_IntParMap)
  {
    m_IntParMap[iter.first] = iter.second;
  }

  for (const auto &iter : saveparams->m_StringParMap)
  {
    m_StringParMap[iter.first] = iter.second;
  }
  return;
}

void PHParameters::SaveToNodeTree(PHCompositeNode *topNode, const std::string &nodename)
{
  // write itself since this class is fine with saving by root
  PdbParameterMap *nodeparams = findNode::getClass<PdbParameterMap>(topNode, nodename);
  if (!nodeparams)
  {
    nodeparams = new PdbParameterMap();
    PHIODataNode<PdbParameterMap> *newnode = new PHIODataNode<PdbParameterMap>(nodeparams, nodename);
    topNode->addNode(newnode);
  }
  else
  {
    nodeparams->Reset();  // just clear previous content in case variables were deleted
  }
  CopyToPdbParameterMap(nodeparams);
  return;
}

void PHParameters::UpdateNodeTree(PHCompositeNode *topNode, const std::string &nodename)
{
  PdbParameterMap *nodeparams = findNode::getClass<PdbParameterMap>(topNode, nodename);
  if (!nodeparams)
  {
    std::cout << PHWHERE << " could not find PdbParameterMap " << nodename
              << " which must exist" << std::endl;
    gSystem->Exit(1);
  }
  CopyToPdbParameterMap(nodeparams);
  return;
}

void PHParameters::SaveToNodeTree(PHCompositeNode *topNode, const std::string &nodename, const int detid)
{
  // write itself since this class is fine with saving by root
  PdbParameterMapContainer *nodeparamcontainer = findNode::getClass<PdbParameterMapContainer>(topNode, nodename);
  if (!nodeparamcontainer)
  {
    nodeparamcontainer = new PdbParameterMapContainer();
    PHIODataNode<PdbParameterMapContainer> *newnode =
        new PHIODataNode<PdbParameterMapContainer>(nodeparamcontainer, nodename);
    topNode->addNode(newnode);
  }
  PdbParameterMap *nodeparams = nodeparamcontainer->GetParametersToModify(detid);
  if (nodeparams)
  {
    nodeparams->Reset();
  }
  else
  {
    nodeparams = new PdbParameterMap();
    nodeparamcontainer->AddPdbParameterMap(detid, nodeparams);
  }
  CopyToPdbParameterMap(nodeparams);
  return;
}

void PHParameters::UpdateNodeTree(PHCompositeNode *topNode, const std::string &nodename, const int detid)
{
  PdbParameterMapContainer *nodeparamcontainer = findNode::getClass<PdbParameterMapContainer>(topNode, nodename);
  if (!nodeparamcontainer)
  {
    std::cout << PHWHERE << " could not find PdbParameterMapContainer " << nodename
              << " which must exist" << std::endl;
    gSystem->Exit(1);
  }
  PdbParameterMap *nodeparams = nodeparamcontainer->GetParametersToModify(detid);
  if (!nodeparams)
  {
    std::cout << PHWHERE << " could not find PdbParameterMap for detector " << detid
              << " which must exist" << std::endl;
    gSystem->Exit(1);
  }
  CopyToPdbParameterMap(nodeparams);
  return;
}

int PHParameters::WriteToDB()
{
  PdbBankManager *bankManager = PdbBankManager::instance();
  PdbApplication *application = bankManager->getApplication();
  if (!application->startUpdate())
  {
    std::cout << PHWHERE << " Aborting, Database not writable" << std::endl;
    application->abort();
    gSystem->Exit(1);
    exit(1);
  }

  //  Make a bank ID...
  PdbBankID bankID(0);  // lets start at zero
  PHTimeStamp TStart(0);
  PHTimeStamp TStop(0xffffffff);

  std::string tablename = m_Detector + "_geoparams";
  std::transform(tablename.begin(), tablename.end(), tablename.begin(), ::tolower);
  PdbCalBank *NewBank = bankManager->createBank("PdbParameterMapBank", bankID,
                                                "Geometry Parameters", TStart, TStop, tablename);
  if (NewBank)
  {
    NewBank->setLength(1);
    PdbParameterMap *myparm = (PdbParameterMap *) &NewBank->getEntry(0);
    CopyToPdbParameterMap(myparm);
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

int PHParameters::ReadFromDB(const std::string &name, const int detid)
{
  PdbBankManager *bankManager = PdbBankManager::instance();
  PdbApplication *application = bankManager->getApplication();
  if (!application->startRead())
  {
    std::cout << PHWHERE << " Aborting, Database not readable" << std::endl;
    application->abort();
    gSystem->Exit(1);
    exit(1);
  }

  //  Make a bank ID...
  PdbBankID bankID(0);  // lets start at zero
  PHTimeStamp TSearch(10);

  std::string tablename = name + "_geoparams";
  std::transform(tablename.begin(), tablename.end(), tablename.begin(), ::tolower);
  PdbCalBank *NewBank = bankManager->fetchBank("PdbParameterMapContainerBank", bankID, tablename, TSearch);
  if (NewBank)
  {
    PdbParameterMapContainer *myparm = (PdbParameterMapContainer *) &NewBank->getEntry(0);
    FillFrom(myparm, detid);
    delete NewBank;
  }
  else
  {
    std::cout << PHWHERE " Reading from DB failed" << std::endl;
    return -1;
  }
  return 0;
}

int PHParameters::ReadFromDB()
{
  PdbBankManager *bankManager = PdbBankManager::instance();
  PdbApplication *application = bankManager->getApplication();
  if (!application->startRead())
  {
    std::cout << PHWHERE << " Aborting, Database not readable" << std::endl;
    application->abort();
    gSystem->Exit(1);
    exit(1);
  }

  //  Make a bank ID...
  PdbBankID bankID(0);  // lets start at zero
  PHTimeStamp TSearch(10);

  std::string tablename = m_Detector + "_geoparams";
  std::transform(tablename.begin(), tablename.end(), tablename.begin(), ::tolower);
  PdbCalBank *NewBank = bankManager->fetchBank("PdbParameterMapBank", bankID,
                                               tablename, TSearch);
  if (NewBank)
  {
    PdbParameterMap *myparm = (PdbParameterMap *) &NewBank->getEntry(0);
    FillFrom(myparm);
    delete NewBank;
  }
  else
  {
    std::cout << PHWHERE " Reading from DB failed" << std::endl;
    return -1;
  }
  return 0;
}

int PHParameters::ReadFromCDB(const std::string &domain)
{
  std::string url = XploadInterface::instance()->getUrl(domain);
  TFile *f = TFile::Open(url.c_str());
  if (!f)
  {
    std::cout << "could not open " << url << std::endl;
    gSystem->Exit(1);
  }
  PdbParameterMap *myparm = static_cast<PdbParameterMap *>(f->Get("PdbParameterMap"));
  if (!myparm)
  {
    std::cout << "could not get PdbParameterMap from " << url << std::endl;
    gSystem->Exit(1);
  }
  FillFrom(myparm);
  delete myparm;
  delete f;
  return 0;
}

int PHParameters::WriteToCDBFile(const std::string &filename)
{
  PdbParameterMap *myparm = new PdbParameterMap();
  CopyToPdbParameterMap(myparm);
  TFile *f = TFile::Open(filename.c_str(), "recreate");
  myparm->Write();
  delete f;
  delete myparm;
  return 0;
}

int PHParameters::WriteToFile(const std::string &extension, const std::string &dir)
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
  fnamestream << m_Detector << "_geoparams"
              << "-" << bankID.getInternalValue()
              << "-" << TStart.getTics() << "-" << TStop.getTics() << "-" << time(nullptr)
              << "." << extension;
  std::string fname = fnamestream.str();
  std::transform(fname.begin(), fname.end(), fname.begin(), ::tolower);
  fullpath << fname;

  std::cout << "PHParameters::WriteToFile - save to " << fullpath.str() << std::endl;

  PdbParameterMap *myparm = new PdbParameterMap();
  CopyToPdbParameterMap(myparm);
  TFile *f = TFile::Open(fullpath.str().c_str(), "recreate");
  // force xml file writing to use extended precision shown experimentally
  // to not modify input parameters (.17g)
  std::string floatformat = TBufferXML::GetFloatFormat();
  TBufferXML::SetFloatFormat("%.17g");  // for IEEE 754 double
  myparm->Write();
  delete f;
  // restore previous xml float format
  TBufferXML::SetFloatFormat(floatformat.c_str());
  std::cout << "sleeping 1 second to prevent duplicate inserttimes" << std::endl;
  sleep(1);
  return 0;
}

int PHParameters::ReadFromFile(const std::string &name, const std::string &extension, const int detid, const int issuper, const std::string &dir)
{
  PHTimeStamp TSearch(10);
  PdbBankID bankID(0);
  std::ostringstream fnamestream;
  fnamestream << name << "_geoparams"
              << "-" << bankID.getInternalValue();
  std::string fileprefix = fnamestream.str();
  std::transform(fileprefix.begin(), fileprefix.end(), fileprefix.begin(),
                 ::tolower);
  std::filesystem::path targetDir(dir);

  std::filesystem::recursive_directory_iterator diriter(targetDir), eod;
  boost::char_separator<char> sep("-.");
  std::map<unsigned int, std::string> calibfiles;
  BOOST_FOREACH (std::filesystem::path const &i, std::make_pair(diriter, eod))
  {
    if (is_regular_file(i))
    {
      // leaf() gives the filename without path,
      // the string.compare(0...) checks if the filename starts with fileprefix
      // if not coninue
      std::string basename = i.filename().string();
      if (basename.compare(0, fileprefix.size(), fileprefix) != 0)
      {
        continue;
      }
      // extension() contains the . - like .xml, so we
      // just compare the extensions instead of !=
      // and check that the size makes sense
      if (i.extension().string().find(extension) == std::string::npos || i.extension().string().size() != extension.size() + 1)
      {
        continue;
      }
      boost::tokenizer<boost::char_separator<char> > tok(basename, sep);
      boost::tokenizer<boost::char_separator<char> >::iterator iter =
          tok.begin();
      ++iter;  // that skips the file prefix excluding bank id
      ++iter;  // that skips the bank id we checked already as part of the filename
      PHTimeStamp TStart(ConvertStringToUint(*iter));
      if (TSearch < TStart)
      {
        continue;
      }
      ++iter;
      PHTimeStamp TStop(ConvertStringToUint(*iter));
      if (TSearch >= TStop)
      {
        continue;
      }
      ++iter;
      calibfiles[ConvertStringToUint(*iter)] = i.string();
    }
  }
  if (calibfiles.empty())
  {
    std::cout << "No calibration file like " << dir << "/" << fileprefix << " found" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  std::cout << "PHParameters::ReadFromFile - Reading from File: " << (calibfiles.rbegin())->second << " ... ";
  std::string fname = (calibfiles.rbegin())->second;
  TFile *f = TFile::Open(fname.c_str());
  if (issuper)
  {
    PdbParameterMapContainer *myparm = static_cast<PdbParameterMapContainer *>(f->Get("PdbParameterMapContainer"));
    assert(myparm);

    if (myparm->GetParameters(detid) == nullptr)
    {
      std::cout << "Missing PdbParameterMapContainer Detector Id " << detid << std::endl;
      gSystem->Exit(1);
    }

    std::cout << "Received PdbParameterMapContainer Detector Id " << detid << " with (Hash = 0x" << std::hex << myparm->GetParameters(detid)->get_hash() << std::dec << ")" << std::endl;

    FillFrom(myparm, detid);
    delete myparm;
  }
  else
  {
    PdbParameterMap *myparm = static_cast<PdbParameterMap *>(f->Get("PdbParameterMap"));
    assert(myparm);
    std::cout << "Received PdbParameterMap with (Hash = 0x" << std::hex << myparm->get_hash() << std::dec << ")" << std::endl;

    FillFrom(myparm);
    delete myparm;
  }
  delete f;

  return 0;
}

int PHParameters::ReadFromCDBFile(const std::string &url)
{
  TFile *f = TFile::Open(url.c_str());
  PdbParameterMap *myparm = static_cast<PdbParameterMap *>(f->Get("PdbParameterMap"));
  FillFrom(myparm);
  delete myparm;
  delete f;
  return 0;
}

void PHParameters::CopyToPdbParameterMap(PdbParameterMap *myparm)
{
  for (std::map<const std::string, double>::const_iterator iter = m_DoubleParMap.begin();
       iter != m_DoubleParMap.end(); ++iter)
  {
    myparm->set_double_param(iter->first, iter->second);
  }
  for (std::map<const std::string, int>::const_iterator iter = m_IntParMap.begin();
       iter != m_IntParMap.end(); ++iter)
  {
    myparm->set_int_param(iter->first, iter->second);
  }
  for (std::map<const std::string, std::string>::const_iterator iter = m_StringParMap.begin();
       iter != m_StringParMap.end(); ++iter)
  {
    myparm->set_string_param(iter->first, iter->second);
  }
}

unsigned int
PHParameters::ConvertStringToUint(const std::string &str) const
{
  unsigned int tics;
  try
  {
    tics = boost::lexical_cast<unsigned int>(str);
  }
  catch (boost::bad_lexical_cast const &)
  {
    std::cout << "Cannot extract timestamp from " << str << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  return tics;
}
