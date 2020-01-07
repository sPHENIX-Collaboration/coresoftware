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
//#include <TBufferFile.h>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>
#include <boost/stacktrace.hpp>
#include <boost/tokenizer.hpp>

// this is an ugly hack, the gcc optimizer has a bug which
// triggers the uninitialized variable warning which
// stops compilation because of our -Werror
#include <boost/version.hpp>  // to get BOOST_VERSION
#if (__GNUC__ == 4 && __GNUC_MINOR__ == 4 && BOOST_VERSION == 105700)
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma message "ignoring bogus gcc warning in boost header lexical_cast.hpp"
#include <boost/lexical_cast.hpp>
#pragma GCC diagnostic warning "-Wuninitialized"
#else
#include <boost/lexical_cast.hpp>
#endif

#include <unistd.h>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>

using namespace std;

PHParameters::PHParameters(const PHParameters &params, const std::string &name)
{
  set_name(name);
  FillFrom(&params);
}

PHParameters::~PHParameters()
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
  cout << PHWHERE << " integer parameter " << name
       << " does not exist (forgot to set?)" << endl;
  cout << "Here is the stacktrace: " << endl;
  cout << boost::stacktrace::stacktrace();
  cout << "Check the stacktrace for the guilty party (typically #2)" << endl;
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
  cout << "int parameters: " << endl;
  for (map<const string, int>::const_iterator iter = m_IntParMap.begin();
       iter != m_IntParMap.end(); ++iter)
  {
    cout << iter->first << ": " << iter->second << endl;
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
  cout << PHWHERE << " double parameter " << name
       << " does not exist (forgot to set?)" << endl;
  cout << "Here is the stacktrace: " << endl;
  cout << boost::stacktrace::stacktrace();
  cout << "Check the stacktrace for the guilty party (typically #2)" << endl;

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

void PHParameters::Print(Option_t *option) const
{
  cout << "Parameters for " << m_Detector << endl;
  printint();
  printdouble();
  printstring();
  return;
}

size_t
PHParameters::get_hash() const
{
  size_t seed = 0;

  for (dMap::const_iterator iter = m_DoubleParMap.begin();
       iter != m_DoubleParMap.end(); ++iter)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, iter->first);
    boost::hash_combine(seed, iter->second);
    //      cout << iter->first << ": " << iter->second <<" -> "<<seed<< endl;
  }

  for (iMap::const_iterator iter = m_IntParMap.begin();
       iter != m_IntParMap.end(); ++iter)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, iter->first);
    boost::hash_combine(seed, iter->second);
    //      cout << iter->first << ": " << iter->second <<" -> "<<seed<< endl;
  }

  for (strMap::const_iterator iter = m_StringParMap.begin();
       iter != m_StringParMap.end(); ++iter)
  {
    //      size_t seed = 0;
    boost::hash_combine(seed, iter->first);
    boost::hash_combine(seed, iter->second);
    //      cout << iter->first << ": " << iter->second <<" -> "<<seed<< endl;
  }

  return seed;
}

void PHParameters::printdouble() const
{
  cout << "double parameters: " << endl;
  for (map<const string, double>::const_iterator iter = m_DoubleParMap.begin();
       iter != m_DoubleParMap.end(); ++iter)
  {
    cout << iter->first << ": " << iter->second << endl;
  }
  return;
}

void PHParameters::set_string_param(const std::string &name, const string &str)
{
  m_StringParMap[name] = str;
}

string
PHParameters::get_string_param(const std::string &name) const
{
  if (m_StringParMap.find(name) != m_StringParMap.end())
  {
    return m_StringParMap.find(name)->second;
  }
  cout << PHWHERE << " string parameter " << name
       << " does not exist (forgot to set?)" << endl;
  cout << "Here is the stacktrace: " << endl;
  cout << boost::stacktrace::stacktrace();
  cout << "Check the stacktrace for the guilty party (typically #2)" << endl;
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
  cout << "string parameters: " << endl;
  for (map<const string, string>::const_iterator iter = m_StringParMap.begin();
       iter != m_StringParMap.end(); ++iter)
  {
    cout << iter->first << ": " << iter->second << endl;
  }
  return;
}

void PHParameters::FillFrom(const PdbParameterMap *saveparams)
{
  assert(saveparams);

  pair<std::map<const std::string, double>::const_iterator,
       std::map<const std::string, double>::const_iterator>
      begin_end_d =
          saveparams->get_dparam_iters();
  for (map<const std::string, double>::const_iterator iter = begin_end_d.first;
       iter != begin_end_d.second; ++iter)
  {
    m_DoubleParMap[iter->first] = iter->second;
  }
  pair<std::map<const std::string, int>::const_iterator,
       std::map<const std::string, int>::const_iterator>
      begin_end_i =
          saveparams->get_iparam_iters();
  for (map<const std::string, int>::const_iterator iter = begin_end_i.first;
       iter != begin_end_i.second; ++iter)
  {
    m_IntParMap[iter->first] = iter->second;
  }
  pair<std::map<const std::string, string>::const_iterator,
       std::map<const std::string, string>::const_iterator>
      begin_end_s =
          saveparams->get_cparam_iters();
  for (map<const std::string, string>::const_iterator iter = begin_end_s.first;
       iter != begin_end_s.second; ++iter)
  {
    m_StringParMap[iter->first] = iter->second;
  }

  return;
}

void PHParameters::FillFrom(const PdbParameterMapContainer *saveparamcontainer, const int layer)
{
  //  assert(saveparamcontainer != NULL);

  const PdbParameterMap *saveparams = saveparamcontainer->GetParameters(layer);
  if (!saveparams)
  {
    return;
  }
  pair<std::map<const std::string, double>::const_iterator,
       std::map<const std::string, double>::const_iterator>
      begin_end_d =
          saveparams->get_dparam_iters();
  for (map<const std::string, double>::const_iterator iter = begin_end_d.first;
       iter != begin_end_d.second; ++iter)
  {
    m_DoubleParMap[iter->first] = iter->second;
  }
  pair<std::map<const std::string, int>::const_iterator,
       std::map<const std::string, int>::const_iterator>
      begin_end_i =
          saveparams->get_iparam_iters();
  for (map<const std::string, int>::const_iterator iter = begin_end_i.first;
       iter != begin_end_i.second; ++iter)
  {
    m_IntParMap[iter->first] = iter->second;
  }
  pair<std::map<const std::string, string>::const_iterator,
       std::map<const std::string, string>::const_iterator>
      begin_end_s =
          saveparams->get_cparam_iters();
  for (map<const std::string, string>::const_iterator iter = begin_end_s.first;
       iter != begin_end_s.second; ++iter)
  {
    m_StringParMap[iter->first] = iter->second;
  }

  return;
}

void PHParameters::FillFrom(const PHParameters *saveparams)
{
  assert(saveparams);

  for (dMap::const_iterator iter = saveparams->m_DoubleParMap.begin();
       iter != saveparams->m_DoubleParMap.end(); ++iter)
    m_DoubleParMap[iter->first] = iter->second;

  for (iMap::const_iterator iter = saveparams->m_IntParMap.begin();
       iter != saveparams->m_IntParMap.end(); ++iter)
    m_IntParMap[iter->first] = iter->second;

  for (strMap::const_iterator iter = saveparams->m_StringParMap.begin();
       iter != saveparams->m_StringParMap.end(); ++iter)
    m_StringParMap[iter->first] = iter->second;

  return;
}

void PHParameters::SaveToNodeTree(PHCompositeNode *topNode, const string &nodename)
{
  // write itself since this class is fine with saving by root
  PdbParameterMap *nodeparams = findNode::getClass<PdbParameterMap>(topNode,
                                                                    nodename);
  if (!nodeparams)
  {
    nodeparams = new PdbParameterMap();
    PHIODataNode<PdbParameterMap> *newnode =
        new PHIODataNode<PdbParameterMap>(nodeparams, nodename);
    topNode->addNode(newnode);
  }
  else
  {
    nodeparams->Reset();  // just clear previous content in case variables were deleted
  }
  CopyToPdbParameterMap(nodeparams);
  return;
}

void PHParameters::SaveToNodeTree(PHCompositeNode *topNode, const string &nodename, const int layer)
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
  PdbParameterMap *nodeparams = nodeparamcontainer->GetParametersToModify(layer);
  if (nodeparams)
  {
    nodeparams->Reset();
  }
  else
  {
    nodeparams = new PdbParameterMap();
    nodeparamcontainer->AddPdbParameterMap(layer, nodeparams);
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
    cout << PHWHERE << " Aborting, Database not writable" << endl;
    application->abort();
    gSystem->Exit(1);
    exit(1);
  }

  //  Make a bank ID...
  PdbBankID bankID(0);  // lets start at zero
  PHTimeStamp TStart(0);
  PHTimeStamp TStop(0xffffffff);

  string tablename = m_Detector + "_geoparams";
  std::transform(tablename.begin(), tablename.end(), tablename.begin(),
                 ::tolower);
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
    cout << PHWHERE " Committing to DB failed" << endl;
    return -1;
  }
  return 0;
}

int PHParameters::ReadFromDB(const string &name, const int layer)
{
  PdbBankManager *bankManager = PdbBankManager::instance();
  PdbApplication *application = bankManager->getApplication();
  if (!application->startRead())
  {
    cout << PHWHERE << " Aborting, Database not readable" << endl;
    application->abort();
    gSystem->Exit(1);
    exit(1);
  }

  //  Make a bank ID...
  PdbBankID bankID(0);  // lets start at zero
  PHTimeStamp TSearch(10);

  string tablename = name + "_geoparams";
  std::transform(tablename.begin(), tablename.end(), tablename.begin(),
                 ::tolower);
  PdbCalBank *NewBank = bankManager->fetchBank("PdbParameterMapContainerBank", bankID,
                                               tablename, TSearch);
  if (NewBank)
  {
    PdbParameterMapContainer *myparm = (PdbParameterMapContainer *) &NewBank->getEntry(0);
    FillFrom(myparm, layer);
    delete NewBank;
  }
  else
  {
    cout << PHWHERE " Reading from DB failed" << endl;
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
    cout << PHWHERE << " Aborting, Database not readable" << endl;
    application->abort();
    gSystem->Exit(1);
    exit(1);
  }

  //  Make a bank ID...
  PdbBankID bankID(0);  // lets start at zero
  PHTimeStamp TSearch(10);

  string tablename = m_Detector + "_geoparams";
  std::transform(tablename.begin(), tablename.end(), tablename.begin(),
                 ::tolower);
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
    cout << PHWHERE " Reading from DB failed" << endl;
    return -1;
  }
  return 0;
}

int PHParameters::WriteToFile(const string &extension, const string &dir)
{
  ostringstream fullpath;
  ostringstream fnamestream;
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
              << "-" << TStart.getTics() << "-" << TStop.getTics() << "-" << time(0)
              << "." << extension;
  string fname = fnamestream.str();
  std::transform(fname.begin(), fname.end(), fname.begin(), ::tolower);
  fullpath << fname;

  cout << "PHParameters::WriteToFile - save to " << fullpath.str() << endl;

  PdbParameterMap *myparm = new PdbParameterMap();
  CopyToPdbParameterMap(myparm);
  TFile *f = TFile::Open(fullpath.str().c_str(), "recreate");
  // force xml file writing to use extended precision shown experimentally
  // to not modify input parameters (.17g)
  string floatformat = TBufferXML::GetFloatFormat();
  TBufferXML::SetFloatFormat("%.17g");  // for IEEE 754 double
  myparm->Write();
  delete f;
  // restore previous xml float format
  TBufferXML::SetFloatFormat(floatformat.c_str());
  cout << "sleeping 1 second to prevent duplicate inserttimes" << endl;
  sleep(1);
  return 0;
}

int PHParameters::ReadFromFile(const string &name, const string &extension, const int layer, const int issuper, const string &dir)
{
  PHTimeStamp TSearch(10);
  PdbBankID bankID(0);
  ostringstream fnamestream;
  fnamestream << name << "_geoparams"
              << "-" << bankID.getInternalValue();
  string fileprefix = fnamestream.str();
  std::transform(fileprefix.begin(), fileprefix.end(), fileprefix.begin(),
                 ::tolower);
  boost::filesystem::path targetDir(dir);

  boost::filesystem::recursive_directory_iterator iter(targetDir), eod;
  boost::char_separator<char> sep("-.");
  map<unsigned int, string> calibfiles;
  BOOST_FOREACH (boost::filesystem::path const &i, make_pair(iter, eod))
  {
    if (is_regular_file(i))
    {
      // boost leaf() gives the filename without path,
      // this checks if the filename starts with fileprefix
      // (start pos of substring=0), if not coninue
      string basename = i.filename().string();
      if (basename.find(fileprefix) != 0)
      {
        continue;
      }
      // extension() contains the . - like .xml, so we
      // just compare the extensions instead of !=
      // and check that the size makes sense
      if (i.extension().string().find(extension) == string::npos || i.extension().string().size() != extension.size() + 1)
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
    cout << "No calibration file like " << dir << "/" << fileprefix << " found" << endl;
    gSystem->Exit(1);
  }
  cout << "PHParameters::ReadFromFile - Reading from File: " << (calibfiles.rbegin())->second << " ... ";
  string fname = (calibfiles.rbegin())->second;
  TFile *f = TFile::Open(fname.c_str());
  if (issuper)
  {
    PdbParameterMapContainer *myparm = static_cast<PdbParameterMapContainer *>(f->Get("PdbParameterMapContainer"));
    assert(myparm);

    if (myparm->GetParameters(layer) == nullptr)
      cout << "Missing PdbParameterMapContainer layer #" << layer << endl;
    assert(myparm->GetParameters(layer));

    cout << "Received PdbParameterMapContainer layer #" << layer << " with (Hash = 0x" << std::hex << myparm->GetParameters(layer)->get_hash() << std::dec << ")" << endl;

    FillFrom(myparm, layer);
    delete myparm;
  }
  else
  {
    PdbParameterMap *myparm = static_cast<PdbParameterMap *>(f->Get("PdbParameterMap"));
    assert(myparm);
    cout << "Received PdbParameterMap with (Hash = 0x" << std::hex << myparm->get_hash() << std::dec << ")" << endl;

    FillFrom(myparm);
    delete myparm;
  }
  delete f;

  return 0;
}

void PHParameters::CopyToPdbParameterMap(PdbParameterMap *myparm)
{
  for (map<const string, double>::const_iterator iter = m_DoubleParMap.begin();
       iter != m_DoubleParMap.end(); ++iter)
  {
    myparm->set_double_param(iter->first, iter->second);
  }
  for (map<const string, int>::const_iterator iter = m_IntParMap.begin();
       iter != m_IntParMap.end(); ++iter)
  {
    myparm->set_int_param(iter->first, iter->second);
  }
  for (map<const string, string>::const_iterator iter = m_StringParMap.begin();
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
    cout << "Cannot extract timestamp from " << str << endl;
    gSystem->Exit(1);
    exit(1);
  }
  return tics;
}
