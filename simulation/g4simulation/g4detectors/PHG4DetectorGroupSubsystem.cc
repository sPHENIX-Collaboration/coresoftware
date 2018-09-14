#include "PHG4DetectorGroupSubsystem.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <pdbcalbase/PdbParameterMap.h>
#include <pdbcalbase/PdbParameterMapContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TSystem.h>

#include <boost/format.hpp>

#include <iostream>
#include <sstream>

using namespace std;

PHG4DetectorGroupSubsystem::PHG4DetectorGroupSubsystem(const std::string &name, const int lyr)
  : PHG4Subsystem(name)
  , m_ParamsContainer(nullptr)
  , m_ParamsContainerDefault(new PHParametersContainer(Name()))
  , m_SaveTopNode(nullptr)
  , m_OverlapCheckFlag(false)
  , m_Layer(lyr)
  , m_UseDBFlag(0)
  , m_BeginRunExecutedFlag(0)
  , m_FileType(PHG4DetectorGroupSubsystem::none)
  , m_SuperDetector("NONE")
  , m_CalibFileDir("./")
{
  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str().c_str());
}

int PHG4DetectorGroupSubsystem::Init(PHCompositeNode *topNode)
{
  m_SaveTopNode = topNode;
  m_ParamsContainerDefault->set_name(Name());
  int iret = InitSubsystem(topNode);
  return iret;
}

int PHG4DetectorGroupSubsystem::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));

  string g4geonodename = "G4GEO_";
  string paramnodename = "G4GEOPARAM_";
  string calibdetname;
  int isSuperDetector = 0;
  if (m_SuperDetector != "NONE")
  {
    g4geonodename += SuperDetector();
    m_ParamsContainer = findNode::getClass<PHParametersContainer>(parNode, g4geonodename);
    if (!m_ParamsContainer)
    {
      PHNodeIterator parIter(parNode);
      PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", SuperDetector()));
      if (!DetNode)
      {
        DetNode = new PHCompositeNode(SuperDetector());
        parNode->addNode(DetNode);
      }
      m_ParamsContainer = new PHParametersContainer(m_SuperDetector);
      DetNode->addNode(new PHDataNode<PHParametersContainer>(m_ParamsContainer, g4geonodename));
    }
    paramnodename += m_SuperDetector;
    calibdetname = m_SuperDetector;
    isSuperDetector = 1;
  }
  else
  {
    m_ParamsContainer = new PHParametersContainer(Name());
    g4geonodename += m_ParamsContainer->Name();
    parNode->addNode(new PHDataNode<PHParametersContainer>(m_ParamsContainer, g4geonodename));
    paramnodename += m_ParamsContainer->Name();
    calibdetname = m_ParamsContainer->Name();
  }

  PHParametersContainer::ConstRange begin_end = m_ParamsContainerDefault->GetAllParameters();
  for (PHParametersContainer::ConstIterator iter = begin_end.first; iter != begin_end.second; ++iter)
  {
    m_ParamsContainer->AddPHParameters(iter->first, iter->second);
  }
  // the content has been handed off to the param container on the node tree
  // clear our internal map of parameters and delete it to avoid it being used accidentally
  m_ParamsContainerDefault->clear();
  delete m_ParamsContainerDefault;
  m_ParamsContainerDefault = nullptr;
  // ASSUMPTION: if we read from DB and/or file we don't want the stuff from
  // the node tree
  // We leave the defaults intact in case there is no entry for
  // those in the object read from the DB or file
  // Order: read first DB, then calib file if both are enabled
  if (ReadDB() || get_filetype() != PHG4DetectorGroupSubsystem::none)
  {
    if (ReadDB())
    {
      ReadParamsFromDB(calibdetname, isSuperDetector);
    }
    if (get_filetype() != PHG4DetectorGroupSubsystem::none)
    {
      ReadParamsFromFile(calibdetname, get_filetype(), isSuperDetector);
    }
  }
  else
  {
    // if not filled from file or DB, check if we have a node containing those calibrations
    // on the node tree and load them (the embedding wants to use the
    // parameters saved on the previous pass)
    PdbParameterMapContainer *nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, paramnodename);
    if (nodeparams)
    {
      m_ParamsContainer->FillFrom(nodeparams);
    }
  }
  // parameters set in the macro always override whatever is read from
  // the node tree, DB or file
  UpdateParametersWithMacro();
  // save updated persistant copy on node tree
  PHCompositeNode *RunDetNode = runNode;
  if (m_SuperDetector != "NONE")
  {
    PHNodeIterator runIter(runNode);
    RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!RunDetNode)
    {
      RunDetNode = new PHCompositeNode(SuperDetector());
      runNode->addNode(RunDetNode);
    }
  }
  m_ParamsContainer->SaveToNodeTree(RunDetNode, paramnodename);
  int iret = InitRunSubsystem(topNode);
  if (Verbosity() > 0)
  {
    PdbParameterMapContainer *nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, paramnodename);
    cout << Name() << endl;
    nodeparams->print();
  }
  m_BeginRunExecutedFlag = 1;
  return iret;
}

void PHG4DetectorGroupSubsystem::SuperDetector(const std::string &name)
{
  m_SuperDetector = name;
  return;
}

double PHG4DetectorGroupSubsystem::get_double_param(const int detid, const std::string &name) const
{
  const PHParameters *params = m_ParamsContainer->GetParameters(detid);
  if (params)
  {
    return params->get_double_param(name);
  }
  else
  {
    cout << PHWHERE << " no parameters for detid " << detid << endl;
    gSystem->Exit(1);
    exit(1);
  }
}

int PHG4DetectorGroupSubsystem::get_int_param(const int detid, const std::string &name) const
{
  const PHParameters *params = m_ParamsContainer->GetParameters(detid);
  if (params)
  {
    return params->get_int_param(name);
  }
  else
  {
    cout << PHWHERE << " no parameters for detid " << detid << endl;
    gSystem->Exit(1);
    exit(1);
  }
}

string
PHG4DetectorGroupSubsystem::get_string_param(const int detid, const std::string &name) const
{
  const PHParameters *params = m_ParamsContainer->GetParameters(detid);
  if (params)
  {
    return params->get_string_param(name);
  }
  else
  {
    cout << PHWHERE << " no parameters for detid " << detid << endl;
    gSystem->Exit(1);
    exit(1);
  }
}

void PHG4DetectorGroupSubsystem::set_double_param(const int detid, const std::string &name, const double dval)
{
  auto iter = m_DefaultDoubleParamsMap.find(detid);
  if (iter == m_DefaultDoubleParamsMap.end())
  {
    cout << "called like set_double_param(" << detid << ", \""
         << name << "\", " << dval << ")" << endl;
    cout << "detid " << detid << " not implemented" << endl;
    cout << "implemented detector ids: " << endl;
    for (auto &iter2 : m_DefaultDoubleParamsMap)
    {
      cout << "detid: " << iter2.first << endl;
    }
    return;
  }
  if (iter->second.find(name) == iter->second.end())
  {
    cout << "double parameter " << name << " not implemented for detid "
         << detid << endl;
    cout << "implemented double parameters are:" << endl;
    for (auto &iter2 : iter->second)
    {
      cout << iter2.first << endl;
    }
    return;
  }
  // here we know we have entries for the detector id and the variable name exists
  // in the defaults, so now lets set it
  // with C++11 insert returns a pair of an iterator to the element and
  // a boolean if the object was inserted or if it already exist (in which case it
  // does not get replaced). We do not check this because we do not care if a new map
  // was inserted or not. All we need is the iterator to it
  map<const std::string, double> newdmap;
  auto ret = m_MacroDoubleParamsMap.insert(make_pair(detid, newdmap));
  // here we use the operator [] rather than insert because we
  // want to create a new entry if [name] does not exist. If it does
  // exist we want to overwrite it (so even if a parameter is set twice,
  // the last setting in the macro is used). Using insert would preserve the first
  // parameter setting
  ret.first->second[name] = dval;
  return;
}

void PHG4DetectorGroupSubsystem::set_int_param(const int detid, const std::string &name, const int ival)
{
  auto iter = m_DefaultIntegerParamsMap.find(detid);
  if (iter == m_DefaultIntegerParamsMap.end())
  {
    cout << "called like set_int_param(" << detid << ", \""
         << name << "\", " << ival << ")" << endl;
    cout << "detid " << detid << " not implemented" << endl;
    cout << "implemented detector ids: " << endl;
    for (auto &iter2 : m_DefaultIntegerParamsMap)
    {
      cout << "detid: " << iter2.first << endl;
    }
    return;
  }
  if (iter->second.find(name) == iter->second.end())
  {
    cout << "int parameter " << name << " not implemented for detid"
         << detid << endl;
    cout << "implemented int parameters are:" << endl;
    for (auto &iter2 : iter->second)
    {
      cout << iter2.first << endl;
    }
    return;
  }
  // here we know we have entries for the detector id and the variable name exists
  // in the defaults, so now lets set it
  map<const std::string, int> newintmap;
  auto ret = m_MacroIntegerParamsMap.insert(make_pair(detid, newintmap));
  // here we use the operator [] rather than insert because we
  // want to create a new entry if [name] does not exist. If it does
  // exist we want to overwrite it (so even if a parameter is set twice,
  // the last setting in the macro is used). Using insert would preserve the first
  // parameter setting
  ret.first->second[name] = ival;
}

void PHG4DetectorGroupSubsystem::set_string_param(const int detid, const std::string &name, const string &sval)
{
  auto iter = m_DefaultStringParamsMap.find(detid);
  if (iter == m_DefaultStringParamsMap.end())
  {
    cout << "called like set_string_param(" << detid << ", \""
         << name << "\", " << sval << ")" << endl;
    cout << "detid " << detid << " not implemented" << endl;
    cout << "implemented detector ids: " << endl;
    for (auto &iter2 : m_DefaultStringParamsMap)
    {
      cout << "detid: " << iter2.first << endl;
    }
    return;
  }
  if (iter->second.find(name) == iter->second.end())
  {
    cout << "string parameter " << name << " not implemented for detid "
         << detid << endl;
    cout << "implemented string parameters are:" << endl;
    for (auto &iter2 : iter->second)
    {
      cout << iter2.first << endl;
    }
    return;
  }
  // here we know we have entries for the detector id and the variable name exists
  // in the defaults, so now lets set it

  // with C++11 insert returns a pair of an iterator to the element and
  // a boolean if the object was inserted or if it already exist (in which case it
  // does not get replaced). We do not check this because we do not care if a new map
  // was inserted or not. All we need is the iterator to it
  map<const std::string, string> newdmap;
  auto ret = m_MacroStringParamsMap.insert(make_pair(detid, newdmap));
  // here we use the operator [] rather than insert because we
  // want to create a new entry if [name] does not exist. If it does
  // exist we want to overwrite it (so even if a parameter is set twice,
  // the last setting in the macro is used). Using insert would preserve the first
  // parameter setting
  ret.first->second[name] = sval;
  return;
}

void PHG4DetectorGroupSubsystem::UpdateParametersWithMacro()
{
  for (auto &iter1 : m_MacroDoubleParamsMap)
  {
    PHParameters *params = GetParamsContainer()->GetParametersToModify(iter1.first);
    for (auto &iter2 : iter1.second)
    {
      params->set_double_param(iter2.first, iter2.second);
    }
  }

  for (auto &iter1 : m_MacroIntegerParamsMap)
  {
    PHParameters *params = GetParamsContainer()->GetParametersToModify(iter1.first);
    for (auto &iter2 : iter1.second)
    {
      params->set_int_param(iter2.first, iter2.second);
    }
  }

  for (auto &iter1 : m_MacroStringParamsMap)
  {
    PHParameters *params = GetParamsContainer()->GetParametersToModify(iter1.first);
    for (auto &iter2 : iter1.second)
    {
      params->set_string_param(iter2.first, iter2.second);
    }
  }
  return;
}

// all set_default_XXX_param methods work the same way
// They use the returned pair <iterator, bool> for c++11 map.insert. The boolean tells us
// if a new element was inserted (true) or an iterator to the exisiting one is returned (false)
// map.insert does not overwrite existing entries. So what we are doing here is
// get an iterator to the double/int/string map for a given detector (for us it does not
// matter if it already exists or a new one is inserted, we just need an iterator to it)
// Then we use map.insert to insert the new double/int/string parameter into the
// double/int/string map. If
// the return boolean is false, it means an entry for this parameter already exists
// which is just bad (means you called set_default_XXX_param for the same parameter
// twice in your code), so tell the user to fix the code and exit

void PHG4DetectorGroupSubsystem::set_default_double_param(const int detid, const std::string &name, const double dval)
{
  map<const std::string, double> newdoublemap;
  auto ret = m_DefaultDoubleParamsMap.insert(make_pair(detid, newdoublemap));
  auto ret2 = ret.first->second.insert(make_pair(name, dval));
  if (ret2.second == false)
  {
    cout << PHWHERE << "Default double Parameter " << name
         << " for detid " << detid << " already set to "
         << ret.first->second[name] << " will not overwrite with " << dval << endl;
    cout << "Means: You are calling set_default_double_param twice for the same parameter" << endl;
    cout << "Please make up your mind and call it only once using the correct default" << endl;
    gSystem->Exit(1);
    exit(1);
  }
  return;
}

void PHG4DetectorGroupSubsystem::set_default_int_param(const int detid, const std::string &name, const int ival)
{
  map<const std::string, int> newintmap;
  auto ret = m_DefaultIntegerParamsMap.insert(make_pair(detid, newintmap));
  auto ret2 = ret.first->second.insert(make_pair(name, ival));
  if (ret2.second == false)
  {
    cout << PHWHERE << "Default integer Parameter " << name
         << " for detid " << detid << " already set to "
         << ret.first->second[name] << " will not overwrite with " << ival << endl;
    cout << "Means: You are calling set_default_int_param twice for the same parameter" << endl;
    cout << "Please make up your mind and call it only once using the correct default" << endl;
    gSystem->Exit(1);
    exit(1);
  }
  return;
}

void PHG4DetectorGroupSubsystem::set_default_string_param(const int detid, const std::string &name, const string &sval)
{
  map<const std::string, string> newstringmap;
  auto ret = m_DefaultStringParamsMap.insert(make_pair(detid, newstringmap));
  auto ret2 = ret.first->second.insert(make_pair(name, sval));
  if (ret2.second == false)
  {
    cout << PHWHERE << "Default String Parameter " << name
         << " for detid " << detid << " already set to "
         << ret.first->second[name] << " will not overwrite with " << sval << endl;
    cout << "Means: You are calling set_default_string_param twice for the same parameter" << endl;
    cout << "Please make up your mind and call it only once using the correct default" << endl;
    gSystem->Exit(1);
    exit(1);
  }
  return;
}

void PHG4DetectorGroupSubsystem::InitializeParameters()
{
  for (auto &iter : m_LayerSet)
  {
    set_default_int_param(iter, "absorberactive", 0);
    set_default_int_param(iter, "absorbertruth", 0);
    set_default_int_param(iter, "blackhole", 0);
  }
  SetDefaultParameters();  // call method from specific subsystem
  // now load those parameters to our params class
  for (auto iter1 : m_DefaultDoubleParamsMap)
  {
    PHParameters *detidparams = m_ParamsContainerDefault->GetParametersToModify(iter1.first);
    if (!detidparams)
    {
      detidparams = new PHParameters((boost::format("%s_%d") % Name() % iter1.first).str());
      m_ParamsContainerDefault->AddPHParameters(iter1.first, detidparams);
    }
    for (auto &iter2 : iter1.second)
    {
      detidparams->set_double_param(iter2.first, iter2.second);
    }
  }

  for (auto &iter1 : m_DefaultIntegerParamsMap)
  {
    PHParameters *detidparams = m_ParamsContainerDefault->GetParametersToModify(iter1.first);
    if (!detidparams)
    {
      detidparams = new PHParameters((boost::format("%s_%d") % Name() % iter1.first).str());
      m_ParamsContainerDefault->AddPHParameters(iter1.first, detidparams);
    }
    for (auto &iter2 : iter1.second)
    {
      detidparams->set_int_param(iter2.first, iter2.second);
    }
  }

  for (auto &iter1 : m_DefaultStringParamsMap)
  {
    PHParameters *detidparams = m_ParamsContainerDefault->GetParametersToModify(iter1.first);
    if (!detidparams)
    {
      detidparams = new PHParameters((boost::format("%s_%d") % Name() % iter1.first).str());
      m_ParamsContainerDefault->AddPHParameters(iter1.first, detidparams);
    }
    for (auto &iter2 : iter1.second)
    {
      detidparams->set_string_param(iter2.first, iter2.second);
    }
  }
}

int PHG4DetectorGroupSubsystem::SaveParamsToDB()
{
  int iret = 0;
  assert(m_ParamsContainer);
  iret = m_ParamsContainer->WriteToDB();
  if (iret)
  {
    cout << "problem committing to DB" << endl;
  }
  return iret;
}

int PHG4DetectorGroupSubsystem::ReadParamsFromDB(const string &name, const int issuper)
{
  int iret = 1;
  // if (issuper)
  //   {
  //     iret = params->ReadFromDB(name,layer);
  //   }
  // else
  //   {
  //     iret = params->ReadFromDB();
  //   }
  if (iret)
  {
    cout << "problem reading from DB" << endl;
  }
  return iret;
}

int PHG4DetectorGroupSubsystem::SaveParamsToFile(const PHG4DetectorGroupSubsystem::FILE_TYPE ftyp)
{
  string extension;
  switch (ftyp)
  {
  case xml:
    extension = "xml";
    break;
  case root:
    extension = "root";
    break;
  default:
    cout << PHWHERE << "filetype " << ftyp << " not implemented" << endl;
    exit(1);
  }
  int iret = 0;
  assert(m_ParamsContainer);
  iret = m_ParamsContainer->WriteToFile(extension, m_CalibFileDir);
  if (iret)
  {
    cout << "problem saving to " << extension << " file " << endl;
  }
  return iret;
}

int PHG4DetectorGroupSubsystem::ReadParamsFromFile(const string &name, const PHG4DetectorGroupSubsystem::FILE_TYPE ftyp, const int issuper)
{
  string extension;
  switch (ftyp)
  {
  case xml:
    extension = "xml";
    break;
  case root:
    extension = "root";
    break;
  default:
    cout << PHWHERE << "filetype " << ftyp << " not implemented" << endl;
    exit(1);
  }
  int iret = 1;
  //  int iret = params->ReadFromFile(name, extension, layer, issuper, m_CalibFileDir);
  if (iret)
  {
    cout << "problem reading from " << extension << " file " << endl;
  }
  return iret;
}

void PHG4DetectorGroupSubsystem::SetActive(const int detid, const int i)
{
  set_int_param(detid, "active", i);
}

void PHG4DetectorGroupSubsystem::SetActive(const int i)
{
  for (auto &detid : m_LayerSet)
  {
    set_int_param(detid, "active", i);
  }
}

void PHG4DetectorGroupSubsystem::SetAbsorberActive(const int detid, const int i)
{
  set_int_param(detid, "absorberactive", i);
}

void PHG4DetectorGroupSubsystem::SetAbsorberActive(const int i)
{
  for (auto &detid : m_LayerSet)
  {
    set_int_param(detid, "absorberactive", i);
  }
}

void PHG4DetectorGroupSubsystem::BlackHole(const int detid, const int i)
{
  set_int_param(detid, "blackhole", i);
}

void PHG4DetectorGroupSubsystem::BlackHole(const int i)
{
  for (auto &detid : m_LayerSet)
  {
    set_int_param(detid, "blackhole", i);
  }
}

void PHG4DetectorGroupSubsystem::SetAbsorberTruth(const int detid, const int i)
{
  set_int_param(detid, "absorbertruth", i);
}

void PHG4DetectorGroupSubsystem::SetAbsorberTruth(const int i)
{
  for (auto &detid : m_LayerSet)
  {
    set_int_param(detid, "absorbertruth", i);
  }
}

void PHG4DetectorGroupSubsystem::PrintDefaultParams() const
{
  cout << "Default Parameters: " << endl;
  cout << "int values: " << endl;
  for (auto &iter1 : m_DefaultIntegerParamsMap)
  {
    cout << "Detector id: " << iter1.first << endl;
    for (auto &iter2 : iter1.second)
    {
      cout << iter2.first << " : " << iter2.second << endl;
    }
  }
  cout << "double values: " << endl;
  for (auto &iter1 : m_DefaultDoubleParamsMap)
  {
    cout << "Detector id: " << iter1.first << endl;
    for (auto &iter2 : iter1.second)
    {
      cout << iter2.first << " : " << iter2.second << endl;
    }
  }
  cout << "string values: " << endl;
  for (auto &iter1 : m_DefaultStringParamsMap)
  {
    cout << "Detector id: " << iter1.first << endl;
    for (auto &iter2 : iter1.second)
    {
      cout << iter2.first << " : " << iter2.second << endl;
    }
  }
  return;
}

void PHG4DetectorGroupSubsystem::PrintMacroParams() const
{
  cout << "Macro Parameters: " << endl;
  cout << "int values: " << endl;
  for (auto &iter1 : m_MacroIntegerParamsMap)
  {
    cout << "Detector id: " << iter1.first << endl;
    for (auto &iter2 : iter1.second)
    {
      cout << iter2.first << " : " << iter2.second << endl;
    }
  }
  cout << "double values: " << endl;
  for (auto &iter1 : m_MacroDoubleParamsMap)
  {
    cout << "Detector id: " << iter1.first << endl;
    for (auto &iter2 : iter1.second)
    {
      cout << iter2.first << " : " << iter2.second << endl;
    }
  }
  cout << "string values: " << endl;
  for (auto &iter1 : m_MacroStringParamsMap)
  {
    cout << "Detector id: " << iter1.first << endl;
    for (auto &iter2 : iter1.second)
    {
      cout << iter2.first << " : " << iter2.second << endl;
    }
  }
  return;
}
