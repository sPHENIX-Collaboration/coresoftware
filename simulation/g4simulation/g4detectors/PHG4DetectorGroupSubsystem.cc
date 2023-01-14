#include "PHG4DetectorGroupSubsystem.h"

#include <g4main/PHG4Subsystem.h>  // for PHG4Subsystem

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <pdbcalbase/PdbParameterMapContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TSystem.h>

#include <boost/format.hpp>

// boost stacktrace has shadowed variables
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/stacktrace.hpp>
#pragma GCC diagnostic pop

#include <cassert>  // for assert
#include <cstdlib>  // for exit
#include <iostream>
#include <sstream>

PHG4DetectorGroupSubsystem::PHG4DetectorGroupSubsystem(const std::string &name, const int lyr)
  : PHG4Subsystem(name)
  , m_ParamsContainerDefault(new PHParametersContainer(Name()))
  , m_Layer(lyr)
{
  // put the layer into the name so we get unique names
  // for multiple layers
  std::ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str());
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

  std::string g4geonodename = "G4GEO_";
  std::string paramnodename = "G4GEOPARAM_";
  std::string calibdetname;
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
  for (PHParametersContainer::ConstIterator piter = begin_end.first; piter != begin_end.second; ++piter)
  {
    m_ParamsContainer->AddPHParameters(piter->first, piter->second);
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
  // define the materials for the detector
  // at this point all flags are known so materials set in the macro can
  // be implemented here
  DefineMaterials();
  int iret = InitRunSubsystem(topNode);
  m_ParamsContainer->UpdateNodeTree(RunDetNode, paramnodename);
  if (Verbosity() > 0)
  {
    PdbParameterMapContainer *nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, paramnodename);
    std::cout << Name() << std::endl;
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
    std::cout << PHWHERE << " no parameters for detid " << detid << std::endl;
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
    std::cout << PHWHERE << " no parameters for detid " << detid << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
}

std::string
PHG4DetectorGroupSubsystem::get_string_param(const int detid, const std::string &name) const
{
  const PHParameters *params = m_ParamsContainer->GetParameters(detid);
  if (params)
  {
    return params->get_string_param(name);
  }
  else
  {
    std::cout << PHWHERE << " no parameters for detid " << detid << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
}

void PHG4DetectorGroupSubsystem::set_double_param(const int detid, const std::string &name, const double dval)
{
  auto iter = m_DefaultDoubleParamsMap.find(detid);
  if (iter == m_DefaultDoubleParamsMap.end())
  {
    std::cout << "called like set_double_param(" << detid << ", \""
              << name << "\", " << dval << ")" << std::endl;
    std::cout << "detid " << detid << " not implemented" << std::endl;
    std::cout << "implemented detector ids: " << std::endl;
    for (auto &iter2 : m_DefaultDoubleParamsMap)
    {
      std::cout << "detid: " << iter2.first << std::endl;
    }
    return;
  }
  if (iter->second.find(name) == iter->second.end())
  {
    std::cout << "double parameter " << name << " not implemented for detid "
              << detid << std::endl;
    std::cout << "implemented double parameters are:" << std::endl;
    for (auto &iter2 : iter->second)
    {
      std::cout << iter2.first << std::endl;
    }
    return;
  }
  // here we know we have entries for the detector id and the variable name exists
  // in the defaults, so now lets set it
  // with C++11 insert returns a pair of an iterator to the element and
  // a boolean if the object was inserted or if it already exist (in which case it
  // does not get replaced). We do not check this because we do not care if a new map
  // was inserted or not. All we need is the iterator to it
  std::map<const std::string, double> newdmap;
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
    std::cout << "called like set_int_param(" << detid << ", \""
              << name << "\", " << ival << ")" << std::endl;
    std::cout << "detid " << detid << " not implemented" << std::endl;
    std::cout << "implemented detector ids: " << std::endl;
    for (auto &iter2 : m_DefaultIntegerParamsMap)
    {
      std::cout << "detid: " << iter2.first << std::endl;
    }
    return;
  }
  if (iter->second.find(name) == iter->second.end())
  {
    std::cout << "int parameter " << name << " not implemented for detid"
              << detid << std::endl;
    std::cout << "implemented int parameters are:" << std::endl;
    for (auto &iter2 : iter->second)
    {
      std::cout << iter2.first << std::endl;
    }
    return;
  }
  // here we know we have entries for the detector id and the variable name exists
  // in the defaults, so now lets set it
  std::map<const std::string, int> newintmap;
  auto ret = m_MacroIntegerParamsMap.insert(make_pair(detid, newintmap));
  // here we use the operator [] rather than insert because we
  // want to create a new entry if [name] does not exist. If it does
  // exist we want to overwrite it (so even if a parameter is set twice,
  // the last setting in the macro is used). Using insert would preserve the first
  // parameter setting
  ret.first->second[name] = ival;
}

void PHG4DetectorGroupSubsystem::set_string_param(const int detid, const std::string &name, const std::string &sval)
{
  auto iter = m_DefaultStringParamsMap.find(detid);
  if (iter == m_DefaultStringParamsMap.end())
  {
    std::cout << "called like set_string_param(" << detid << ", \""
              << name << "\", " << sval << ")" << std::endl;
    std::cout << "detid " << detid << " not implemented" << std::endl;
    std::cout << "implemented detector ids: " << std::endl;
    for (auto &iter2 : m_DefaultStringParamsMap)
    {
      std::cout << "detid: " << iter2.first << std::endl;
    }
    return;
  }
  if (iter->second.find(name) == iter->second.end())
  {
    std::cout << "string parameter " << name << " not implemented for detid "
              << detid << std::endl;
    std::cout << "implemented string parameters are:" << std::endl;
    for (auto &iter2 : iter->second)
    {
      std::cout << iter2.first << std::endl;
    }
    return;
  }
  // here we know we have entries for the detector id and the variable name exists
  // in the defaults, so now lets set it

  // with C++11 insert returns a pair of an iterator to the element and
  // a boolean if the object was inserted or if it already exist (in which case it
  // does not get replaced). We do not check this because we do not care if a new map
  // was inserted or not. All we need is the iterator to it
  std::map<const std::string, std::string> newdmap;
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
  std::map<const std::string, double> newdoublemap;
  auto ret = m_DefaultDoubleParamsMap.insert(make_pair(detid, newdoublemap));
  auto ret2 = ret.first->second.insert(make_pair(name, dval));
  if (ret2.second == false)
  {
    std::cout << PHWHERE << "Default double Parameter " << name
              << " for detid " << detid << " already set to "
              << ret.first->second[name] << " will not overwrite with " << dval << std::endl;
    std::cout << "Means: You are calling set_default_double_param twice for the same parameter" << std::endl;
    std::cout << "Please make up your mind and call it only once using the correct default" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  return;
}

void PHG4DetectorGroupSubsystem::set_default_int_param(const int detid, const std::string &name, const int ival)
{
  std::map<const std::string, int> newintmap;
  auto ret = m_DefaultIntegerParamsMap.insert(make_pair(detid, newintmap));
  auto ret2 = ret.first->second.insert(make_pair(name, ival));
  if (ret2.second == false)
  {
    std::cout << PHWHERE << "Default integer Parameter " << name
              << " for detid " << detid << " already set to "
              << ret.first->second[name] << " will not overwrite with " << ival << std::endl;
    std::cout << "Means: You are calling set_default_int_param twice for the same parameter" << std::endl;
    std::cout << "Please make up your mind and call it only once using the correct default" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  return;
}

void PHG4DetectorGroupSubsystem::set_default_string_param(const int detid, const std::string &name, const std::string &sval)
{
  std::map<const std::string, std::string> newstringmap;
  auto ret = m_DefaultStringParamsMap.insert(make_pair(detid, newstringmap));
  auto ret2 = ret.first->second.insert(make_pair(name, sval));
  if (ret2.second == false)
  {
    std::cout << PHWHERE << "Default String Parameter " << name
              << " for detid " << detid << " already set to "
              << ret.first->second[name] << " will not overwrite with " << sval << std::endl;
    std::cout << "Means: You are calling set_default_string_param twice for the same parameter" << std::endl;
    std::cout << "Please make up your mind and call it only once using the correct default" << std::endl;
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
    set_default_int_param(iter, "supportactive", 0);
  }
  SetDefaultParameters();  // call method from specific subsystem
  // now load those parameters to our params class
  for (const auto &iter1 : m_DefaultDoubleParamsMap)
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
    std::cout << "problem committing to DB" << std::endl;
  }
  return iret;
}

int PHG4DetectorGroupSubsystem::ReadParamsFromDB(const std::string & /*name*/, const int /*issuper*/)
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
  //  if (iret)
  std::cout << boost::stacktrace::stacktrace();
  std::cout << std::endl
            << "DO NOT PANIC - this is not a segfault" << std::endl;
  std::cout << "This method is a dummy, tell the offline gurus about it and give this stack trace" << std::endl;
  {
    std::cout << "problem reading from DB" << std::endl;
  }
  gSystem->Exit(1);
  return iret;
}

int PHG4DetectorGroupSubsystem::SaveParamsToFile(const PHG4DetectorGroupSubsystem::FILE_TYPE ftyp)
{
  std::string extension;
  switch (ftyp)
  {
  case xml:
    extension = "xml";
    break;
  case root:
    extension = "root";
    break;
  default:
    std::cout << PHWHERE << "filetype " << ftyp << " not implemented" << std::endl;
    exit(1);
  }
  int iret = 0;
  assert(m_ParamsContainer);
  iret = m_ParamsContainer->WriteToFile(extension, m_CalibFileDir);
  if (iret)
  {
    std::cout << "problem saving to " << extension << " file " << std::endl;
  }
  return iret;
}

int PHG4DetectorGroupSubsystem::ReadParamsFromFile(const std::string & /*name*/, const PHG4DetectorGroupSubsystem::FILE_TYPE ftyp, const int /*issuper*/)
{
  std::string extension;
  switch (ftyp)
  {
  case xml:
    extension = "xml";
    break;
  case root:
    extension = "root";
    break;
  default:
    std::cout << PHWHERE << "filetype " << ftyp << " not implemented" << std::endl;
    exit(1);
  }
  int iret = 1;
  //  int iret = params->ReadFromFile(name, extension, layer, issuper, m_CalibFileDir);
  //  if (iret)
  std::cout << boost::stacktrace::stacktrace();
  std::cout << std::endl
            << "DO NOT PANIC - this is not a segfault" << std::endl;
  std::cout << "This method is a dummy, tell the offline gurus about it and give this stack trace" << std::endl;
  std::cout << "problem reading from " << extension << " file " << std::endl;
  std::cout << "problem reading from DB" << std::endl;
  gSystem->Exit(1);
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

void PHG4DetectorGroupSubsystem::SetSupportActive(const int detid, const int i)
{
  set_int_param(detid, "supportactive", i);
}

void PHG4DetectorGroupSubsystem::SetSupportActive(const int i)
{
  for (auto &detid : m_LayerSet)
  {
    set_int_param(detid, "supportactive", i);
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
  std::cout << "Default Parameters: " << std::endl;
  std::cout << "int values: " << std::endl;
  for (auto &iter1 : m_DefaultIntegerParamsMap)
  {
    std::cout << "Detector id: " << iter1.first << std::endl;
    for (auto &iter2 : iter1.second)
    {
      std::cout << iter2.first << " : " << iter2.second << std::endl;
    }
  }
  std::cout << "double values: " << std::endl;
  for (auto &iter1 : m_DefaultDoubleParamsMap)
  {
    std::cout << "Detector id: " << iter1.first << std::endl;
    for (auto &iter2 : iter1.second)
    {
      std::cout << iter2.first << " : " << iter2.second << std::endl;
    }
  }
  std::cout << "string values: " << std::endl;
  for (auto &iter1 : m_DefaultStringParamsMap)
  {
    std::cout << "Detector id: " << iter1.first << std::endl;
    for (auto &iter2 : iter1.second)
    {
      std::cout << iter2.first << " : " << iter2.second << std::endl;
    }
  }
  return;
}

void PHG4DetectorGroupSubsystem::PrintMacroParams() const
{
  std::cout << "Macro Parameters: " << std::endl;
  std::cout << "int values: " << std::endl;
  for (auto &iter1 : m_MacroIntegerParamsMap)
  {
    std::cout << "Detector id: " << iter1.first << std::endl;
    for (auto &iter2 : iter1.second)
    {
      std::cout << iter2.first << " : " << iter2.second << std::endl;
    }
  }
  std::cout << "double values: " << std::endl;
  for (auto &iter1 : m_MacroDoubleParamsMap)
  {
    std::cout << "Detector id: " << iter1.first << std::endl;
    for (auto &iter2 : iter1.second)
    {
      std::cout << iter2.first << " : " << iter2.second << std::endl;
    }
  }
  std::cout << "string values: " << std::endl;
  for (auto &iter1 : m_MacroStringParamsMap)
  {
    std::cout << "Detector id: " << iter1.first << std::endl;
    for (auto &iter2 : iter1.second)
    {
      std::cout << iter2.first << " : " << iter2.second << std::endl;
    }
  }
  return;
}
