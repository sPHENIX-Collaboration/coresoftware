#include "PHG4DetectorSubsystem.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <pdbcalbase/PdbParameterMapContainer.h>

#include <g4main/PHG4Subsystem.h>  // for PHG4Subsystem

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <cstdlib>  // for exit, NULL
#include <iostream>
#include <utility>  // for pair

PHG4DetectorSubsystem::PHG4DetectorSubsystem(const std::string &name, const int lyr)
  : PHG4Subsystem(name)
  , params(new PHParameters(Name()))
  , layer(lyr)
{
  // put the layer into the name so we get unique names
  // for multiple layers
  Name(name + "_" + std::to_string(lyr));
}

int PHG4DetectorSubsystem::Init(PHCompositeNode *topNode)
{
  savetopNode = topNode;
  params->set_name(Name());
  int iret = InitSubsystem(topNode);
  return iret;
}

int PHG4DetectorSubsystem::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));

  std::string g4geonodename = "G4GEO_";
  std::string paramnodename = "G4GEOPARAM_";
  std::string calibdetname;
  int isSuperDetector = 0;
  if (superdetector != "NONE")
  {
    g4geonodename += SuperDetector();
    paramscontainer = findNode::getClass<PHParametersContainer>(parNode, g4geonodename);
    if (!paramscontainer)
    {
      PHNodeIterator parIter(parNode);
      PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", SuperDetector()));
      if (!DetNode)
      {
        DetNode = new PHCompositeNode(SuperDetector());
        parNode->addNode(DetNode);
      }
      paramscontainer = new PHParametersContainer(superdetector);
      DetNode->addNode(new PHDataNode<PHParametersContainer>(paramscontainer, g4geonodename));
    }
    paramscontainer->AddPHParameters(layer, params);
    paramnodename += superdetector;
    calibdetname = superdetector;
    isSuperDetector = 1;
  }
  else
  {
    g4geonodename += params->Name();
    parNode->addNode(new PHDataNode<PHParameters>(params, g4geonodename));
    paramnodename += params->Name();
    calibdetname = params->Name();
  }

  // ASSUMPTION: if we read from DB and/or file we don't want the stuff from
  // the node tree
  // We leave the defaults intact in case there is no entry for
  // those in the object read from the DB or file
  // Order: read first DB, then calib file if both are enabled
  if (ReadDB() || get_filetype() != PHG4DetectorSubsystem::none)
  {
    if (!m_Domain.empty())
    {
      ReadParamsFromCDB(m_Domain);
    }
    else
    {
      if (ReadDB())
      {
        ReadParamsFromDB(calibdetname, isSuperDetector);
      }
      if (get_filetype() != PHG4DetectorSubsystem::none)
      {
        ReadParamsFromFile(calibdetname, get_filetype(), isSuperDetector);
      }
    }
  }
  else
  {
    PdbParameterMapContainer *nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, paramnodename);
    if (nodeparams)
    {
      params->FillFrom(nodeparams, layer);
    }
  }
  // parameters set in the macro always override whatever is read from
  // the node tree, DB or file
  UpdateParametersWithMacro();
  // save updated persistant copy on node tree
  PHCompositeNode *RunDetNode = runNode;
  if (superdetector != "NONE")
  {
    PHNodeIterator runIter(runNode);
    RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!RunDetNode)
    {
      RunDetNode = new PHCompositeNode(SuperDetector());
      runNode->addNode(RunDetNode);
    }
  }
  // put parameters on the node tree so the subsystem can read them from there (it does not have to)
  params->SaveToNodeTree(RunDetNode, paramnodename, layer);
  // define the materials for the detector
  // at this point all flags are known so materials set in the macro can
  // be implemented here
  DefineMaterials();
  int iret = InitRunSubsystem(topNode);
  // update parameters on node tree in case the subsystem has changed them
  // In the end the parameters on the node tree must reflect what was actually used
  params->UpdateNodeTree(RunDetNode, paramnodename, layer);
  if (Verbosity() > 0)
  {
    PdbParameterMapContainer *nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, paramnodename);
    std::cout << Name() << std::endl;
    nodeparams->print();
  }
  beginrunexecuted = 1;
  return iret;
}

void PHG4DetectorSubsystem::SuperDetector(const std::string &name)
{
  superdetector = name;
  return;
}

void PHG4DetectorSubsystem::set_double_param(const std::string &name, const double dval)
{
  if (default_double.find(name) == default_double.end())
  {
    std::cout << "double parameter " << name << " not implemented" << std::endl;
    std::cout << "implemented double parameters are:" << std::endl;
    for (std::map<const std::string, double>::const_iterator iter = default_double.begin(); iter != default_double.end(); ++iter)
    {
      std::cout << iter->first << std::endl;
    }
    return;
  }
  dparams[name] = dval;
}

double
PHG4DetectorSubsystem::get_double_param(const std::string &name) const
{
  return params->get_double_param(name);
}

void PHG4DetectorSubsystem::set_int_param(const std::string &name, const int ival)
{
  if (default_int.find(name) == default_int.end())
  {
    std::cout << "integer parameter " << name << " not implemented" << std::endl;
    std::cout << "implemented integer parameters are:" << std::endl;
    for (std::map<const std::string, int>::const_iterator iter = default_int.begin(); iter != default_int.end(); ++iter)
    {
      std::cout << iter->first << std::endl;
    }
    return;
  }
  iparams[name] = ival;
}

int PHG4DetectorSubsystem::get_int_param(const std::string &name) const
{
  return params->get_int_param(name);
}

void PHG4DetectorSubsystem::set_string_param(const std::string &name, const std::string &sval)
{
  if (default_string.find(name) == default_string.end())
  {
    std::cout << "string parameter " << name << " not implemented" << std::endl;
    std::cout << "implemented string parameters are:" << std::endl;
    for (std::map<const std::string, std::string>::const_iterator iter = default_string.begin(); iter != default_string.end(); ++iter)
    {
      std::cout << iter->first << std::endl;
    }
    return;
  }
  cparams[name] = sval;
}

std::string
PHG4DetectorSubsystem::get_string_param(const std::string &name) const
{
  return params->get_string_param(name);
}

void PHG4DetectorSubsystem::UpdateParametersWithMacro()
{
  for (std::map<const std::string, double>::const_iterator iter = dparams.begin(); iter != dparams.end(); ++iter)
  {
    params->set_double_param(iter->first, iter->second);
  }
  for (std::map<const std::string, int>::const_iterator iter = iparams.begin(); iter != iparams.end(); ++iter)
  {
    params->set_int_param(iter->first, iter->second);
  }
  for (std::map<const std::string, std::string>::const_iterator iter = cparams.begin(); iter != cparams.end(); ++iter)
  {
    params->set_string_param(iter->first, iter->second);
  }
  return;
}

void PHG4DetectorSubsystem::set_default_double_param(const std::string &name, const double dval)
{
  if (default_double.find(name) == default_double.end())
  {
    default_double[name] = dval;
  }
  else
  {
    std::cout << "trying to overwrite default double " << name << " "
              << default_double[name] << " with " << dval << std::endl;
    exit(1);
  }
  return;
}

void PHG4DetectorSubsystem::set_default_int_param(const std::string &name, const int ival)
{
  if (default_int.find(name) == default_int.end())
  {
    default_int[name] = ival;
  }
  else
  {
    std::cout << "trying to overwrite default int " << name << " "
              << default_int[name] << " with " << ival << std::endl;
    exit(1);
  }
  return;
}

void PHG4DetectorSubsystem::set_default_string_param(const std::string &name, const std::string &sval)
{
  if (default_string.find(name) == default_string.end())
  {
    default_string[name] = sval;
  }
  else
  {
    std::cout << "trying to overwrite default string " << name << " "
              << default_string[name] << " with " << sval << std::endl;
    exit(1);
  }
  return;
}

void PHG4DetectorSubsystem::InitializeParameters()
{
  set_default_int_param("absorberactive", 0);
  set_default_int_param("absorbertruth", 0);
  set_default_int_param("active", 0);
  set_default_int_param("blackhole", 0);
  set_default_int_param("supportactive", 0);

  SetDefaultParameters();  // call method from specific subsystem
  // now load those parameters to our params class
  for (std::map<const std::string, double>::const_iterator iter = default_double.begin(); iter != default_double.end(); ++iter)
  {
    params->set_double_param(iter->first, iter->second);
  }
  for (std::map<const std::string, int>::const_iterator iter = default_int.begin(); iter != default_int.end(); ++iter)
  {
    params->set_int_param(iter->first, iter->second);
  }
  for (std::map<const std::string, std::string>::const_iterator iter = default_string.begin(); iter != default_string.end(); ++iter)
  {
    params->set_string_param(iter->first, iter->second);
  }
}

int PHG4DetectorSubsystem::SaveParamsToDB()
{
  int iret = 0;
  if (paramscontainer)
  {
    iret = paramscontainer->WriteToDB();
  }
  else
  {
    iret = params->WriteToDB();
  }
  if (iret)
  {
    std::cout << "problem committing to DB" << std::endl;
  }
  return iret;
}

int PHG4DetectorSubsystem::ReadParamsFromDB(const std::string &name, const int issuper)
{
  int iret = 0;
  if (issuper)
  {
    iret = params->ReadFromDB(name, layer);
  }
  else
  {
    iret = params->ReadFromDB();
  }
  if (iret)
  {
    std::cout << "problem reading from DB" << std::endl;
  }
  return iret;
}

int PHG4DetectorSubsystem::SaveParamsToFile(const PHG4DetectorSubsystem::FILE_TYPE ftyp)
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
  if (paramscontainer)
  {
    iret = paramscontainer->WriteToFile(extension, calibfiledir);
  }
  else
  {
    iret = params->WriteToFile(extension, calibfiledir);
  }
  if (iret)
  {
    std::cout << "problem saving to " << extension << " file " << std::endl;
  }
  return iret;
}

int PHG4DetectorSubsystem::ReadParamsFromFile(const std::string &name, const PHG4DetectorSubsystem::FILE_TYPE ftyp, const int issuper)
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
  int iret = params->ReadFromFile(name, extension, layer, issuper, calibfiledir);
  if (iret)
  {
    std::cout << "problem reading from " << extension << " file " << std::endl;
  }
  return iret;
}

void PHG4DetectorSubsystem::SetActive(const int i)
{
  iparams["active"] = i;
}

void PHG4DetectorSubsystem::SetAbsorberActive(const int i)
{
  iparams["absorberactive"] = i;
}

void PHG4DetectorSubsystem::SetSupportActive(const int i)
{
  iparams["supportactive"] = i;
}

void PHG4DetectorSubsystem::BlackHole(const int i)
{
  iparams["blackhole"] = i;
}

void PHG4DetectorSubsystem::SetAbsorberTruth(const int i)
{
  iparams["absorbertruth"] = i;
}

int PHG4DetectorSubsystem::ReadParamsFromCDB(const std::string &domain)
{
  params->ReadFromCDB(domain);
  return 0;
}
