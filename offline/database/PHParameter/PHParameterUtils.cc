#include "PHParameterUtils.h"
#include "PHParameters.h"

#include <pdbcalbase/PdbParameterMap.h>

#include <ffamodules/CDBInterface.h>

#include <TFile.h>
#include <TSystem.h>

#include <iostream>

void PHParameterUtils::FillPHParametersFromCDB(PHParameters &params, const std::string &domain)
{
  std::string url =  CDBInterface::instance()->getUrl(domain);
  params.ReadFromCDBFile(url);
  return;
}
