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
  TFile *f = TFile::Open(url.c_str());
  if (!f)
  {
    std::cout << "could not open " << url
	      << " for domain " << domain << std::endl;
    gSystem->Exit(1);
  }
  PdbParameterMap *myparm = static_cast<PdbParameterMap *>(f->Get("PdbParameterMap"));
  if (!myparm)
  {
    std::cout << "could not get PdbParameterMap from " << url << std::endl;
    gSystem->Exit(1);
  }
  params.FillFrom(myparm);
  delete myparm;
  delete f;
  std::cout << "filling params from cdb " << domain << std::endl;
  params.Print();
}
