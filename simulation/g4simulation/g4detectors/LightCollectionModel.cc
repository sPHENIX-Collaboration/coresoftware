#include "LightCollectionModel.h"

#include <fun4all/Fun4AllServer.h>

#include <ffamodules/CDBInterface.h>

#include <phool/recoConsts.h>

#include <TAxis.h>  // for TAxis
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TObject.h>  // for TObject
#include <TSystem.h>

#include <cassert>
#include <iostream>

LightCollectionModel::~LightCollectionModel()
{
  delete data_grid_light_guide_efficiency;
  delete data_grid_fiber_trans;
}

void LightCollectionModel::load_data_from_CDB(
    const std::string &domain,
    const std::string &histogram_light_guide_model,
    const std::string &histogram_fiber_model)
{
  recoConsts *rc = recoConsts::instance();
  std::string url = CDBInterface::instance()->getUrl(domain);
  if (url.empty())
  {
    std::cout << "No calibration for domain " << domain << " for timestamp " << rc->get_uint64Flag("TIMESTAMP") << std::endl;
    gSystem->Exit(1);
  }
  TFile *fin = TFile::Open(url.c_str());
  if (!fin)
  {
    std::cout << "could not open " << url << std::endl;
    gSystem->Exit(1);
  }
  delete data_grid_light_guide_efficiency;
  data_grid_light_guide_efficiency = dynamic_cast<TH2 *>(fin->Get(histogram_light_guide_model.c_str()));
  assert(data_grid_light_guide_efficiency);
  data_grid_light_guide_efficiency->SetDirectory(nullptr);
  delete data_grid_fiber_trans;
  data_grid_fiber_trans = dynamic_cast<TH1 *>(fin->Get(histogram_fiber_model.c_str()));
  assert(data_grid_fiber_trans);
  data_grid_fiber_trans->SetDirectory(nullptr);
  delete fin;
}

void LightCollectionModel::load_data_file(
    const std::string &input_file,
    const std::string &histogram_light_guide_model,
    const std::string &histogram_fiber_model)
{
  TFile *fin = TFile::Open(input_file.c_str());

  assert(fin);
  assert(fin->IsOpen());

  delete data_grid_light_guide_efficiency;
  data_grid_light_guide_efficiency = dynamic_cast<TH2 *>(fin->Get(histogram_light_guide_model.c_str()));
  assert(data_grid_light_guide_efficiency);
  data_grid_light_guide_efficiency->SetDirectory(nullptr);

  delete data_grid_fiber_trans;
  data_grid_fiber_trans = dynamic_cast<TH1 *>(fin->Get(histogram_fiber_model.c_str()));
  assert(data_grid_fiber_trans);
  data_grid_fiber_trans->SetDirectory(nullptr);

  delete fin;
}

double LightCollectionModel::get_light_guide_efficiency(const double x_fraction, const double y_fraction)
{
  assert(data_grid_light_guide_efficiency);
  assert(x_fraction >= 0);
  assert(x_fraction <= 1);
  assert(y_fraction >= 0);
  assert(y_fraction <= 1);

  const double eff = data_grid_light_guide_efficiency->Interpolate(x_fraction,
                                                                   y_fraction);

  if (!data_grid_light_guide_efficiency_verify)
  {
    data_grid_light_guide_efficiency_verify = new TH2F("data_grid_light_guide_efficiency_verify",
                                                       "light collection efficiency as used in LightCollectionModel;x positio fraction;y position fraction",  //
                                                       100, 0., 1., 100, 0., 1.);
    Fun4AllServer::instance()->registerHisto(data_grid_light_guide_efficiency_verify);
  }

  data_grid_light_guide_efficiency_verify->SetBinContent(                        //
      data_grid_light_guide_efficiency_verify->GetXaxis()->FindBin(x_fraction),  //
      data_grid_light_guide_efficiency_verify->GetYaxis()->FindBin(y_fraction),  //
      eff                                                                        //
  );

  return eff;
}

double LightCollectionModel::get_fiber_transmission(const double z_distance)
{
  assert(data_grid_fiber_trans);

  const double eff = data_grid_fiber_trans->Interpolate(z_distance);
  if (!data_grid_fiber_trans_verify)
  {
    data_grid_fiber_trans_verify = new TH1F("data_grid_fiber_trans",
                                            "SCSF-78 Fiber Transmission as used in LightCollectionModel;position in fiber (cm);Effective transmission",
                                            100, -15, 15);
    Fun4AllServer::instance()->registerHisto(data_grid_fiber_trans_verify);
  }

  data_grid_fiber_trans_verify->SetBinContent(                        //
      data_grid_fiber_trans_verify->GetXaxis()->FindBin(z_distance),  //
      eff                                                             //
  );

  return eff;
}
