/**
 * \file TpcSpaceChargeMatrixInversion1D.cc
 * \brief performs 1D space charge distortion reconstruction using tracks
 * \author Xudong Yu <xyu3@bnl.gov>
 */

#include "TpcSpaceChargeMatrixInversion1D.h"
#include "TpcSpaceChargeMatrixContainer1D.h"
#include "TpcSpaceChargeReconstructionHelper.h"

#include <frog/FROG.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <memory>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

namespace
{
  // layer range
  static constexpr int m_layermin = 7;
  static constexpr int m_layermax = 55;

  // r range
  static constexpr float m_radiusmin = 20;
  static constexpr float m_radiusmax = 78;

}  // namespace

//_____________________________________________________________________
TpcSpaceChargeMatrixInversion1D::TpcSpaceChargeMatrixInversion1D(const std::string& name)
  : Fun4AllBase(name)
{
}

//_____________________________________________________________________
bool TpcSpaceChargeMatrixInversion1D::add_from_file(const std::string& shortfilename, const std::string& objectname)
{
  // get filename from frog
  FROG frog;
  const auto filename = frog.location(shortfilename);
  if (Verbosity())
  {
    std::cout<<"adding "<<filename<<std::endl;
  }

  // open TFile
  std::unique_ptr<TFile> inputfile(TFile::Open(filename));
  if (!inputfile)
  {
    std::cout << "TpcSpaceChargeMatrixInversion1D::add_from_file - could not open file " << filename << std::endl;
    return false;
  }

  // load object from input file
  std::unique_ptr<TpcSpaceChargeMatrixContainer> source_layer_negz(dynamic_cast<TpcSpaceChargeMatrixContainer*>(inputfile->Get((objectname+"_layer_negz").c_str())));
  std::unique_ptr<TpcSpaceChargeMatrixContainer> source_layer_posz(dynamic_cast<TpcSpaceChargeMatrixContainer*>(inputfile->Get((objectname+"_layer_posz").c_str())));
  std::unique_ptr<TpcSpaceChargeMatrixContainer> source_radius_negz(dynamic_cast<TpcSpaceChargeMatrixContainer*>(inputfile->Get((objectname+"_radius_negz").c_str())));
  std::unique_ptr<TpcSpaceChargeMatrixContainer> source_radius_posz(dynamic_cast<TpcSpaceChargeMatrixContainer*>(inputfile->Get((objectname+"_radius_posz").c_str())));

  if (!source_layer_negz || !source_layer_posz || !source_radius_negz || !source_radius_posz)
  {
    std::cout << "TpcSpaceChargeMatrixInversion1D::add_from_file - could not find object name " << objectname << " in file " << filename << std::endl;
    return false;
  }

  // add object
  return add(*source_layer_negz.get(),*source_layer_posz.get(),*source_radius_negz.get(),*source_radius_posz.get());
}

//_____________________________________________________________________
bool TpcSpaceChargeMatrixInversion1D::add(const TpcSpaceChargeMatrixContainer& source_layer_negz, const TpcSpaceChargeMatrixContainer& source_layer_posz, const TpcSpaceChargeMatrixContainer& source_radius_negz, const TpcSpaceChargeMatrixContainer& source_radius_posz)
{
  // check internal container, create if necessary
  if (!m_matrix_container_layer_negz)
  {
    m_matrix_container_layer_negz.reset(new TpcSpaceChargeMatrixContainer1D);

    // get grid dimensions from source
    int layerbins = 0;
    source_layer_negz.get_grid_dimensions(layerbins);
    if (Verbosity())
    {
      std::cout<<"negz layerbins = "<<layerbins<<std::endl;
    }

    // assign
    m_matrix_container_layer_negz->set_grid_dimensions(layerbins);
  }

  if (!m_matrix_container_layer_posz)
  {
    m_matrix_container_layer_posz.reset(new TpcSpaceChargeMatrixContainer1D);

    // get grid dimensions from source
    int layerbins = 0;
    source_layer_posz.get_grid_dimensions(layerbins);
    if (Verbosity())
    {
      std::cout<<"posz layerbins = "<<layerbins<<std::endl;
    }

    // assign
    m_matrix_container_layer_posz->set_grid_dimensions(layerbins);
  }

  if (!m_matrix_container_radius_negz)
  {
    m_matrix_container_radius_negz.reset(new TpcSpaceChargeMatrixContainer1D);

    // get grid dimensions from source
    int radiusbins = 0;
    source_radius_negz.get_grid_dimensions(radiusbins);
    if (Verbosity())
    {
      std::cout<<"negz radiusbins = "<<radiusbins<<std::endl;
    }

    // assign
    m_matrix_container_radius_negz->set_grid_dimensions(radiusbins);
  }

  if (!m_matrix_container_radius_posz)
  {
    m_matrix_container_radius_posz.reset(new TpcSpaceChargeMatrixContainer1D);

    // get grid dimensions from source
    int radiusbins = 0;
    source_radius_posz.get_grid_dimensions(radiusbins);
    if (Verbosity())
    {
      std::cout<<"posz radiusbins = "<<radiusbins<<std::endl;
    }

    // assign
    m_matrix_container_radius_posz->set_grid_dimensions(radiusbins);
  }

  // add content
  return m_matrix_container_layer_negz->add(source_layer_negz) && m_matrix_container_layer_posz->add(source_layer_posz) && m_matrix_container_radius_negz->add(source_radius_negz) && m_matrix_container_radius_posz->add(source_radius_posz);
}

//_____________________________________________________________________
void TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections()
{
  // get grid dimensions from matrix container
  int layerbins_negz = 0;
  int layerbins_posz = 0;
  m_matrix_container_layer_negz->get_grid_dimensions(layerbins_negz);
  m_matrix_container_layer_posz->get_grid_dimensions(layerbins_posz);

  int radiusbins_negz = 0;
  int radiusbins_posz = 0;
  m_matrix_container_radius_negz->get_grid_dimensions(radiusbins_negz);
  m_matrix_container_radius_posz->get_grid_dimensions(radiusbins_posz);

  // create output histograms
  std::unique_ptr<TH1> hentries_layer_negz(new TH1F("hentries_layer_negz", "hentries_layer_negz;layer;Entries", layerbins_negz, m_layermin, m_layermax));
  std::unique_ptr<TH1> hphi_layer_negz(new TH1F("hDistortionP_layer_negz", "hDistortionP_layer_negz;layer;#Delta#phi (rad)", layerbins_negz, m_layermin, m_layermax));
  std::unique_ptr<TH1> hz_layer_negz(new TH1F("hDistortionZ_layer_negz", "hDistortionZ_layer_negz;layer;#DeltaZ (cm)", layerbins_negz, m_layermin, m_layermax));
  std::unique_ptr<TH1> hr_layer_negz(new TH1F("hDistortionR_layer_negz", "hDistortionR_layer_negz;layer;#DeltaR (cm)", layerbins_negz, m_layermin, m_layermax));

  std::unique_ptr<TH1> hentries_layer_posz(new TH1F("hentries_layer_posz", "hentries_layer_posz;layer;Entries", layerbins_posz, m_layermin, m_layermax));
  std::unique_ptr<TH1> hphi_layer_posz(new TH1F("hDistortionP_layer_posz", "hDistortionP_layer_posz;layer;#Delta#phi (rad)", layerbins_posz, m_layermin, m_layermax));
  std::unique_ptr<TH1> hz_layer_posz(new TH1F("hDistortionZ_layer_posz", "hDistortionZ_layer_posz;layer;#DeltaZ (cm)", layerbins_posz, m_layermin, m_layermax));
  std::unique_ptr<TH1> hr_layer_posz(new TH1F("hDistortionR_layer_posz", "hDistortionR_layer_posz;layer;#DeltaR (cm)", layerbins_posz, m_layermin, m_layermax));

  std::unique_ptr<TH1> hentries_radius_negz(new TH1F("hentries_radius_negz", "hentries_radius_negz;radius;Entries", radiusbins_negz, m_radiusmin, m_radiusmax));
  std::unique_ptr<TH1> hphi_radius_negz(new TH1F("hDistortionP_radius_negz", "hDistortionP_radius_negz;radius;#Delta#phi (rad)", radiusbins_negz, m_radiusmin, m_radiusmax));
  std::unique_ptr<TH1> hz_radius_negz(new TH1F("hDistortionZ_radius_negz", "hDistortionZ_radius_negz;radius;#DeltaZ (cm)", radiusbins_negz, m_radiusmin, m_radiusmax));
  std::unique_ptr<TH1> hr_radius_negz(new TH1F("hDistortionR_radius_negz", "hDistortionR_radius_negz;radius;#DeltaR (cm)", radiusbins_negz, m_radiusmin, m_radiusmax));

  std::unique_ptr<TH1> hentries_radius_posz(new TH1F("hentries_radius_posz", "hentries_radius_posz;radius;Entries", radiusbins_posz, m_radiusmin, m_radiusmax));
  std::unique_ptr<TH1> hphi_radius_posz(new TH1F("hDistortionP_radius_posz", "hDistortionP_radius_posz;radius;#Delta#phi (rad)", radiusbins_posz, m_radiusmin, m_radiusmax));
  std::unique_ptr<TH1> hz_radius_posz(new TH1F("hDistortionZ_radius_posz", "hDistortionZ_radius_posz;radius;#DeltaZ (cm)", radiusbins_posz, m_radiusmin, m_radiusmax));
  std::unique_ptr<TH1> hr_radius_posz(new TH1F("hDistortionR_radius_posz", "hDistortionR_radius_posz;radius;#DeltaR (cm)", radiusbins_posz, m_radiusmin, m_radiusmax));

  // matrix convenience definition
  /* number of coordinates must match that of the matrix container */
  static constexpr int ncoord = 3;
  using matrix_t = Eigen::Matrix<float, ncoord, ncoord>;
  using column_t = Eigen::Matrix<float, ncoord, 1>;

  // loop over bins
  for (int ilayer = 0; ilayer < layerbins_negz; ++ilayer)
  {
    // get cell index
    const auto icell = m_matrix_container_layer_negz->get_cell_index(ilayer);

    // minimum number of entries per bin
    static constexpr int min_cluster_count = 10;
    const auto cell_entries = m_matrix_container_layer_negz->get_entries(icell);
    if (cell_entries < min_cluster_count)
    {
      continue;
    }

    // build eigen matrices from container
    matrix_t lhs;
    for (int i = 0; i < ncoord; ++i)
    {
      for (int j = 0; j < ncoord; ++j)
      {
        lhs(i, j) = m_matrix_container_layer_negz->get_lhs(icell, i, j);
        if ( Verbosity()>2 )
        {
          std::cout<<"layer negz lhs("<<i<<","<<j<<") = "<<lhs(i, j)<<std::endl;
        }
      }
    }

    column_t rhs;
    for (int i = 0; i < ncoord; ++i)
    {
      rhs(i) = m_matrix_container_layer_negz->get_rhs(icell, i);
    }

    if (Verbosity())
    {
      // print matrices and entries
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer negz inverting bin " << ilayer << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer negz entries: " << cell_entries << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer negz lhs: \n"
                << lhs << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer negz rhs: \n"
                << rhs << std::endl;
    }

    // calculate result using linear solving
    const auto cov = lhs.inverse();
    auto partialLu = lhs.partialPivLu();
    const auto result = partialLu.solve(rhs);

    // fill histograms
    hentries_layer_negz->SetBinContent(ilayer + 1, cell_entries);

    hphi_layer_negz->SetBinContent(ilayer + 1, result(0));
    hphi_layer_negz->SetBinError(ilayer + 1, std::sqrt(cov(0, 0)));

    hz_layer_negz->SetBinContent(ilayer + 1, result(1));
    hz_layer_negz->SetBinError(ilayer + 1, std::sqrt(cov(1, 1)));

    hr_layer_negz->SetBinContent(ilayer + 1, result(2));
    hr_layer_negz->SetBinError(ilayer + 1, std::sqrt(cov(2, 2)));

    if (Verbosity())
    {
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer negz drphi: " << result(0) << " +/- " << std::sqrt(cov(0, 0)) << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer negz dz: " << result(1) << " +/- " << std::sqrt(cov(1, 1)) << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer negz dr: " << result(2) << " +/- " << std::sqrt(cov(2, 2)) << std::endl;
      std::cout << std::endl;
    }
    
  }

  for (int ilayer = 0; ilayer < layerbins_posz; ++ilayer)
  {
    // get cell index
    const auto icell = m_matrix_container_layer_posz->get_cell_index(ilayer);

    // minimum number of entries per bin
    static constexpr int min_cluster_count = 2;
    const auto cell_entries = m_matrix_container_layer_posz->get_entries(icell);
    if (cell_entries < min_cluster_count)
    {
      continue;
    }

    // build eigen matrices from container
    matrix_t lhs;
    for (int i = 0; i < ncoord; ++i)
    {
      for (int j = 0; j < ncoord; ++j)
      {
        lhs(i, j) = m_matrix_container_layer_posz->get_lhs(icell, i, j);
        if ( Verbosity()>2 )
        {
          std::cout<<"posz lhs("<<i<<","<<j<<") = "<<lhs(i, j)<<std::endl;
        }
      }
    }

    column_t rhs;
    for (int i = 0; i < ncoord; ++i)
    {
      rhs(i) = m_matrix_container_layer_posz->get_rhs(icell, i);
    }

    if (Verbosity())
    {
      // print matrices and entries
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer posz inverting bin " << ilayer << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer posz entries: " << cell_entries << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer posz lhs: \n"
                << lhs << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer posz rhs: \n"
                << rhs << std::endl;
    }

    // calculate result using linear solving
    const auto cov = lhs.inverse();
    auto partialLu = lhs.partialPivLu();
    const auto result = partialLu.solve(rhs);

    // fill histograms
    hentries_layer_posz->SetBinContent(ilayer + 1, cell_entries);

    hphi_layer_posz->SetBinContent(ilayer + 1, result(0));
    hphi_layer_posz->SetBinError(ilayer + 1, std::sqrt(cov(0, 0)));

    hz_layer_posz->SetBinContent(ilayer + 1, result(1));
    hz_layer_posz->SetBinError(ilayer + 1, std::sqrt(cov(1, 1)));

    hr_layer_posz->SetBinContent(ilayer + 1, result(2));
    hr_layer_posz->SetBinError(ilayer + 1, std::sqrt(cov(2, 2)));

    if (Verbosity())
    {
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer posz drphi: " << result(0) << " +/- " << std::sqrt(cov(0, 0)) << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer posz dz: " << result(1) << " +/- " << std::sqrt(cov(1, 1)) << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - layer posz dr: " << result(2) << " +/- " << std::sqrt(cov(2, 2)) << std::endl;
      std::cout << std::endl;
    }
    
  }


  // loop over bins
  for (int iradius = 0; iradius < radiusbins_negz; ++iradius)
  {
    // get cell index
    const auto icell = m_matrix_container_radius_negz->get_cell_index(iradius);

    // minimum number of entries per bin
    static constexpr int min_cluster_count = 2;
    const auto cell_entries = m_matrix_container_radius_negz->get_entries(icell);
    if (cell_entries < min_cluster_count)
    {
      continue;
    }

    // build eigen matrices from container
    matrix_t lhs;
    for (int i = 0; i < ncoord; ++i)
    {
      for (int j = 0; j < ncoord; ++j)
      {
        lhs(i, j) = m_matrix_container_radius_negz->get_lhs(icell, i, j);
        if ( Verbosity()>2 )
        {
          std::cout<<"radius negz lhs("<<i<<","<<j<<") = "<<lhs(i, j)<<std::endl;
        }
      }
    }

    column_t rhs;
    for (int i = 0; i < ncoord; ++i)
    {
      rhs(i) = m_matrix_container_radius_negz->get_rhs(icell, i);
    }

    if (Verbosity())
    {
      // print matrices and entries
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius negz inverting bin " << iradius << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius negz entries: " << cell_entries << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius negz lhs: \n"
                << lhs << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius negz rhs: \n"
                << rhs << std::endl;
    }

    // calculate result using linear solving
    const auto cov = lhs.inverse();
    auto partialLu = lhs.partialPivLu();
    const auto result = partialLu.solve(rhs);

    // fill histograms
    hentries_radius_negz->SetBinContent(iradius + 1, cell_entries);

    hphi_radius_negz->SetBinContent(iradius + 1, result(0));
    hphi_radius_negz->SetBinError(iradius + 1, std::sqrt(cov(0, 0)));

    hz_radius_negz->SetBinContent(iradius + 1, result(1));
    hz_radius_negz->SetBinError(iradius + 1, std::sqrt(cov(1, 1)));

    hr_radius_negz->SetBinContent(iradius + 1, result(2));
    hr_radius_negz->SetBinError(iradius + 1, std::sqrt(cov(2, 2)));

    if (Verbosity())
    {
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius negz drphi: " << result(0) << " +/- " << std::sqrt(cov(0, 0)) << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius negz dz: " << result(1) << " +/- " << std::sqrt(cov(1, 1)) << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius negz dr: " << result(2) << " +/- " << std::sqrt(cov(2, 2)) << std::endl;
      std::cout << std::endl;
    }
    
  }

  for (int iradius = 0; iradius < radiusbins_posz; ++iradius)
  {
    // get cell index
    const auto icell = m_matrix_container_radius_posz->get_cell_index(iradius);

    // minimum number of entries per bin
    static constexpr int min_cluster_count = 2;
    const auto cell_entries = m_matrix_container_radius_posz->get_entries(icell);
    if (cell_entries < min_cluster_count)
    {
      continue;
    }

    // build eigen matrices from container
    matrix_t lhs;
    for (int i = 0; i < ncoord; ++i)
    {
      for (int j = 0; j < ncoord; ++j)
      {
        lhs(i, j) = m_matrix_container_radius_posz->get_lhs(icell, i, j);
        if ( Verbosity()>2 )
        {
          std::cout<<"posz lhs("<<i<<","<<j<<") = "<<lhs(i, j)<<std::endl;
        }
      }
    }

    column_t rhs;
    for (int i = 0; i < ncoord; ++i)
    {
      rhs(i) = m_matrix_container_radius_posz->get_rhs(icell, i);
    }

    if (Verbosity())
    {
      // print matrices and entries
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius posz inverting bin " << iradius << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius posz entries: " << cell_entries << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius posz lhs: \n"
                << lhs << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius posz rhs: \n"
                << rhs << std::endl;
    }

    // calculate result using linear solving
    const auto cov = lhs.inverse();
    auto partialLu = lhs.partialPivLu();
    const auto result = partialLu.solve(rhs);

    // fill histograms
    hentries_radius_posz->SetBinContent(iradius + 1, cell_entries);

    hphi_radius_posz->SetBinContent(iradius + 1, result(0));
    hphi_radius_posz->SetBinError(iradius + 1, std::sqrt(cov(0, 0)));

    hz_radius_posz->SetBinContent(iradius + 1, result(1));
    hz_radius_posz->SetBinError(iradius + 1, std::sqrt(cov(1, 1)));

    hr_radius_posz->SetBinContent(iradius + 1, result(2));
    hr_radius_posz->SetBinError(iradius + 1, std::sqrt(cov(2, 2)));

    if (Verbosity())
    {
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius posz drphi: " << result(0) << " +/- " << std::sqrt(cov(0, 0)) << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius posz dz: " << result(1) << " +/- " << std::sqrt(cov(1, 1)) << std::endl;
      std::cout << "TpcSpaceChargeMatrixInversion1D::calculate_distortion_corrections - radius posz dr: " << result(2) << " +/- " << std::sqrt(cov(2, 2)) << std::endl;
      std::cout << std::endl;
    }
    
  }

  // split histograms in two along z axis and write
  // also write histograms suitable for space charge reconstruction
  auto process_histogram = [](TH1* h_negz, TH1* h_posz, const TString& name)
  {
/*
    const auto& result = TpcSpaceChargeReconstructionHelper::split(h);
    std::unique_ptr<TH1> hneg(std::get<0>(result));
    std::unique_ptr<TH1> hpos(std::get<1>(result));

    return std::make_tuple(
        TpcSpaceChargeReconstructionHelper::add_guarding_bins(hneg.get(), name + "_negz"),
        TpcSpaceChargeReconstructionHelper::add_guarding_bins(hpos.get(), name + "_posz"));
*/

    return std::make_tuple(
        TpcSpaceChargeReconstructionHelper::add_guarding_bins(h_negz, name + "_negz"),
        TpcSpaceChargeReconstructionHelper::add_guarding_bins(h_posz, name + "_posz"));
  };

  // apply finishing transformations to histograms and save in container
  if (!m_dcc_average_layer)
  {
    m_dcc_average_layer.reset(new TpcDistortionCorrectionContainer);
  }

  std::tie(m_dcc_average_layer->m_hentries[0], m_dcc_average_layer->m_hentries[1]) = process_histogram(hentries_layer_negz.get(), hentries_layer_posz.get(), "hentries");
  std::tie(m_dcc_average_layer->m_hDRint[0], m_dcc_average_layer->m_hDRint[1]) = process_histogram(hr_layer_negz.get(), hr_layer_posz.get(), "hIntDistortionR");
  std::tie(m_dcc_average_layer->m_hDPint[0], m_dcc_average_layer->m_hDPint[1]) = process_histogram(hphi_layer_negz.get(), hphi_layer_posz.get(), "hIntDistortionP");
  std::tie(m_dcc_average_layer->m_hDZint[0], m_dcc_average_layer->m_hDZint[1]) = process_histogram(hz_layer_negz.get(), hz_layer_posz.get(), "hIntDistortionZ");

  std::tie(h_rdphi_layer_negz, h_rdphi_layer_posz) = process_histogram(transform_dphi_rdphi_layer(hphi_layer_negz.get(), "hIntDistortionRP"), transform_dphi_rdphi_layer(hphi_layer_posz.get(), "hIntDistortionRP"), "hIntDistortionRP");


  if (!m_dcc_average_radius)
  {
    m_dcc_average_radius.reset(new TpcDistortionCorrectionContainer);
  }

  std::tie(m_dcc_average_radius->m_hentries[0], m_dcc_average_radius->m_hentries[1]) = process_histogram(hentries_radius_negz.get(), hentries_radius_posz.get(), "hentries");
  std::tie(m_dcc_average_radius->m_hDRint[0], m_dcc_average_radius->m_hDRint[1]) = process_histogram(hr_radius_negz.get(), hr_radius_posz.get(), "hIntDistortionR");
  std::tie(m_dcc_average_radius->m_hDPint[0], m_dcc_average_radius->m_hDPint[1]) = process_histogram(hphi_radius_negz.get(), hphi_radius_posz.get(), "hIntDistortionP");
  std::tie(m_dcc_average_radius->m_hDZint[0], m_dcc_average_radius->m_hDZint[1]) = process_histogram(hz_radius_negz.get(), hz_radius_posz.get(), "hIntDistortionZ");

  std::tie(h_rdphi_radius_negz, h_rdphi_radius_posz) = process_histogram(transform_dphi_rdphi_radius(hphi_radius_negz.get(), "hIntDistortionRP"), transform_dphi_rdphi_radius(hphi_radius_posz.get(), "hIntDistortionRP"), "hIntDistortionRP");

}

//_____________________________________________________________________
void TpcSpaceChargeMatrixInversion1D::save_distortion_corrections(const std::string& filename)
{
  if (!m_dcc_average_layer)
  {
    std::cout << "TpcSpaceChargeMatrixInversion1D::save_distortion_corrections - invalid distortion correction container." << std::endl;
    return;
  }

  if (!m_dcc_average_radius)
  {
    std::cout << "TpcSpaceChargeMatrixInversion1D::save_distortion_corrections - invalid distortion correction container." << std::endl;
    return;
  }
  // save everything to root file
  std::cout << "TpcSpaceChargeMatrixInversion1D::save_distortions - writing histograms to " << filename << std::endl;
  std::unique_ptr<TFile> outputfile_layer(TFile::Open((filename+"_layer.root").c_str(), "RECREATE"));
  outputfile_layer->cd();

  for (const auto& h_list : {m_dcc_average_layer->m_hentries, m_dcc_average_layer->m_hDRint, m_dcc_average_layer->m_hDPint, m_dcc_average_layer->m_hDZint})
  {
    for (const auto& h : h_list)
    {
      if (h)
      {
        h->Write(h->GetName());
      }
    }
  }

  if (h_rdphi_layer_negz)
  {
    h_rdphi_layer_negz->Write(h_rdphi_layer_negz->GetName());
  }
  if (h_rdphi_layer_posz)
  {
    h_rdphi_layer_posz->Write(h_rdphi_layer_posz->GetName());
  }

  // close TFile
  outputfile_layer->Close();

  std::unique_ptr<TFile> outputfile_radius(TFile::Open((filename+"_radius.root").c_str(), "RECREATE"));
  outputfile_radius->cd();

  for (const auto& h_list : {m_dcc_average_radius->m_hentries, m_dcc_average_radius->m_hDRint, m_dcc_average_radius->m_hDPint, m_dcc_average_radius->m_hDZint})
  {
    for (const auto& h : h_list)
    {
      if (h)
      {
        h->Write(h->GetName());
      }
    }
  }

  if (h_rdphi_radius_negz)
  {
    h_rdphi_radius_negz->Write(h_rdphi_radius_negz->GetName());
  }
  if (h_rdphi_radius_posz)
  {
    h_rdphi_radius_posz->Write(h_rdphi_radius_posz->GetName());
  }

  // close TFile
  outputfile_radius->Close();
}

//___________________________________________________________________________
TH1* TpcSpaceChargeMatrixInversion1D::transform_dphi_rdphi_layer(const TH1* source, const TString& name)
{
  // create new histogram
  auto hout = new TH1F(name, name,
                       source->GetXaxis()->GetNbins(), source->GetXaxis()->GetXmin(), source->GetXaxis()->GetXmax());

  // update axis legend
  hout->GetXaxis()->SetTitle(source->GetXaxis()->GetTitle());
  //hout->GetYaxis()->SetTitle(source->GetYaxis()->GetTitle());
  hout->GetYaxis()->SetTitle("r#Delta#phi (cm)");

  // copy content
  const auto layerbins = source->GetXaxis()->GetNbins();

  loadTpcRadius();

  // fill center
  for (int ilayer = 0; ilayer < layerbins; ++ilayer)
  {
    int layer = source->GetXaxis()->GetBinCenter(ilayer + 1);
    float r = TpcRadiusMap[layer];

    hout->SetBinContent(ilayer + 1, r*source->GetBinContent(ilayer + 1));
    hout->SetBinError(ilayer + 1, r*source->GetBinError(ilayer + 1));
    if (Verbosity())
    {
      std::cout<<"TpcSpaceChargeMatrixInversion1D::transform_dphi_rdphi_layer - "<<source->GetName()<<" layer "<<layer<<" , dphi "<<source->GetBinContent(ilayer + 1)<<" -> r "<<r<<" , rdphi "<<hout->GetBinContent(ilayer + 1)<<std::endl;
    }
  }

  return hout;
}

//___________________________________________________________________________
TH1* TpcSpaceChargeMatrixInversion1D::transform_dphi_rdphi_radius(const TH1* source, const TString& name)
{
  // create new histogram
  auto hout = new TH1F(name, name,
                       source->GetXaxis()->GetNbins(), source->GetXaxis()->GetXmin(), source->GetXaxis()->GetXmax());

  // update axis legend
  hout->GetXaxis()->SetTitle(source->GetXaxis()->GetTitle());
  //hout->GetYaxis()->SetTitle(source->GetYaxis()->GetTitle());
  hout->GetYaxis()->SetTitle("r#Delta#phi (cm)");

  // copy content
  const auto radiusbins = source->GetXaxis()->GetNbins();

  // fill center
  for (int iradius = 0; iradius < radiusbins; ++iradius)
  {
    float r = source->GetXaxis()->GetBinCenter(iradius + 1);

    hout->SetBinContent(iradius + 1, r*source->GetBinContent(iradius + 1));
    hout->SetBinError(iradius + 1, r*source->GetBinError(iradius + 1));
    if (Verbosity())
    {
      std::cout<<"TpcSpaceChargeMatrixInversion1D::transform_dphi_rdphi_radius - "<<source->GetName()<<" radius "<<r<<" , dphi "<<source->GetBinContent(iradius + 1)<<" -> "<<" rdphi "<<hout->GetBinContent(iradius + 1)<<std::endl;
    }
  }

  return hout;
}
