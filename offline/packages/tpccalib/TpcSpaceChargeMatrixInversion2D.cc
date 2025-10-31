/**
 * \file TpcSpaceChargeMatrixInversion2D.cc
 * \brief performs 2D space charge distortion reconstruction using tracks
 * \author Xudong Yu <xyu3@bnl.gov>
 */

#include "TpcSpaceChargeMatrixInversion2D.h"
#include "TpcSpaceChargeMatrixContainer2D.h"
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

namespace
{
  // phi range
  static constexpr float m_phimin = 0;
  static constexpr float m_phimax = 2. * M_PI;

  // TODO: could try to get the r and z range from TPC geometry
  // r range
  static constexpr float m_rmin = 20;
  static constexpr float m_rmax = 78;

  // z range
  static constexpr float m_zmin = -105.5;
  static constexpr float m_zmax = 105.5;
}  // namespace

//_____________________________________________________________________
TpcSpaceChargeMatrixInversion2D::TpcSpaceChargeMatrixInversion2D(const std::string& name)
  : Fun4AllBase(name)
{
}

//_____________________________________________________________________
void TpcSpaceChargeMatrixInversion2D::load_cm_distortion_corrections(const std::string& filename)
{
  std::cout << "TpcSpaceChargeMatrixInversion2D::load_cm_distortion_corrections - loading " << filename << std::endl;

  // open TFile
  auto distortion_tfile = TFile::Open(filename.c_str());
  if (!distortion_tfile)
  {
    std::cout << "TpcSpaceChargeMatrixInversion2D::load_cm_distortion_corrections - cannot open " << filename << std::endl;
    exit(1);
  }

  // make sure container exists
  if (!m_dcc_cm)
  {
    m_dcc_cm.reset(new TpcDistortionCorrectionContainer);
  }

  const std::array<const std::string, 2> extension = {{"_negz", "_posz"}};
  for (int j = 0; j < 2; ++j)
  {
    m_dcc_cm->m_hDPint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionP") + extension[j]).c_str()));
    assert(m_dcc_cm->m_hDPint[j]);
    m_dcc_cm->m_hDRint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionR") + extension[j]).c_str()));
    assert(m_dcc_cm->m_hDRint[j]);
    m_dcc_cm->m_hDZint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionZ") + extension[j]).c_str()));
    assert(m_dcc_cm->m_hDZint[j]);
  }

  // check histogram dimension. We expect 2D histograms
  m_dcc_cm->m_dimensions = m_dcc_cm->m_hDPint[0]->GetDimension();
  assert(m_dcc_cm->m_dimensions == 2);
}

//_____________________________________________________________________
void TpcSpaceChargeMatrixInversion2D::load_average_distortion_corrections(const std::string& filename)
{
  std::cout << "TpcSpaceChargeMatrixInversion2D::load_average_distortion_corrections - loading " << filename << std::endl;

  // open TFile
  auto distortion_tfile = TFile::Open(filename.c_str());
  if (!distortion_tfile)
  {
    std::cout << "TpcSpaceChargeMatrixInversion2D::load_average_distortion_corrections - cannot open " << filename << std::endl;
    exit(1);
  }

  // make sure container exists
  if (!m_dcc_average)
  {
    m_dcc_average.reset(new TpcDistortionCorrectionContainer);
  }

  const std::array<const std::string, 2> extension = {{"_negz", "_posz"}};
  for (int j = 0; j < 2; ++j)
  {
    m_dcc_average->m_hDPint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionP") + extension[j]).c_str()));
    assert(m_dcc_average->m_hDPint[j]);
    m_dcc_average->m_hDRint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionR") + extension[j]).c_str()));
    assert(m_dcc_average->m_hDRint[j]);
    m_dcc_average->m_hDZint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionZ") + extension[j]).c_str()));
    assert(m_dcc_average->m_hDZint[j]);
  }

  // check histogram dimension. We expect 2D histograms
  m_dcc_average->m_dimensions = m_dcc_average->m_hDPint[0]->GetDimension();
  assert(m_dcc_average->m_dimensions == 3);
}

//_____________________________________________________________________
bool TpcSpaceChargeMatrixInversion2D::add_from_file(const std::string& shortfilename, const std::string& objectname)
{
  // get filename from frog
  FROG frog;
  const auto filename = frog.location(shortfilename);

  // open TFile
  std::unique_ptr<TFile> inputfile(TFile::Open(filename));
  if (!inputfile)
  {
    std::cout << "TpcSpaceChargeMatrixInversion2D::add_from_file - could not open file " << filename << std::endl;
    return false;
  }

  // load object from input file
  std::unique_ptr<TpcSpaceChargeMatrixContainer> source(dynamic_cast<TpcSpaceChargeMatrixContainer*>(inputfile->Get(objectname.c_str())));
  if (!source)
  {
    std::cout << "TpcSpaceChargeMatrixInversion2D::add_from_file - could not find object name " << objectname << " in file " << filename << std::endl;
    return false;
  }

  // add object
  return add(*source.get());
}

//_____________________________________________________________________
bool TpcSpaceChargeMatrixInversion2D::add(const TpcSpaceChargeMatrixContainer& source)
{
  // check internal container, create if necessary
  if (!m_matrix_container)
  {
    m_matrix_container.reset(new TpcSpaceChargeMatrixContainer2D);

    // get grid dimensions from source
    int pbins = 0;
    int rbins = 0;
    int zbins = 0;
    source.get_grid_dimensions(pbins, rbins, zbins);

    // assign
    m_matrix_container->set_grid_dimensions(pbins, rbins, zbins);
  }

  // add content
  return m_matrix_container->add(source);
}

//_____________________________________________________________________
void TpcSpaceChargeMatrixInversion2D::calculate_distortion_corrections()
{
  // get grid dimensions from matrix container
  int pbins = 0;
  int rbins = 0;
  int zbins = 0;
  m_matrix_container->get_grid_dimensions(pbins, rbins, zbins);

  // create output histograms
  std::unique_ptr<TH2> hentries(new TH2F("hentries_rec", "hentries_rec", rbins, m_rmin, m_rmax, zbins, m_zmin, m_zmax));
  std::unique_ptr<TH2> hphi(new TH2F("hDistortionP_rec", "hDistortionP_rec", rbins, m_rmin, m_rmax, zbins, m_zmin, m_zmax));
  std::unique_ptr<TH2> hz(new TH2F("hDistortionZ_rec", "hDistortionZ_rec", rbins, m_rmin, m_rmax, zbins, m_zmin, m_zmax));
  std::unique_ptr<TH2> hr(new TH2F("hDistortionR_rec", "hDistortionR_rec", rbins, m_rmin, m_rmax, zbins, m_zmin, m_zmax));

  // set axis labels
  for (const auto& h : {hentries.get(), hphi.get(), hz.get(), hr.get()})
  {
    h->GetXaxis()->SetTitle("r (cm)");
    h->GetYaxis()->SetTitle("z (cm)");
  }

  // matrix convenience definition
  /* number of coordinates must match that of the matrix container */
  static constexpr int ncoord = 3;
  using matrix_t = Eigen::Matrix<float, ncoord, ncoord>;
  using column_t = Eigen::Matrix<float, ncoord, 1>;

  // loop over bins
  for (int ir = 0; ir < rbins; ++ir)
  {
    for (int iz = 0; iz < zbins; ++iz)
    {
      // get cell index
      const auto icell = m_matrix_container->get_cell_index(ir, iz);

      // minimum number of entries per bin
      //static constexpr int min_cluster_count = 50;
      const auto cell_entries = m_matrix_container->get_entries(icell);
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
          lhs(i, j) = m_matrix_container->get_lhs(icell, i, j);
        }
      }

      column_t rhs;
      for (int i = 0; i < ncoord; ++i)
      {
        rhs(i) = m_matrix_container->get_rhs(icell, i);
      }

      if (Verbosity())
      {
        // print matrices and entries
        std::cout << "TpcSpaceChargeMatrixInversion2D::calculate_distortion_corrections - inverting bin " << iz << ", " << ir << std::endl;
        std::cout << "TpcSpaceChargeMatrixInversion2D::calculate_distortion_corrections - entries: " << cell_entries << std::endl;
        std::cout << "TpcSpaceChargeMatrixInversion2D::calculate_distortion_corrections - lhs: \n"
                  << lhs << std::endl;
        std::cout << "TpcSpaceChargeMatrixInversion2D::calculate_distortion_corrections - rhs: \n"
                  << rhs << std::endl;
      }

      // calculate result using linear solving
      const auto cov = lhs.inverse();
      auto partialLu = lhs.partialPivLu();
      const auto result = partialLu.solve(rhs);

      // fill histograms
      hentries->SetBinContent(ir + 1, iz + 1, cell_entries);

      hphi->SetBinContent(ir + 1, iz + 1, result(0));
      hphi->SetBinError(ir + 1, iz + 1, std::sqrt(cov(0, 0)));

      hz->SetBinContent(ir + 1, iz + 1, result(1));
      hz->SetBinError(ir + 1, iz + 1, std::sqrt(cov(1, 1)));

      hr->SetBinContent(ir + 1, iz + 1, result(2));
      hr->SetBinError(ir + 1, iz + 1, std::sqrt(cov(2, 2)));

      if (Verbosity())
      {
        std::cout << "TpcSpaceChargeMatrixInversion2D::calculate_distortion_corrections - drphi: " << result(0) << " +/- " << std::sqrt(cov(0, 0)) << std::endl;
        std::cout << "TpcSpaceChargeMatrixInversion2D::calculate_distortion_corrections - dz: " << result(1) << " +/- " << std::sqrt(cov(1, 1)) << std::endl;
        std::cout << "TpcSpaceChargeMatrixInversion2D::calculate_distortion_corrections - dr: " << result(2) << " +/- " << std::sqrt(cov(2, 2)) << std::endl;
        std::cout << std::endl;
      }
    }
  }

  // split histograms in two along z axis and write
  // also write histograms suitable for space charge reconstruction
  auto process_histogram = [](TH2* h, const TString& name, const int bins, const float min, const float max)
  {
    const auto& h3 = TpcSpaceChargeReconstructionHelper::expand(h, bins, min, max);
    const auto& result = TpcSpaceChargeReconstructionHelper::split(h3);
    std::unique_ptr<TH3> hneg(std::get<0>(result));
    std::unique_ptr<TH3> hpos(std::get<1>(result));

    return std::make_tuple(
        TpcSpaceChargeReconstructionHelper::add_guarding_bins(hneg.get(), name + "_negz"),
        TpcSpaceChargeReconstructionHelper::add_guarding_bins(hpos.get(), name + "_posz"));
  };

  // apply finishing transformations to histograms and save in container
  if (!m_dcc_average)
  {
    m_dcc_average.reset(new TpcDistortionCorrectionContainer);
  }

  std::tie(m_dcc_average->m_hentries[0], m_dcc_average->m_hentries[1]) = process_histogram(hentries.get(), "hentries", pbins, m_phimin, m_phimax);
  std::tie(m_dcc_average->m_hDRint[0], m_dcc_average->m_hDRint[1]) = process_histogram(hr.get(), "hIntDistortionR", pbins, m_phimin, m_phimax);
  std::tie(m_dcc_average->m_hDPint[0], m_dcc_average->m_hDPint[1]) = process_histogram(hphi.get(), "hIntDistortionP", pbins, m_phimin, m_phimax);
  std::tie(m_dcc_average->m_hDZint[0], m_dcc_average->m_hDZint[1]) = process_histogram(hz.get(), "hIntDistortionZ", pbins, m_phimin, m_phimax);
}

//_____________________________________________________________________
void TpcSpaceChargeMatrixInversion2D::extrapolate_distortion_corrections()
{
  if (!m_dcc_average)
  {
    std::cout << "TpcSpaceChargeMatrixInversion2D::extrapolate_distortion_corrections - invalid distortion correction container." << std::endl;
    return;
  }

  std::string name[2] = {"_negz", "_posz"};

  // handle sides independently
  for (int i = 0; i < 2; ++i)
  {
    if (!m_dcc_average->m_hDRint[i])
    {
      std::cout << "TpcSpaceChargeMatrixInversion2D::extrapolate_distortion_corrections - invalid histograms m_hDRint" << std::endl;
      continue;
    }

    // create relevant masks. They are used to define the cells in which distortion correction interpolation must be performed
    // this is a mask matching TPOT acceptance
    std::unique_ptr<TH3> hmask(static_cast<TH3*>(m_dcc_average->m_hDRint[i]->Clone("hmask")));
    TpcSpaceChargeReconstructionHelper::create_tpot_mask(hmask.get());

    // this is a mask matching TPOT acceptance, with small z interpolation between Micromegas modules
    std::unique_ptr<TH3> hmask_extrap_z(static_cast<TH3*>(hmask->Clone("hmask_extrap_z")));
    TpcSpaceChargeReconstructionHelper::extrapolate_z(hmask_extrap_z.get(), hmask.get());

    /*
     * this is a mask matching TPOT acceptance, with small z interpolation between Micromegas modules,
     * copied from sector to sector (no scaling)
     */
    std::unique_ptr<TH3> hmask_extrap_p(static_cast<TH3*>(hmask_extrap_z->Clone("hmask_extrap_p")));
    TpcSpaceChargeReconstructionHelper::extrapolate_phi1(hmask_extrap_p.get(), nullptr, hmask_extrap_z.get());

    /*
     * this is a mask matching TPOT acceptance along z, extrapolated over phi over 2pi.
     * Still empty are the bins from readout plane to the outermost micromegas module
     */
    std::unique_ptr<TH3> hmask_extrap_p2(static_cast<TH3*>(hmask_extrap_p->Clone("hmask_extrap_p2")));
    TpcSpaceChargeReconstructionHelper::extrapolate_phi2(hmask_extrap_p2.get(), hmask_extrap_p.get());

    /*
     * decide on which side of the histograms the pad plane is set.
     * this is used to guide the last z interpolation from readout plane to the outermost micromegas module
     */
    const uint8_t side = (i == 0) ? TpcSpaceChargeReconstructionHelper::Side_negative : TpcSpaceChargeReconstructionHelper::Side_positive;

    // labmda function to process a given histogram, and return the updated one
    auto process_histogram = [&hmask, &hmask_extrap_z, &hmask_extrap_p, &hmask_extrap_p2, side](TH3* h, TH2* h_cm)
    {
      // perform z extrapolation
      TpcSpaceChargeReconstructionHelper::extrapolate_z(h, hmask.get());

      // perform first phi extrapolation, sector to sector, with normalization from central membrane
      TpcSpaceChargeReconstructionHelper::extrapolate_phi1(h, h_cm, hmask_extrap_z.get());

      // perform second phi extrapolation, between sectors, using masks
      TpcSpaceChargeReconstructionHelper::extrapolate_phi2(h, hmask_extrap_p.get());

      // perform second z interpolation from readout plane to outermost micromegas module using masks
      TpcSpaceChargeReconstructionHelper::extrapolate_z2(h, hmask_extrap_p2.get(), side);
    };

    if (!m_dcc_cm)
    {
      m_dcc_cm.reset(new TpcDistortionCorrectionContainer);
    }

    process_histogram(static_cast<TH3*>(m_dcc_average->m_hDRint[i]), static_cast<TH2*>(m_dcc_cm->m_hDRint[i]));
    process_histogram(static_cast<TH3*>(m_dcc_average->m_hDPint[i]), static_cast<TH2*>(m_dcc_cm->m_hDPint[i]));
    process_histogram(static_cast<TH3*>(m_dcc_average->m_hDZint[i]), static_cast<TH2*>(m_dcc_cm->m_hDZint[i]));
    hDRPint[i] = transform_dphi_rdphi(static_cast<TH3*>(m_dcc_average->m_hDPint[i]), "hIntDistortionRP" + name[i]);
  }
}

//_____________________________________________________________________
void TpcSpaceChargeMatrixInversion2D::save_distortion_corrections(const std::string& filename)
{
  if (!m_dcc_average)
  {
    std::cout << "TpcSpaceChargeMatrixInversion2D::save_distortion_corrections - invalid distortion correction container." << std::endl;
    return;
  }

  // save everything to root file
  std::cout << "TpcSpaceChargeMatrixInversion2D::save_distortions - writing histograms to " << filename << std::endl;
  std::unique_ptr<TFile> outputfile(TFile::Open(filename.c_str(), "RECREATE"));
  outputfile->cd();

  for (const auto& h_list : {m_dcc_average->m_hentries, m_dcc_average->m_hDRint, m_dcc_average->m_hDPint, m_dcc_average->m_hDZint})
  {
    for (const auto& h : h_list)
    {
      if (h)
      {
        h->Write(h->GetName());
      }
    }
  }

  for (const auto& h : hDRPint)
  {
    if (h)
    {
      h->Write(h->GetName());
    }
  }

  // close TFile
  outputfile->Close();
}

//___________________________________________________________________________
TH3* TpcSpaceChargeMatrixInversion2D::transform_dphi_rdphi(const TH3* source, const TString& name)
{
  // create new histogram
  auto hout = new TH3F(name, name,
                       source->GetXaxis()->GetNbins(), source->GetXaxis()->GetXmin(), source->GetXaxis()->GetXmax(),
		       source->GetYaxis()->GetNbins(), source->GetYaxis()->GetXmin(), source->GetYaxis()->GetXmax(),
		       source->GetZaxis()->GetNbins(), source->GetZaxis()->GetXmin(), source->GetZaxis()->GetXmax());

  // update axis legend
  hout->GetXaxis()->SetTitle(source->GetXaxis()->GetTitle());
  hout->GetYaxis()->SetTitle(source->GetYaxis()->GetTitle());
  hout->GetZaxis()->SetTitle(source->GetZaxis()->GetTitle());
  hout->GetZaxis()->SetTitle("r#Delta#phi (cm)");

  // copy content
  const auto pbins = source->GetXaxis()->GetNbins();
  const auto rbins = source->GetYaxis()->GetNbins();
  const auto zbins = source->GetZaxis()->GetNbins();

  // fill center
  for (int ip = 0; ip < pbins; ++ip)
  {
    for (int ir = 0; ir < rbins; ++ir)
    {
      float r = source->GetYaxis()->GetBinCenter(ir + 1);
      for (int iz = 0; iz < zbins; ++iz)
      {
        hout->SetBinContent(ip + 1, ir + 1, iz + 1, r*source->GetBinContent(ip + 1, ir + 1, iz + 1));
        hout->SetBinError(ip + 1, ir + 1, iz + 1, r*source->GetBinError(ip + 1, ir + 1, iz + 1));
        if (Verbosity())
        {
          std::cout<<"TpcSpaceChargeMatrixInversion2D::transform_dphi_rdphi - "<<source->GetName()
		  <<" p "<<source->GetXaxis()->GetBinCenter(ip + 1)
		  <<" r "<<source->GetYaxis()->GetBinCenter(ir + 1)
		  <<" z "<<source->GetZaxis()->GetBinCenter(iz + 1)
		  <<" dphi "<<source->GetBinContent(ip + 1, ir + 1, iz + 1)
		  <<" -> "
		  <<" rdphi "<<hout->GetBinContent(ip + 1, ir + 1, iz + 1)
		  <<std::endl;
        }
      }
    } 
  }

  return hout;
}
