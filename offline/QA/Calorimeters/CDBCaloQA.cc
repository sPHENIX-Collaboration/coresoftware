#include "CDBCaloQA.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/TriggerAnalyzer.h>

#include <cdbobjects/CDBTTree.h>

#include <centrality/CentralityInfo.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <zdcinfo/Zdcinfo.h>

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/RunHeader.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/RunnumberRange.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TString.h>

#include <CLHEP/Vector/ThreeVector.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <tuple>
#include <vector>

namespace
{
  const std::array<std::string, CDBCaloQA::kNLayer> layername{"CEMC", "HCALIN", "HCALOUT"};
  const std::array<std::string, CDBCaloQA::kNLayer> layertag{"cemc", "ihcal", "ohcal"};
  const std::array<int, CDBCaloQA::kNLayer> netabins{96, 24, 24};
  const std::array<std::string, CDBCaloQA::kNTrig> trigtag{"MB", "Photon", "Jet"};
  const std::array<double, CDBCaloQA::kNLayer> ebinw{0.005, 0.002, 0.010};
  const std::array<double, CDBCaloQA::kNLayer> elo{-0.2, -0.15, -0.3};
  const std::array<double, CDBCaloQA::kNLayer> ehi{0.3, 0.2, 0.5};
}  // namespace

CDBCaloQA::CDBCaloQA(const std::string &name, const std::string &outfile)
  : SubsysReco(name)
  , m_outfile(outfile)
{
}

CDBCaloQA::~CDBCaloQA()
{
  delete m_hm;
  delete m_trigana;
  delete m_zdccalib;
}

int CDBCaloQA::InitRun(PHCompositeNode *topNode)
{
  RunHeader *runhdr = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (runhdr)
  {
    m_runnumber = runhdr->get_RunNumber();
  }
  if (m_runnumber >= RunnumberRange::RUN2AUAU_FIRST && m_runnumber <= RunnumberRange::RUN2AUAU_LAST)
  {
    m_species = "AuAu";
  }
  else if (m_runnumber >= RunnumberRange::RUN3AUAU_FIRST && m_runnumber <= RunnumberRange::RUN3AUAU_LAST)
  {
    m_species = "AuAu";
  }
  else if (m_runnumber >= RunnumberRange::RUN3OO_FIRST && m_runnumber <= RunnumberRange::RUN3OO_LAST)
  {
    m_species = "OO";
  }
  else
  {
    m_species = "pp";
  }
  if (!m_species_user.empty())
  {
    m_species = m_species_user;
  }
  m_mbbit = (m_runnumber >= RunnumberRange::RUN2PP_FIRST && m_runnumber <= RunnumberRange::RUN2PP_LAST) ? 10 : 14;
  if (m_mbbit_user >= 0)
  {
    m_mbbit = m_mbbit_user;
  }
  std::cout << Name() << ": run " << m_runnumber << ", species " << m_species << ", MB trigger bit " << m_mbbit << std::endl;

  if (!m_booked)
  {
    m_trigana = new TriggerAnalyzer();
    m_booked = true;
    if (m_select)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
    createHistos();
  }
  if (m_select)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  if (!m_mc && !m_zdccalib)
  {
    const std::string url = CDBInterface::instance()->getUrl("data_driven_zdc_calib");
    if (!url.empty())
    {
      m_zdccalib = new CDBTTree(url);
      m_zdccalib->LoadCalibrations();
    }
  }
  h_runnumber->SetBinContent(1, m_runnumber);
  return Fun4AllReturnCodes::EVENT_OK;
}

TH1 *CDBCaloQA::book1(const std::string &name, const std::string &title, int nb, double lo, double hi)
{
  TH1 *h = new TH1D(("h_" + Name() + "_" + name).c_str(), title.c_str(), nb, lo, hi);
  m_hm->registerHisto(h);
  return h;
}

TH2 *CDBCaloQA::book2(const std::string &name, const std::string &title, int nx, double xlo, double xhi, int ny, double ylo, double yhi)
{
  TH2 *h = new TH2D(("h2_" + Name() + "_" + name).c_str(), title.c_str(), nx, xlo, xhi, ny, ylo, yhi);
  m_hm->registerHisto(h);
  return h;
}

void CDBCaloQA::createHistos()
{
  m_hm = new Fun4AllHistoManager(Name() + "_HISTOS");
  m_hm->setOutfileName(m_outfile);

  float qmax = 25.;
  float qbinw = 0.01;
  float esummax = 100.;
  if (m_species == "OO")
  {
    qmax = 50.;
    esummax = 400.;
  }
  else if (m_species == "AuAu")
  {
    qmax = 200.;
    esummax = 2000.;
  }
  const float qlo = -1.;
  const int nq = std::lround((qmax - qlo) / qbinw);

  h_evt = book1("evt", ";;events", 4, 0, 4);
  h_evt->GetXaxis()->SetBinLabel(1, "all");
  for (int t = 0; t < kNTrig; t++)
  {
    h_evt->GetXaxis()->SetBinLabel(t + 2, trigtag[t].c_str());
  }
  h_runnumber = book1("runnumber", ";;run number", 1, 0, 1);
  h_vtxz = book1("vtxz", "GlobalVertex z, all events;z (cm);events", 200, -100, 100);
  h_noise_nevt = book1("noise_nevt", "events for noise: MB, |v_{z}| < 100 cm;;events", 1, 0, 1);

  h2_mbd_q_pmt = book2("mbd_q_pmt", "MBD PMT charge, valid TDC, MB;PMT;q", 128, 0, 128, nq, qlo, qmax);
  h2_mbd_t_pmt = book2("mbd_t_pmt", "MBD PMT time q>0, MB;PMT;t (ns)", 128, 0, 128, 500, -25, 25);
  h_mbd_qavg_south = book1("mbd_qavg_south", "MBD south average PMT charge, MB;<q>;events", 500, 0, qmax);
  h_mbd_qavg_north = book1("mbd_qavg_north", "MBD north average PMT charge, MB;<q>;events", 500, 0, qmax);
  h_mbd_npmt_south = book1("mbd_npmt_south", "MBD south PMTs q>0.4, MB;N_{PMT};events", 65, -0.5, 64.5);
  h_mbd_npmt_north = book1("mbd_npmt_north", "MBD north PMTs q>0.4, MB;N_{PMT};events", 65, -0.5, 64.5);
  h_mbd_zvtx = book1("mbd_zvtx", "MBD z-vertex, MB;z (cm);events", 300, -150, 150);
  h_mbd_t0 = book1("mbd_t0", "MBD t0, MB;t0 (ns);events", 200, -20, 20);

  float zdcmax = 1000.;
  float sepdqmax = 15.;
  float sepdqbinw = 0.05;
  float sepdsummax = 300.;
  if (m_species == "OO")
  {
    zdcmax = 2500.;
    sepdqmax = 30.;
    sepdsummax = 1000.;
  }
  else if (m_species == "AuAu")
  {
    zdcmax = 6000.;
    sepdqmax = 50.;
    sepdqbinw = 0.1;
    sepdsummax = 5000.;
  }
  const int nzdc = std::lround((zdcmax + 50.) / 2.5);
  const std::string zu = "GeV";
  h_zdc_e_south = book1("zdc_e_south", ("ZDC south energy, ZDC trigger;E (" + zu + ");events").c_str(), nzdc, -50., zdcmax);
  h_zdc_e_north = book1("zdc_e_north", ("ZDC north energy, ZDC trigger;E (" + zu + ");events").c_str(), nzdc, -50., zdcmax);
  for (int a = 0; a < 2; a++)
  {
    const std::string arm = a ? "north" : "south";
    for (int m = 0; m < 3; m++)
    {
      h_zdc_mod_e[3 * a + m] = book1(Form("zdc_mod_e_%s_%d", arm.c_str(), m + 1),
                                     Form("ZDC %s compartment %d energy, ZDC trigger;E (%s);events", arm.c_str(), m + 1, zu.c_str()), nzdc, -50., zdcmax);
    }
  }
  h_zdc_e_south_mb = book1("zdc_e_south_mb", ("ZDC south energy, MB;E (" + zu + ");events").c_str(), std::lround(zdcmax), 0, zdcmax);
  h_zdc_e_north_mb = book1("zdc_e_north_mb", ("ZDC north energy, MB;E (" + zu + ");events").c_str(), std::lround(zdcmax), 0, zdcmax);
  const double ohmax = (m_species == "pp") ? 20. : (m_species == "OO" ? 40. : 400.);
  const std::string zsel = (m_species == "pp") ? "MB" : "MB |z| < 10 cm";
  for (int a = 0; !m_mc && a < 2; a++)
  {
    const std::string arm = a ? "north" : "south";
    h2_zdc_e_ohcal[a] = book2("zdc_e_ohcal_" + arm, ("ZDC " + arm + " energy vs OHCal energy sum, " + zsel + ";#SigmaE_{OHCal} (GeV);E (" + zu + ")").c_str(), 200, 0, ohmax, 250, 0, zdcmax);
    h2_zdc_e_t0[a] = book2("zdc_e_t0_" + arm, ("ZDC " + arm + " energy vs MBD t0, " + zsel + ";t0 (ns);E (" + zu + ")").c_str(), 40, -10, 10, 250, 0, zdcmax);
    p2_zdc_t0_ohcal[a] = new TProfile2D(("p2_" + Name() + "_zdc_t0_ohcal_" + arm).c_str(), ("ZDC " + arm + " mean energy, " + zsel + ";t0 (ns);#SigmaE_{OHCal} (GeV)").c_str(), 40, -10, 10, 40, 0, ohmax);
    m_hm->registerHisto(p2_zdc_t0_ohcal[a]);
  }
  h2_zdc_e_ns_mb = book2("zdc_e_ns_mb", ("ZDC north vs south energy, MB;E_{south} (" + zu + ");E_{north} (" + zu + ")").c_str(), 200, 0, zdcmax, 200, 0, zdcmax);
  h2_zdc_e_ns = book2("zdc_e_ns", ("ZDC north vs south energy, ZDC trigger;E_{south} (" + zu + ");E_{north} (" + zu + ")").c_str(), 200, 0, zdcmax, 200, 0, zdcmax);
  const int nsq = std::lround((sepdqmax + 1.) / sepdqbinw);
  h2_sepd_q_ch = book2("sepd_q_ch", "sEPD tile charge, MB;channel;q (nMIP)", 744, -0.5, 743.5, nsq, -1., sepdqmax);
  h_sepd_qsum_south = book1("sepd_qsum_south", "sEPD south charge sum per event, MB;#Sigmaq (nMIP);events", 500, 0, sepdsummax);
  h_sepd_qsum_north = book1("sepd_qsum_north", "sEPD north charge sum per event, MB;#Sigmaq (nMIP);events", 500, 0, sepdsummax);
  h2_sepd_qsum_ns = book2("sepd_qsum_ns", "sEPD north vs south charge sum, MB;#Sigmaq_{south} (nMIP);#Sigmaq_{north} (nMIP)",
                          200, 0, sepdsummax, 200, 0, sepdsummax);
  h_sepd_nhit_south = book1("sepd_nhit_south", "sEPD south tiles with q > 0.5 nMIP, MB;N_{tile};events", 373, -0.5, 372.5);
  h_sepd_nhit_north = book1("sepd_nhit_north", "sEPD north tiles with q > 0.5 nMIP, MB;N_{tile};events", 373, -0.5, 372.5);

  for (int l = 0; l < kNLayer; l++)
  {
    const std::string &tag = layertag[l];
    const std::string &ln = layername[l];
    const int ne = std::lround((ehi[l] - elo[l]) / ebinw[l]) + 1;
    const double elo_c = elo[l] - 0.5 * ebinw[l];
    const double ehi_c = elo_c + ne * ebinw[l];
    h_tow_e[l] = book1(tag + "_e", (ln + " good tower energy, MB;E (GeV);towers").c_str(), ne, elo_c, ehi_c);
    h_tow_e_noZSXCalib[l] = book1(tag + "_e_noZSXCalib", (ln + " good tower energy, no ZSCrossCalib, MB;E (GeV);towers").c_str(), ne, elo_c, ehi_c);
    h_tow_esum[l] = book1(tag + "_esum", (ln + " good tower energy sum per event, MB;#SigmaE (GeV);events").c_str(), 1100, -0.1 * esummax, esummax);
    const double fracmax = (m_species == "pp") ? 0.05 : (m_species == "OO" ? 0.3 : 1.0);
    const double fbinw = (l == kCEMC) ? 1e-4 : 1. / 1536.;
    const int nfrac = std::lround(fracmax / fbinw) + 1;
    h_tow_occ[l] = book1(tag + "_occ", (ln + " fraction of good towers above threshold, MB;fraction;events").c_str(), nfrac, -0.5 * fbinw, (nfrac - 0.5) * fbinw);
    p_tow_occ_mean[l] = new TProfile(("p_" + Name() + "_" + tag + "_occ_mean").c_str(), (ln + " mean fraction of good towers above threshold, MB").c_str(), 1, 0, 1);
    m_hm->registerHisto(p_tow_occ_mean[l]);
    if (l != kCEMC)
    {
      p_tsc_e[l] = new TProfile(("p_" + Name() + "_" + tag + "_tsc_e").c_str(),
                                Form("%s good tower mean energy in %.1f-%.1f GeV, MB;tower (i#eta*64+i#phi);<E> (GeV)", ln.c_str(), m_tsclo[l], m_tschi[l]),
                                1536, -0.5, 1535.5);
      p_tsc_occ[l] = new TProfile(("p_" + Name() + "_" + tag + "_tsc_occ").c_str(),
                                  Form("%s good tower fraction of events with E > %.1f GeV, MB;tower (i#eta*64+i#phi);fraction", ln.c_str(), m_tsclo[l]),
                                  1536, -0.5, 1535.5);
      m_hm->registerHisto(p_tsc_e[l]);
      m_hm->registerHisto(p_tsc_occ[l]);
    }
    h_tow_nonzs_frac[l] = book1(tag + "_nonzs_frac", (ln + " fraction of non-ZS towers, MB;fraction;events").c_str(), nfrac, -0.5 * fbinw, (nfrac - 0.5) * fbinw);
    h_tow_time[l] = book1(tag + "_time", (ln + " tower time above threshold, MB;t (samples);towers").c_str(), 200, -10, 10);
    h_tow_time_mbdt0[l] = book1(tag + "_time_mbdt0", (ln + " tower time - MBD t0 above threshold, MB;t - t0 (ns);towers").c_str(), 200, -100, 100);
    h_noise_e_fine[l] = book1(tag + "_noise_e_fine", (ln + " good tower energy, MB |v_{z}|<100 cm;E (GeV);towers").c_str(), 2000, -0.25, 0.25);
    h2_noise_e_fine_ieta[l] = book2(tag + "_noise_e_fine_ieta", (ln + " good tower energy vs i#eta, MB |v_{z}|<100 cm;i#eta;E (GeV)").c_str(),
                                    netabins[l], -0.5, netabins[l] - 0.5, 2000, -0.25, 0.25);
  }

  for (int t = 0; t < kNTrig; t++)
  {
    h2_clus_e_eta[t] = book2("clus_e_eta_" + trigtag[t], ("EMCal cluster, " + trigtag[t] + ";E_{core} (GeV);#eta").c_str(), 300, 0, 30, 48, -1.2, 1.2);
    h2_jet_pt_eta[t] = book2("jet_pt_eta_" + trigtag[t], ("R=0.4 jet, " + trigtag[t] + ";p_{T} (GeV);#eta").c_str(), 80, 0, 80, 48, -1.2, 1.2);
  }
  if (m_species != "pp")
  {
    h_cent = book1("cent", "centrality, MB trigger and MB event selection;centrality (%);events", 100, 0, 100);
    h_cent_mbd_qsum = book1("cent_mbd_qsum", "MBD charge sum, MB trigger and MB event selection;MBD charge sum (N+S);events",
                            1000, 0, m_species == "OO" ? 500. : 3000.);
  }
  h_pi0_mass = book1("pi0_mass", "diphoton mass, MB;M_{#gamma#gamma} (GeV);pairs", 120, 0, 1.2);
}

void CDBCaloQA::captureCdbUrls()
{
  std::stringstream ss;
  std::streambuf *old = std::cout.rdbuf(ss.rdbuf());
  CDBInterface::instance()->Print();
  std::cout.rdbuf(old);
  std::string line;
  while (std::getline(ss, line))
  {
    const size_t d0 = line.find("domain: ");
    const size_t u0 = line.find(", url: ");
    const size_t t0 = line.find(", timestamp: ");
    if (d0 == std::string::npos || u0 == std::string::npos || t0 == std::string::npos)
    {
      continue;
    }
    m_cdburls[line.substr(d0 + 8, u0 - d0 - 8)] = line.substr(u0 + 7, t0 - u0 - 7);
  }
}

std::vector<int> CDBCaloQA::zdcBits(TriggerAnalyzer *trigana)
{
  std::vector<int> bits;
  for (int b = 0; b < 64; b++)
  {
    if (trigana->getTriggerName(b).find("ZDC") != std::string::npos && trigana->getTriggerPrescale(b) >= 0)
    {
      bits.push_back(b);
    }
  }
  return bits;
}

std::vector<int> CDBCaloQA::mbdBits(TriggerAnalyzer *trigana)
{
  std::vector<int> bits;
  for (int b = 0; b < 64; b++)
  {
    if (trigana->getTriggerName(b).find("MBD N&S") != std::string::npos)
    {
      bits.push_back(b);
    }
  }
  return bits;
}

bool CDBCaloQA::zdcEvent(TriggerAnalyzer *trigana, const std::vector<int> &zdcbits, const std::vector<int> &vetobits)
{
  bool zdc = false;
  for (int b : zdcbits)
  {
    zdc = zdc || trigana->didTriggerFire(b);
  }
  for (int b : vetobits)
  {
    zdc = zdc && !trigana->checkRawTrigger(b);
  }
  return zdc;
}

int CDBCaloQA::process_event(PHCompositeNode *topNode)
{
  if (!m_select && !m_cdbcaptured)
  {
    captureCdbUrls();
    m_cdbcaptured = true;
  }
  if (m_mc)
  {
    m_mbname = "MC, all events";
  }
  else
  {
    m_trigana->decodeTriggers(topNode);
  }
  if (m_mbname.empty())
  {
    if (m_mbbit_user < 0 && m_species == "AuAu" && m_mbbit == 14 && m_trigana->getTriggerPrescale(14) < 0)
    {
      m_mbbit = 10;
      std::cout << Name() << ": bit 14 is off in run " << m_runnumber << ", MB trigger bit 10" << std::endl;
    }
    m_mbname = m_trigana->getTriggerName(m_mbbit);
    const int ps = m_trigana->getTriggerPrescale(m_mbbit);
    if (m_mbbit_user < 0 && m_mbbit != 12 && ps > m_mbor_prescale && m_trigana->getTriggerPrescale(12) >= 0)
    {
      m_mbbit_or = 12;
      const std::string n12 = m_trigana->getTriggerName(m_mbbit_or);
      const size_t cut = n12.rfind(", ");
      m_mbname += " OR " + (cut == std::string::npos ? n12 : n12.substr(cut + 2));
    }
    std::cout << Name() << ": MB trigger bit " << m_mbbit << " prescale " << ps
              << (m_mbbit_or >= 0 ? ", bit 12 also MB" : "") << ": " << m_mbname << std::endl;
  }
  if (!m_mc && !m_zdcdecided)
  {
    m_zdcbits = zdcBits(m_trigana);
    for (int b : m_zdcbits)
    {
      m_zdcname += (m_zdcname.empty() ? "" : " OR ") + m_trigana->getTriggerName(b) + " (bit " + std::to_string(b) + ")";
    }
    if (m_species != "pp" && !m_zdcbits.empty())
    {
      m_mbdbits = mbdBits(m_trigana);
      m_zdcname += ", no MBD N&S";
    }
    std::cout << Name() << ": ZDC histograms from " << (m_zdcname.empty() ? "no events, no ZDC trigger on" : m_zdcname) << std::endl;
    m_zdcdecided = true;
  }

  std::array<bool, kNTrig> fired{};
  for (int t = 0; t < kNTrig; t++)
  {
    if (m_mc)
    {
      fired[t] = (t == kMB);
    }
    else if (t == kMB)
    {
      fired[t] = m_trigana->didTriggerFire(m_mbbit) || (m_mbbit_or >= 0 && m_trigana->didTriggerFire(m_mbbit_or));
    }
    else
    {
      fired[t] = !m_trigname[t].empty() && m_trigana->didTriggerFire(m_trigname[t]);
    }
  }
  fired[kMB] = fired[kMB] && (m_minmb <= 0 || m_nmb < m_minmb);
  fired[kPhoton] = fired[kPhoton] && m_ntrig[kPhoton] < m_mintrig;
  fired[kJet] = fired[kJet] && m_ntrig[kJet] < m_mintrig;
  const bool isMB = fired[kMB];
  const bool isZdc = !m_mc && m_nzdc < m_minzdc && zdcEvent(m_trigana, m_zdcbits, m_mbdbits);
  const bool keep = isMB || isZdc || fired[kPhoton] || fired[kJet];

  m_nevt++;
  if (keep)
  {
    m_nmb += isMB;
    m_nzdc += isZdc;
    m_ntrig[kPhoton] += fired[kPhoton];
    m_ntrig[kJet] += fired[kJet];
  }
  const bool done = (m_minmb > 0 || m_minzdc > 0 || m_mintrig > 0) &&
                    (m_minmb <= 0 || m_nmb >= m_minmb) && (m_zdcbits.empty() || m_nzdc >= m_minzdc) &&
                    (m_mc || m_trigname[kPhoton].empty() || m_ntrig[kPhoton] >= m_mintrig) &&
                    (m_mc || m_trigname[kJet].empty() || m_ntrig[kJet] >= m_mintrig);
  if (!keep)
  {
    return done ? Fun4AllReturnCodes::ABORTRUN : (m_select ? Fun4AllReturnCodes::ABORTEVENT : Fun4AllReturnCodes::EVENT_OK);
  }
  if (m_select)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  h_evt->Fill(0);
  for (int t = 0; t < kNTrig; t++)
  {
    if (fired[t])
    {
      h_evt->Fill(t + 1);
    }
  }

  m_vtxz = std::numeric_limits<float>::quiet_NaN();
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (vertexmap && !vertexmap->empty())
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
    if (vtx)
    {
      m_vtxz = vtx->get_z();
      h_vtxz->Fill(m_vtxz);
    }
  }

  m_noise_evt = isMB && std::isfinite(m_vtxz) && std::abs(m_vtxz) < m_noise_vzcut;
  if (m_noise_evt)
  {
    h_noise_nevt->Fill(0.5);
  }

  processMbd(topNode, isMB);
  if (isZdc)
  {
    processZdc(topNode);
  }
  Zdcinfo *zdcinfo = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
  if (isMB && zdcinfo)
  {
    h2_zdc_e_ns_mb->Fill(zdcinfo->get_zdc_energy(0), zdcinfo->get_zdc_energy(1));
    h_zdc_e_south_mb->Fill(zdcinfo->get_zdc_energy(0));
    h_zdc_e_north_mb->Fill(zdcinfo->get_zdc_energy(1));
    TowerInfoContainer *ohcal = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
    const bool zok = m_species == "pp" || (std::isfinite(m_vtxz) && std::abs(m_vtxz) < 10);
    if (h2_zdc_e_ohcal[0] && ohcal && zok)
    {
      float esum = 0;
      for (unsigned int ch = 0; ch < ohcal->size(); ch++)
      {
        TowerInfo *t = ohcal->get_tower_at_channel(ch);
        if (t->get_isGood())
        {
          esum += t->get_energy();
        }
      }
      for (int a = 0; a < 2; a++)
      {
        const float e = zdcinfo->get_zdc_energy(a);
        h2_zdc_e_ohcal[a]->Fill(esum, e);
        if (std::isfinite(m_mbdt0))
        {
          h2_zdc_e_t0[a]->Fill(m_mbdt0, e);
          p2_zdc_t0_ohcal[a]->Fill(m_mbdt0, esum, e);
        }
      }
    }
  }
  processSepd(topNode, isMB);
  for (int l = 0; isMB && l < kNLayer; l++)
  {
    processTowers(topNode, l);
  }
  processClusters(topNode, fired, isMB);
  processJets(topNode, fired);
  processCentrality(topNode, isMB);

  if (done)
  {
    std::cout << Name() << ": " << m_nevt << " events read, " << m_nmb << " MB, " << m_nzdc << " ZDC-triggered, "
              << m_ntrig[kPhoton] << " photon, " << m_ntrig[kJet] << " jet events, stop" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void CDBCaloQA::processMbd(PHCompositeNode *topNode, bool isMB)
{
  m_mbdt0 = std::numeric_limits<float>::quiet_NaN();
  MbdOut *mbdout = findNode::getClass<MbdOut>(topNode, "MbdOut");
  if (mbdout)
  {
    m_mbdt0 = mbdout->get_t0();
  }
  if (!isMB)
  {
    return;
  }
  if (mbdout)
  {
    h_mbd_zvtx->Fill(mbdout->get_zvtx());
    h_mbd_t0->Fill(m_mbdt0);
  }

  MbdPmtContainer *pmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if (!pmts)
  {
    return;
  }
  std::array<float, 2> qsum{0, 0};
  std::array<int, 2> nhit{0, 0};
  const int npmt = pmts->get_npmt();
  for (int i = 0; i < npmt; i++)
  {
    MbdPmtHit *pmt = pmts->get_pmt(i);
    const float q = pmt->get_q();
    const int arm = i / 64;
    if (std::isfinite(pmt->get_tt()))
    {
      h2_mbd_q_pmt->Fill(i, q);
      if (q > 0)
      {
        h2_mbd_t_pmt->Fill(i, pmt->get_time());
      }
    }
    qsum[arm] += q;
    if (q > 0.4)
    {
      nhit[arm]++;
    }
  }
  h_mbd_qavg_south->Fill(qsum[0] / 64.);
  h_mbd_qavg_north->Fill(qsum[1] / 64.);
  h_mbd_npmt_south->Fill(nhit[0]);
  h_mbd_npmt_north->Fill(nhit[1]);
}

void CDBCaloQA::processZdc(PHCompositeNode *topNode)
{
  Zdcinfo *zdcinfo = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
  if (zdcinfo)
  {
    const float es = zdcinfo->get_zdc_energy(0);
    const float en = zdcinfo->get_zdc_energy(1);
    if (es > 0)
    {
      h_zdc_e_south->Fill(es);
    }
    if (en > 0)
    {
      h_zdc_e_north->Fill(en);
    }
    h2_zdc_e_ns->Fill(es, en);
  }
  TowerInfoContainer *zdc = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_ZDC");
  if (!zdc || !m_zdccalib || zdc->size() < 16)
  {
    return;
  }
  for (int a = 0; a < 2; a++)
  {
    for (int m = 0; m < 3; m++)
    {
      const int ch = 8 * a + 2 * m;
      TowerInfo *tower = zdc->get_tower_at_channel(ch);
      const float t = tower->get_time();
      const float e = (t > 5. && t < 9.) ? tower->get_energy() * m_zdccalib->GetFloatValue(ch, "zdc_calib") : 0.;
      if (e > 0)
      {
        h_zdc_mod_e[3 * a + m]->Fill(e);
      }
    }
  }
}

void CDBCaloQA::processSepd(PHCompositeNode *topNode, bool isMB)
{
  if (!isMB)
  {
    return;
  }
  TowerInfoContainer *sepd = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_SEPD");
  if (!sepd)
  {
    sepd = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_EPD");
  }
  if (!sepd)
  {
    return;
  }
  std::array<float, 2> qsum{0, 0};
  std::array<int, 2> nhit{0, 0};
  const int size = sepd->size();
  for (int ch = 0; ch < size; ch++)
  {
    const float q = sepd->get_tower_at_channel(ch)->get_energy();
    const int arm = TowerInfoDefs::get_epd_arm(sepd->encode_key(ch));
    h2_sepd_q_ch->Fill(ch, q);
    qsum[arm] += q;
    if (q > m_sepd_hitthr)
    {
      nhit[arm]++;
    }
  }
  h_sepd_qsum_south->Fill(qsum[0]);
  h_sepd_qsum_north->Fill(qsum[1]);
  h2_sepd_qsum_ns->Fill(qsum[0], qsum[1]);
  h_sepd_nhit_south->Fill(nhit[0]);
  h_sepd_nhit_north->Fill(nhit[1]);
}

void CDBCaloQA::processTowers(PHCompositeNode *topNode, int layer)
{
  const std::string &ln = layername[layer];
  TowerInfoContainer *calib = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_" + ln);
  TowerInfoContainer *noxc = findNode::getClass<TowerInfoContainer>(topNode, m_noxc_prefix + ln);
  TowerInfoContainer *raw = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_" + ln);
  if (!calib)
  {
    return;
  }

  int ngood = 0;
  int nabove = 0;
  float esum = 0;
  const int size = calib->size();
  for (int ch = 0; ch < size; ch++)
  {
    TowerInfo *tower = calib->get_tower_at_channel(ch);
    if (!tower->get_isGood())
    {
      continue;
    }
    ngood++;
    const float e = tower->get_energy();
    esum += e;
    const int ieta = calib->getTowerEtaBin(calib->encode_key(ch));
    h_tow_e[layer]->Fill(e);
    if (p_tsc_e[layer])
    {
      const int itow = ieta * 64 + calib->getTowerPhiBin(calib->encode_key(ch));
      p_tsc_occ[layer]->Fill(itow, e > m_tsclo[layer] ? 1. : 0.);
      if (e > m_tsclo[layer] && e < m_tschi[layer])
      {
        p_tsc_e[layer]->Fill(itow, e);
      }
    }
    if (m_noise_evt)
    {
      h_noise_e_fine[layer]->Fill(e);
      h2_noise_e_fine_ieta[layer]->Fill(ieta, e);
    }
    if (noxc)
    {
      h_tow_e_noZSXCalib[layer]->Fill(noxc->get_tower_at_channel(ch)->get_energy());
    }
    if (e > m_towthr[layer])
    {
      nabove++;
      const float t = tower->get_time();
      h_tow_time[layer]->Fill(t);
      if (std::isfinite(m_mbdt0))
      {
        h_tow_time_mbdt0[layer]->Fill(t * m_ns_per_sample - m_mbdt0);
      }
    }
  }
  if (ngood > 0)
  {
    h_tow_occ[layer]->Fill(static_cast<float>(nabove) / ngood);
    p_tow_occ_mean[layer]->Fill(0.5, static_cast<float>(nabove) / ngood);
  }
  h_tow_esum[layer]->Fill(esum);

  if (raw && raw->size() > 0)
  {
    int nnonzs = 0;
    const int rsize = raw->size();
    for (int ch = 0; ch < rsize; ch++)
    {
      if (!raw->get_tower_at_channel(ch)->get_isZS())
      {
        nnonzs++;
      }
    }
    h_tow_nonzs_frac[layer]->Fill(static_cast<float>(nnonzs) / rsize);
  }
}

void CDBCaloQA::processClusters(PHCompositeNode *topNode, const std::array<bool, kNTrig> &fired, bool isMB)
{
  RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
  if (!clusters)
  {
    return;
  }
  const float vz = std::isfinite(m_vtxz) ? m_vtxz : 0;
  const CLHEP::Hep3Vector vertex(0, 0, vz);

  const float minClusPt1 = 1.5;
  const float minClusPt2 = 1.0;
  const float maxClusE = 100;
  const float maxChi2 = 4;
  const float maxAlpha = 0.6;
  const bool dopi0 = isMB && std::isfinite(m_vtxz) && std::abs(m_vtxz) < 30;

  std::vector<TLorentzVector> photons;
  RawClusterContainer::ConstRange range = clusters->getClusters();
  for (auto it = range.first; it != range.second; ++it)
  {
    RawCluster *clus = it->second;
    const CLHEP::Hep3Vector ev = RawClusterUtility::GetECoreVec(*clus, vertex);
    const float e = ev.mag();
    const float eta = ev.pseudoRapidity();
    for (int t = 0; t < kNTrig; t++)
    {
      if (fired[t])
      {
        h2_clus_e_eta[t]->Fill(e, eta);
      }
    }
    if (dopi0 && ev.perp() > minClusPt2 && e < maxClusE && clus->get_chi2() < maxChi2)
    {
      TLorentzVector p;
      p.SetPtEtaPhiE(ev.perp(), eta, ev.phi(), e);
      photons.push_back(p);
    }
  }
  for (size_t i = 0; i < photons.size(); i++)
  {
    for (size_t j = i + 1; j < photons.size(); j++)
    {
      if (std::max(photons[i].Pt(), photons[j].Pt()) < minClusPt1)
      {
        continue;
      }
      const float alpha = std::abs(photons[i].E() - photons[j].E()) / (photons[i].E() + photons[j].E());
      if (alpha > maxAlpha)
      {
        continue;
      }
      h_pi0_mass->Fill((photons[i] + photons[j]).M());
    }
  }
}

void CDBCaloQA::processJets(PHCompositeNode *topNode, const std::array<bool, kNTrig> &fired)
{
  JetContainer *jets = findNode::getClass<JetContainer>(topNode, m_jetnode);
  if (!jets)
  {
    return;
  }
  for (unsigned int i = 0; i < jets->size(); i++)
  {
    Jet *jet = jets->get_jet(i);
    for (int t = 0; t < kNTrig; t++)
    {
      if (fired[t])
      {
        h2_jet_pt_eta[t]->Fill(jet->get_pt(), jet->get_eta());
      }
    }
  }
}

void CDBCaloQA::processCentrality(PHCompositeNode *topNode, bool isMB)
{
  if (!h_cent || !isMB)
  {
    return;
  }
  MinimumBiasInfo *mbinfo = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");
  CentralityInfo *cent = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  MbdOut *mbdout = findNode::getClass<MbdOut>(topNode, "MbdOut");
  if (!mbinfo || !cent || !mbdout || !mbinfo->isAuAuMinimumBias())
  {
    return;
  }
  if (!(std::abs(mbdout->get_zvtx()) < m_cent_vzcut))
  {
    return;
  }
  h_cent->Fill(100. * cent->get_centile(CentralityInfo::PROP::mbd_NS));
  h_cent_mbd_qsum->Fill(mbdout->get_q(0) + mbdout->get_q(1));
}

int CDBCaloQA::End(PHCompositeNode * /*topNode*/)
{
  if (m_select)
  {
    std::cout << Name() << ": " << m_nevt << " events read, kept " << m_nmb << " MB, " << m_nzdc << " ZDC-triggered, "
              << m_ntrig[kPhoton] << " photon, " << m_ntrig[kJet] << " jet events" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  m_hm->registerHisto(new TNamed("trig_MB", m_mc ? m_mbname.c_str() : (m_mbbit_or >= 0 ? Form("%s (bits %d, %d)", m_mbname.c_str(), m_mbbit, m_mbbit_or) : Form("%s (bit %d)", m_mbname.c_str(), m_mbbit))));
  m_hm->registerHisto(new TNamed("trig_Photon", m_mc ? "" : m_trigname[kPhoton].c_str()));
  m_hm->registerHisto(new TNamed("trig_Jet", m_mc ? "" : m_trigname[kJet].c_str()));
  m_hm->registerHisto(new TNamed("trig_ZDC", m_zdcname.c_str()));

  captureCdbUrls();
  for (const auto &[domain, url] : m_cdburls)
  {
    m_hm->registerHisto(new TNamed(("cdb_" + domain).c_str(), url.c_str()));
  }
  m_hm->dumpHistos(m_outfile);
  return Fun4AllReturnCodes::EVENT_OK;
}
