#include "MbdDigitization.h"

#include <mbd/MbdDefs.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtHit.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cmath>
#include <iostream>

//____________________________________
MbdDigitization::MbdDigitization(const std::string &name)
  : SubsysReco(name)
  , _tres(0.05)
{
  std::fill(std::begin(f_pmtq), std::end(f_pmtq), 0.);
  std::fill(std::begin(f_pmtt0), std::end(f_pmtt0), std::numeric_limits<float>::quiet_NaN());
  std::fill(std::begin(f_pmtt1), std::end(f_pmtt1), std::numeric_limits<float>::quiet_NaN());
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  m_Seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  gsl_rng_set(m_RandomGenerator, m_Seed);
}

MbdDigitization::~MbdDigitization()
{
  gsl_rng_free(m_RandomGenerator);
  return;
}

//___________________________________
int MbdDigitization::Init(PHCompositeNode *topNode)
{
  // std::cout << PHWHERE << std::endl;
  CreateNodes(topNode);

  _pdg = new TDatabasePDG();  // database of PDG info on particles

  gaussian = new TF1("gaussian", "gaus", 0, 20);
  gaussian->FixParameter(2, _tres);  // set sigma to timing resolution

  return 0;
}

//___________________________________
int MbdDigitization::InitRun(PHCompositeNode *topNode)
{
  GetNodes(topNode);

  return 0;
}

//__________________________________
// Call user instructions for every event
int MbdDigitization::process_event(PHCompositeNode * /*topNode*/)
{
  //**** Initialize Variables

  // PMT data
  float len[MbdDefs::MBD_N_PMT] = {0.};
  float edep[MbdDefs::MBD_N_PMT] = {0.};
  float first_time[MbdDefs::MBD_N_PMT];  // First hit time for each tube
  std::fill_n(first_time, MbdDefs::MBD_N_PMT, 1e12);
  std::fill_n(f_pmtt0, MbdDefs::MBD_N_PMT, 1e12);
  std::fill_n(f_pmtt1, MbdDefs::MBD_N_PMT, 1e12);
  std::fill_n(f_pmtq, MbdDefs::MBD_N_PMT, 0.);

  // Get True Vertex
  // NB: Currently PrimaryVertexIndex is always 1, need to figure out how to handle collision pile-up
  PHG4VtxPoint *vtxp = _truth_container->GetPrimaryVtx(_truth_container->GetPrimaryVertexIndex());
  if (vtxp != nullptr)
  {
    f_vx = vtxp->get_x();
    f_vy = vtxp->get_y();
    f_vz = vtxp->get_z();
    f_vt = vtxp->get_t();

    if (Verbosity())
    {
      std::cout << "VTXP "
                << "\t" << f_vx << "\t" << f_vy << "\t" << f_vz << "\t" << f_vt << std::endl;
    }
  }

  // Go through BBC G4 hits

  TLorentzVector v4;
  unsigned int nhits = 0;

  PHG4HitContainer::ConstRange bbc_hit_range = _bbchits->getHits();
  for (auto hit_iter = bbc_hit_range.first; hit_iter != bbc_hit_range.second; hit_iter++)
  {
    PHG4Hit *this_hit = hit_iter->second;

    unsigned int ch = this_hit->get_layer();  // pmt channel
    // int arm = ch/64;;                         // south=0, north=1

    int trkid = this_hit->get_trkid();
    // if ( trkid>0 && f_evt<20 ) std::cout << "TRKID " << trkid << std::endl;

    PHG4Particle *part = _truth_container->GetParticle(trkid);
    v4.SetPxPyPzE(part->get_px(), part->get_py(), part->get_pz(), part->get_e());

    int pid = part->get_pid();
    TParticlePDG *partinfo = _pdg->GetParticle(pid);
    Double_t charge = -9999.;
    if (partinfo)
    {
      charge = partinfo->Charge() / 3;  // PDG gives charge in 1/3 e
    }
    else if (part->isIon())
    {
      charge = part->get_IonCharge();
    }

    // get the first time
    if (this_hit->get_t(1) < first_time[ch])
    {
      if (fabs(this_hit->get_t(1)) < 106.5)
      {
        first_time[ch] = this_hit->get_t(1) - vtxp->get_t();
        Float_t dt = gsl_ran_gaussian(m_RandomGenerator, _tres);  // get fluctuation in time
        first_time[ch] += dt;
      }
      else
      {
        if (Verbosity())
        {
          std::cout << "BAD " << ch << "\t" << this_hit->get_t(1) << std::endl;
        }
      }
    }

    edep[ch] += this_hit->get_edep();

    // get summed path length for particles that can create CKOV light
    // n.p.e. is determined from path length
    Double_t beta = v4.Beta();
    if (beta > MbdDefs::v_ckov && charge != 0.)
    {
      len[ch] += this_hit->get_path_length();

      _pids[pid] += 1;
    }

    nhits++;
  }

  for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ipmt++)
  {
    // Fill charge and time info
    if (len[ipmt] > 0.)
    {
      if (Verbosity() > 0)
      {
        std::cout << "ipmt " << ipmt << "\t" << len[ipmt] << "\t" << edep[ipmt] << std::endl;
      }

      // Get charge in BBC tube
      float npe = len[ipmt] * (120 / 3.0);                                // we get 120 p.e. per 3 cm
      float dnpe = gsl_ran_gaussian(m_RandomGenerator, std::sqrt(npe));  // get fluctuation in npe

      npe += dnpe;  // apply the fluctuations in npe
      f_pmtq[ipmt] = npe;

      // Now time
      if (first_time[ipmt] < 9999.)
      {
        f_pmtt0[ipmt] = first_time[ipmt] - 8.346;
        f_pmtt1[ipmt] = first_time[ipmt] - 8.346;
      }
      else  // should never happen
      {
        if (Verbosity() != 0)
        {
          std::cout << "ERROR, have hit but no time" << std::endl;
        }
      }

      _bbcpmts->get_pmt(ipmt)->set_pmt(ipmt, f_pmtq[ipmt], f_pmtt0[ipmt], f_pmtt1[ipmt]);

      if (Verbosity() > 0)
      {
        std::cout << "Adding " << ipmt << ", " << f_pmtq[ipmt] << ", " << f_pmtt0[ipmt] << " , " << f_pmtt1[ipmt] << std::endl;
      }
    }
    else
    {
      // empty channel
      _bbcpmts->get_pmt(ipmt)->set_pmt(ipmt, f_pmtq[ipmt], f_pmtt0[ipmt], f_pmtt1[ipmt]);
    }
  }

  return 0;
}

void MbdDigitization::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  PHCompositeNode *bbcNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "MBD"));
  if (!bbcNode)
  {
    bbcNode = new PHCompositeNode("MBD");
    dstNode->addNode(bbcNode);
  }

  //-* contains info for each pmt (nmips, time, etc)
  _bbcpmts = findNode::getClass<MbdPmtContainer>(bbcNode, "MbdPmtContainer");
  if (!_bbcpmts)
  {
    _bbcpmts = new MbdPmtContainerV1();
    PHIODataNode<PHObject> *BbcPmtNode = new PHIODataNode<PHObject>(_bbcpmts, "MbdPmtContainer", "PHObject");
    bbcNode->addNode(BbcPmtNode);
  }
}

//___________________________________
void MbdDigitization::GetNodes(PHCompositeNode *topNode)
{
  // Get the DST objects

  // Truth container
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth_container)
  {
    std::cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << std::endl;
    gSystem->Exit(1);
  }

  // BBC hit container
  _bbchits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_BBC");
  if (!_bbchits)
  {
    std::cout << PHWHERE << " G4HIT_BBC node not found on node tree" << std::endl;
    gSystem->Exit(1);
  }

  /** DST Objects **/

  // MbdPmtContainer
  _bbcpmts = findNode::getClass<MbdPmtContainerV1>(topNode, "MbdPmtContainer");
  if (!_bbcpmts)
  {
    std::cout << PHWHERE << " MbdPmtContainer node not found on node tree" << std::endl;
    gSystem->Exit(1);
  }
}
