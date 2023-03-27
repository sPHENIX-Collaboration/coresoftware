#include "BbcSimReco.h"

#include <bbc/BbcOutV1.h>
#include <bbc/BbcPmtContainerV1.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <ffaobjects/EventHeaderv1.h>

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

#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TTree.h>
#include <TDatabasePDG.h>
#include <TSystem.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cmath>
#include <iostream>

namespace BBCINFO
{
  const Double_t index_refract = 1.4585;
  const Double_t v_ckov = 1.0 / index_refract;  // velocity threshold for CKOV
  const Double_t C = 29.9792458;                // cm/ns

  // kludge where we have the hardcoded positions of the tubes
  // These are the x,y for the south BBC (in cm).
  // The north inverts the x coordinate (x -> -x)

  /* unused causes compiler warning
  const float TubeLoc[64][2] = {
      { -12.2976,	4.26 },
      { -12.2976,	1.42 },
      { -9.83805,	8.52 },
      { -9.83805,	5.68 },
      { -9.83805,	2.84 },
      { -7.37854,	9.94 },
      { -7.37854,	7.1 },
      { -7.37854,	4.26 },
      { -7.37854,	1.42 },
      { -4.91902,	11.36 },
      { -4.91902,	8.52 },
      { -4.91902,	5.68 },
      { -2.45951,	12.78 },
      { -2.45951,	9.94 },
      { -2.45951,	7.1 },
      { 0,	11.36 },
      { 0,	8.52 },
      { 2.45951,	12.78 },
      { 2.45951,	9.94 },
      { 2.45951,	7.1 },
      { 4.91902,	11.36 },
      { 4.91902,	8.52 },
      { 4.91902,	5.68 },
      { 7.37854,	9.94 },
      { 7.37854,	7.1 },
      { 7.37854,	4.26 },
      { 7.37854,	1.42 },
      { 9.83805,	8.52 },
      { 9.83805,	5.68 },
      { 9.83805,	2.84 },
      { 12.2976,	4.26 },
      { 12.2976,	1.42 },
      { 12.2976,	-4.26 },
      { 12.2976,	-1.42 },
      { 9.83805,	-8.52 },
      { 9.83805,	-5.68 },
      { 9.83805,	-2.84 },
      { 7.37854,	-9.94 },
      { 7.37854,	-7.1 },
      { 7.37854,	-4.26 },
      { 7.37854,	-1.42 },
      { 4.91902,	-11.36 },
      { 4.91902,	-8.52 },
      { 4.91902,	-5.68 },
      { 2.45951,	-12.78 },
      { 2.45951,	-9.94 },
      { 2.45951,	-7.1 },
      { 0,	-11.36 },
      { 0,	-8.52 },
      { -2.45951,	-12.78 },
      { -2.45951,	-9.94 },
      { -2.45951,	-7.1 },
      { -4.91902,	-11.36 },
      { -4.91902,	-8.52 },
      { -4.91902,	-5.68 },
      { -7.37854,	-9.94 },
      { -7.37854,	-7.1 },
      { -7.37854,	-4.26 },
      { -7.37854,	-1.42 },
      { -9.83805,	-8.52 },
      { -9.83805,	-5.68 },
      { -9.83805,	-2.84 },
      { -12.2976,	-4.26 },
      { -12.2976,	-1.42 }
  };
  */

}  // namespace BBCINFO

using namespace BBCINFO;

//____________________________________
BbcSimReco::BbcSimReco(const std::string &name)
  : SubsysReco(name)
  , _tres(0.05)
{
  std::fill(std::begin(f_pmtq), std::end(f_pmtq), NAN);
  std::fill(std::begin(f_pmtt0), std::end(f_pmtt0), NAN);
  std::fill(std::begin(f_pmtt1), std::end(f_pmtt1), NAN);
  std::fill(std::begin(f_bbcn), std::end(f_bbcn), 0);
  std::fill(std::begin(f_bbcq), std::end(f_bbcq), NAN);
  std::fill(std::begin(f_bbct), std::end(f_bbct), NAN);
  std::fill(std::begin(f_bbcte), std::end(f_bbcte), NAN);
  hevt_bbct[0] = nullptr;
  hevt_bbct[1] = nullptr;
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  m_Seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  gsl_rng_set(m_RandomGenerator, m_Seed);
}

BbcSimReco::~BbcSimReco()
{
  gsl_rng_free(m_RandomGenerator);
  return;
}

//___________________________________
int BbcSimReco::Init(PHCompositeNode *topNode)
{
  // std::cout << PHWHERE << std::endl;
  CreateNodes(topNode);

  _pdg = new TDatabasePDG();  // database of PDG info on particles

  TString name, title;
  for (int iarm = 0; iarm < 2; iarm++)
  {
    //
    name = "hevt_bbct";
    name += iarm;
    title = "bbc times, arm ";
    title += iarm;
    hevt_bbct[iarm] = new TH1F(name, title, 200, 7.5, 11.5);
    hevt_bbct[iarm]->SetLineColor(4);
  }

  gaussian = new TF1("gaussian", "gaus", 0, 20);
  gaussian->FixParameter(2, _tres);  // set sigma to timing resolution

  return 0;
}

//___________________________________
int BbcSimReco::InitRun(PHCompositeNode *topNode)
{
  GetNodes(topNode);

  return 0;
}

//__________________________________
// Call user instructions for every event
int BbcSimReco::process_event(PHCompositeNode * /*topNode*/)
{
  // GetNodes(topNode);

  f_evt = _evtheader->get_EvtSequence();
  // if(f_evt%100==0) std::cout << PHWHERE << "Events processed: " << f_evt << std::endl;

  //**** Initialize Variables

  // Arm Data
  f_bbcn[0] = 0;
  f_bbcn[1] = 0;
  f_bbcq[0] = 0.;
  f_bbcq[1] = 0.;
  f_bbct[0] = -9999.;
  f_bbct[1] = -9999.;
  f_bbcte[0] = -9999.;
  f_bbcte[1] = -9999.;
  f_bbcz = NAN;
  f_bbct0 = NAN;
  hevt_bbct[0]->Reset();
  hevt_bbct[1]->Reset();

  // PMT data
  float len[128] = {0.};
  float edep[128] = {0.};
  float first_time[128];  // First hit time for each tube
  std::fill_n(first_time, 128, 1e12);
  std::fill_n(f_pmtt0, 128, 1e12);
  std::fill_n(f_pmtt1, 128, 1e12);
  std::fill_n(f_pmtq, 128, 0.);

  // Get True Vertex
  // NB: Currently PrimaryVertexIndex is always 1, need to figure out how to handle collision pile-up
  PHG4VtxPoint *vtxp = _truth_container->GetPrimaryVtx(_truth_container->GetPrimaryVertexIndex());
  if (vtxp != nullptr)
  {
    f_vx = vtxp->get_x();
    f_vy = vtxp->get_y();
    f_vz = vtxp->get_z();
    f_vt = vtxp->get_t();

    if (f_evt < 20)
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
        Float_t dt = gsl_ran_gaussian(m_RandomGenerator, _tres); // get fluctuation in time
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
    if (beta > v_ckov && charge != 0.)
    {
      len[ch] += this_hit->get_path_length();

      _pids[pid] += 1;
    }

    nhits++;
  }

  // process the data from each channel
  for (float &iarm : f_bbct)
  {
    iarm = 0.;
  }

  std::vector<float> hit_times[2];  // times of the hits in each [arm]

  for (int ich = 0; ich < 128; ich++)
  {
    int arm = ich / 64;  // ch 0-63 = south, ch 64-127 = north

    // Fill charge and time info
    if (len[ich] > 0.)
    {
      if (f_evt < 20 && Verbosity() != 0)
      {
        std::cout << "ich " << ich << "\t" << len[ich] << "\t" << edep[ich] << std::endl;
      }

      // Get charge in BBC tube
      float npe = len[ich] * (120 / 3.0);                          // we get 120 p.e. per 3 cm
      float dnpe = gsl_ran_gaussian(m_RandomGenerator, std::sqrt(npe)); // get fluctuation in npe

      npe += dnpe;  // apply the fluctuations in npe

      f_bbcq[arm] += npe;
      f_pmtq[ich] = npe;

      // Now time
      if (first_time[ich] < 9999.)
      {
        // Fill evt histogram
        hevt_bbct[arm]->Fill(first_time[ich]);
        hit_times[arm].push_back(first_time[ich]);

        f_bbct[arm] += first_time[ich];
        // std::cout << "XXX " << ich << "\t" << f_bbct[arm] << "\t" << first_time[ich] << std::endl;

        f_pmtt0[ich] = first_time[ich];
        f_pmtt1[ich] = first_time[ich];
      }
      else  // should never happen
      {
        if (Verbosity() != 0)
        {
          std::cout << "ERROR, have hit but no time" << std::endl;
        }
      }

      _bbcpmts->AddBbcPmt(ich, f_pmtq[ich], f_pmtt0[ich], f_pmtt1[ich]);

      // threshold should be > 0.
      ++f_bbcn[arm];
    }
  }

  // Get best t
  if (f_bbcn[0] > 0 && f_bbcn[1] > 0)
  {
    for (int iarm = 0; iarm < 2; iarm++)
    {
      std::sort(hit_times[iarm].begin(), hit_times[iarm].end());
      float earliest = hit_times[iarm][0];

      gaussian->SetParameter(0, 5);
      gaussian->SetParameter(1, earliest);
      gaussian->SetRange(6, earliest + 5 * 0.05);
      // gaussian->SetParameter(1,hevt_bbct[iarm]->GetMean());
      // gaussian->SetRange(6,hevt_bbct[iarm]->GetMean()+0.125);

      hevt_bbct[iarm]->Fit(gaussian, "BLR");
      hevt_bbct[iarm]->Draw();

      if (f_bbcn[iarm] > 0)
      {
        // f_bbct[iarm] = f_bbct[iarm] / f_bbcn[iarm];
        f_bbct[iarm] = gaussian->GetParameter(1);
        f_bbcte[iarm] = earliest;

        _bbcout->AddBbcNS(iarm, f_bbcn[iarm], f_bbcq[iarm], f_bbct[iarm]);
      }
    }

    // Now calculate zvtx, t0 from best times
    f_bbcz = (f_bbct[0] - f_bbct[1]) * C / 2.0;
    f_bbct0 = (f_bbct[0] + f_bbct[1]) / 2.0;

    _bbcout->set_Vertex(f_bbcz, 0.6);
    _bbcout->set_TimeZero(f_bbct0, 0.05);
  }

  return 0;
}

void BbcSimReco::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
  }

  PHCompositeNode *bbcNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "BBC"));
  if (!bbcNode)
  {
    bbcNode = new PHCompositeNode("BBC");
    dstNode->addNode(bbcNode);
  }

  //-* contains final physics products (nmips, t0, etc)
  BbcOut *bbcout = findNode::getClass<BbcOut>(bbcNode, "BbcOut");
  if (!bbcout)
  {
    bbcout = new BbcOutV1();
    PHIODataNode<PHObject> *BbcOutNode = new PHIODataNode<PHObject>(bbcout, "BbcOut", "PHObject");
    bbcNode->addNode(BbcOutNode);
  }

  //-* contains info for each pmt (nmips, time, etc)
  BbcPmtContainer *bbcpmts = findNode::getClass<BbcPmtContainer>(bbcNode, "BbcPmtContainer");
  if (!bbcpmts)
  {
    bbcpmts = new BbcPmtContainerV1();
    PHIODataNode<PHObject> *BbcPmtNode = new PHIODataNode<PHObject>(bbcpmts, "BbcPmtContainer", "PHObject");
    bbcNode->addNode(BbcPmtNode);
  }
}

//___________________________________
void BbcSimReco::GetNodes(PHCompositeNode *topNode)
{
  // Get the DST objects

  // Truth container
  _truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_truth_container && f_evt < 10)
  {
    std::cout << PHWHERE << " PHG4TruthInfoContainer node not found on node tree" << std::endl;
  }

  // BBC hit container
  _bbchits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_BBC");
  if (!_bbchits && f_evt < 10)
  {
    std::cout << PHWHERE << " G4HIT_BBC node not found on node tree" << std::endl;
  }

  // Event Header info
  _evtheader = findNode::getClass<EventHeaderv1>(topNode, "EventHeader");
  if (!_evtheader && f_evt < 10)
  {
    std::cout << PHWHERE << " G4HIT_BBC node not found on node tree" << std::endl;
  }

  /** DST Objects **/

  // BbcOut data
  _bbcout = findNode::getClass<BbcOut>(topNode, "BbcOut");
  if (!_bbcout && f_evt < 10)
  {
    std::cout << PHWHERE << " BbcOut node not found on node tree" << std::endl;
  }

  // BbcPmtContainer
  _bbcpmts = findNode::getClass<BbcPmtContainer>(topNode, "BbcPmtContainer");
  if (!_bbcpmts && f_evt < 10)
  {
    std::cout << PHWHERE << " BbcPmtContainer node not found on node tree" << std::endl;
  }
}
