// this is the new containers version
// it uses the same MapToPadPlane as the old containers version

#include "PHG4TpcElectronDrift.h"
#include "PHG4TpcDistortion.h"
#include "PHG4TpcPadPlane.h"  // for PHG4TpcPadPlane

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particlev3.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>  // for TrkrHit
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitTruthAssoc.h>  // for TrkrHitTruthA...
#include <trackbase/TrkrHitTruthAssocv1.h>
#include <trackbase/TrkrHitv2.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrTruthTrackv1.h>

#include <trackbase/TrkrTruthTrackContainerv1.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterv4.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterContainerv4.h>

#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>

#include <phparameter/PHParameterInterface.h>  // for PHParameterIn...
#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <pdbcalbase/PdbParameterMapContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>  // for PHDataNode
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TSystem.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_alloc

#include <array>
#include <cassert>
#include <cmath>    // for sqrt, abs, NAN
#include <cstdlib>  // for exit
#include <iostream>
#include <map>      // for _Rb_tree_cons...
#include <utility>  // for pair

#include "TpcClusterBuilder.h"

namespace
{
  template <class T>
  inline constexpr T square(const T &x)
  {
    return x * x;
  }
}  // namespace

PHG4TpcElectronDrift::PHG4TpcElectronDrift(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , temp_hitsetcontainer(new TrkrHitSetContainerv1)
  , single_hitsetcontainer(new TrkrHitSetContainerv1)
{
  InitializeParameters();
  RandomGenerator.reset(gsl_rng_alloc(gsl_rng_mt19937));
  set_seed(PHRandomSeed());
}

//_____________________________________________________________
int PHG4TpcElectronDrift::Init(PHCompositeNode *topNode)
{
  padplane->Init(topNode);
  event_num = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
int PHG4TpcElectronDrift::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }
  auto runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  auto parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  const std::string paramnodename = "G4CELLPARAM_" + detector;
  const std::string geonodename = "G4CELLPAR_" + detector;
  const std::string tpcgeonodename = "G4GEO_" + detector;
  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  if (!g4hit)
  {
    std::cout << Name() << " Could not locate G4HIT node " << hitnodename << std::endl;
    topNode->print();
    gSystem->Exit(1);
    exit(1);
  }
  // new containers
  hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsetcontainer)
  {
    PHNodeIterator dstiter(dstNode);
    auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    hitsetcontainer = new TrkrHitSetContainerv1;
    auto newNode = new PHIODataNode<PHObject>(hitsetcontainer, "TRKR_HITSET", "PHObject");
    DetNode->addNode(newNode);
  }

  hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if (!hittruthassoc)
  {
    PHNodeIterator dstiter(dstNode);
    auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    hittruthassoc = new TrkrHitTruthAssocv1;
    auto newNode = new PHIODataNode<PHObject>(hittruthassoc, "TRKR_HITTRUTHASSOC", "PHObject");
    DetNode->addNode(newNode);
  }

  truthtracks = findNode::getClass<TrkrTruthTrackContainer>(topNode, "TRKR_TRUTHTRACKCONTAINER");
  if (!truthtracks)
  {
    PHNodeIterator dstiter(dstNode);
    auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    truthtracks = new TrkrTruthTrackContainerv1;
    auto newNode = new PHIODataNode<PHObject>(truthtracks, "TRKR_TRUTHTRACKCONTAINER", "PHObject");
    DetNode->addNode(newNode);
  }
  truthclustercontainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_TRUTHCLUSTERCONTAINER");
  if (!truthclustercontainer)
  {
    PHNodeIterator dstiter(dstNode);
    auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    truthclustercontainer = new TrkrClusterContainerv4;
    auto newNode = new PHIODataNode<PHObject>(truthclustercontainer, "TRKR_TRUTHCLUSTERCONTAINER", "PHObject");
    DetNode->addNode(newNode);
  }

  seggeonodename = "CYLINDERCELLGEOM_SVTX";  // + detector;
  seggeo = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, seggeonodename);
  if (!seggeo)
  {
    seggeo = new PHG4TpcCylinderGeomContainer();
    auto newNode = new PHIODataNode<PHObject>(seggeo, seggeonodename, "PHObject");
    runNode->addNode(newNode);
  }

  
  UpdateParametersWithMacro();
  PHNodeIterator runIter(runNode);
  auto RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(detector);
    runNode->addNode(RunDetNode);
  }
  SaveToNodeTree(RunDetNode, paramnodename);

  // save this to the parNode for use
  PHNodeIterator parIter(parNode);
  auto ParDetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", detector));
  if (!ParDetNode)
  {
    ParDetNode = new PHCompositeNode(detector);
    parNode->addNode(ParDetNode);
  }
  PutOnParNode(ParDetNode, geonodename);

  // find Tpc Geo
  PHNodeIterator tpcpariter(ParDetNode);
  auto tpcparams = findNode::getClass<PHParametersContainer>(ParDetNode, tpcgeonodename);
  if (!tpcparams)
  {
    const std::string runparamname = "G4GEOPARAM_" + detector;
    auto tpcpdbparams = findNode::getClass<PdbParameterMapContainer>(RunDetNode, runparamname);
    if (tpcpdbparams)
    {
      tpcparams = new PHParametersContainer(detector);
      if (Verbosity()) tpcpdbparams->print();
      tpcparams->CreateAndFillFrom(tpcpdbparams, detector);
      ParDetNode->addNode(new PHDataNode<PHParametersContainer>(tpcparams, tpcgeonodename));
    }
    else
    {
      std::cout << "PHG4TpcElectronDrift::InitRun - failed to find " << runparamname << " in order to initialize " << tpcgeonodename << ". Aborting run ..." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  assert(tpcparams);

  if (Verbosity()) tpcparams->Print();
  const PHParameters *tpcparam = tpcparams->GetParameters(0);
  assert(tpcparam);
  tpc_length = tpcparam->get_double_param("tpc_length");

  diffusion_long = get_double_param("diffusion_long");
  added_smear_sigma_long = get_double_param("added_smear_long");
  diffusion_trans = get_double_param("diffusion_trans");
  added_smear_sigma_trans = get_double_param("added_smear_trans");
  drift_velocity = get_double_param("drift_velocity");
  // min_time to max_time is the time window for accepting drifted electrons after the trigger
  min_time = 0.0;
  max_time = get_double_param("max_time") + get_double_param("extended_readout_time");
  electrons_per_gev = get_double_param("electrons_per_gev");
  min_active_radius = get_double_param("min_active_radius");
  max_active_radius = get_double_param("max_active_radius");

  if (Verbosity() > 0) {
    std::cout << PHWHERE << " drift velocity " << drift_velocity << " extended_readout_time " << get_double_param("extended_readout_time") << " max time cutoff " << max_time << std::endl;
  }

  auto se = Fun4AllServer::instance();
  dlong = new TH1F("difflong", "longitudinal diffusion", 100, diffusion_long - diffusion_long / 2., diffusion_long + diffusion_long / 2.);
  se->registerHisto(dlong);
  dtrans = new TH1F("difftrans", "transversal diffusion", 100, diffusion_trans - diffusion_trans / 2., diffusion_trans + diffusion_trans / 2.);
  se->registerHisto(dtrans);

  do_ElectronDriftQAHistos = false;  // Whether or not to produce an ElectronDriftQA.root file with useful info
  if (do_ElectronDriftQAHistos)
  {
    hitmapstart = new TH2F("hitmapstart", "g4hit starting X-Y locations", 1560, -78, 78, 1560, -78, 78);
    hitmapend = new TH2F("hitmapend", "g4hit final X-Y locations", 1560, -78, 78, 1560, -78, 78);
    z_startmap = new TH2F("z_startmap", "g4hit starting Z vs. R locations", 2000, -100, 100, 780, 0, 78);
    deltaphi = new TH2F("deltaphi", "Total delta phi; phi (rad);#Delta phi (rad)", 600, -M_PI, M_PI, 1000, -.2, .2);
    deltaRphinodiff = new TH2F("deltaRphinodiff", "Total delta R*phi, no diffusion; r (cm);#Delta R*phi (cm)", 600, 20, 80, 1000, -3, 5);
    deltaphivsRnodiff = new TH2F("deltaphivsRnodiff", "Total delta phi vs. R; phi (rad);#Delta phi (rad)", 600, 20, 80, 1000, -.2, .2);
    deltaz = new TH2F("deltaz", "Total delta z; z (cm);#Delta z (cm)", 1000, 0, 100, 1000, -.5, 5);
    deltaphinodiff = new TH2F("deltaphinodiff", "Total delta phi (no diffusion, only SC distortion); phi (rad);#Delta phi (rad)", 600, -M_PI, M_PI, 1000, -.2, .2);
    deltaphinodist = new TH2F("deltaphinodist", "Total delta phi (no SC distortion, only diffusion); phi (rad);#Delta phi (rad)", 600, -M_PI, M_PI, 1000, -.2, .2);
    deltar = new TH2F("deltar", "Total Delta r; r (cm);#Delta r (cm)", 580, 20, 78, 1000, -3, 5);
    deltarnodiff = new TH2F("deltarnodiff", "Delta r (no diffusion, only SC distortion); r (cm);#Delta r (cm)", 580, 20, 78, 1000, -2, 5);
    deltarnodist = new TH2F("deltarnodist", "Delta r (no SC distortion, only diffusion); r (cm);#Delta r (cm)", 580, 20, 78, 1000, -2, 5);
  }

  if (Verbosity())
  {
    // eval tree only when verbosity is on
    m_outf.reset(new TFile("nt_out.root", "recreate"));
    nt = new TNtuple("nt", "electron drift stuff", "hit:ts:tb:tsig:rad:zstart:zfinal");
    nthit = new TNtuple("nthit", "TrkrHit collecting", "layer:phipad:zbin:neffelectrons");
    ntfinalhit = new TNtuple("ntfinalhit", "TrkrHit collecting", "layer:phipad:zbin:neffelectrons");
    ntpad = new TNtuple("ntpad", "electron by electron pad centroid", "layer:phigem:phiclus:zgem:zclus");
    se->registerHisto(nt);
    se->registerHisto(nthit);
    se->registerHisto(ntpad);
  }

  padplane->InitRun(topNode);
  padplane->CreateReadoutGeometry(topNode, seggeo);

  // print all layers radii
  if (Verbosity())
  {
    const auto range = seggeo->get_begin_end(); 
    std::cout << "PHG4TpcElectronDrift::InitRun - layers: " << std::distance(range.first, range.second) << std::endl;
    int counter = 0;
    for (auto layeriter = range.first; layeriter != range.second; ++layeriter)
    {
      const auto radius = layeriter->second->get_radius();
      std::cout << Form( "%.3f ", radius );
      if( ++counter == 8 )
      {
        counter = 0;
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TpcElectronDrift::process_event(PHCompositeNode *topNode)
{
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if(!m_tGeometry)
  {
    std::cout << PHWHERE << "ActsGeometry not found on node tree. Exiting" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (truth_clusterer == nullptr)  {
    if (Verbosity()) std::cout << " truth clusterer was a null pointer " << std::endl;
    truth_clusterer = new TpcClusterBuilder(truthclustercontainer, m_tGeometry, seggeo);
  } else {
    if (Verbosity()) std::cout << " truth clusterer was NOT a null pointer " << std::endl;
  }


  static constexpr unsigned int print_layer = 18;
  truth_clusterer->is_embedded_track = false;
  std::map<TrkrDefs::hitsetkey,unsigned int> hitset_cnt; // needed for indexing the TrkrClusters into the TrkrClusterContainer
  /* std::map<TrkrDefs::hitsetkey, std::vector<TrkrHitv*>> truthtrack_hits; */


  // tells m_distortionMap which event to look at
  if (m_distortionMap)
  {
    m_distortionMap->load_event(event_num);
  }

  // g4hits
  auto g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename);
  if (!g4hit)
  {
    std::cout << "Could not locate g4 hit node " << hitnodename << std::endl;
    gSystem->Exit(1);
  }
  PHG4TruthInfoContainer *truthinfo =
    findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");


  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE
		<< "ActsGeometry not found on node tree. Exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  unsigned int count_g4hits = 0;
  //  int count_electrons = 0;

  //  double ecollectedhits = 0.0;
  int ncollectedhits = 0;
  double ihit = 0;
  unsigned int dump_interval = 5000;  // dump temp_hitsetcontainer to the node tree after this many g4hits
  unsigned int dump_counter = 0;

  int trkid = -1;

  for (auto hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
  {
    count_g4hits++;
    dump_counter++;

    const double t0 = fmax(hiter->second->get_t(0), hiter->second->get_t(1));
    if (t0 > max_time)
    {
      continue;
    }

    int trkid_new = hiter->second->get_trkid();
    if (trkid != trkid_new)
    {  // starting a new track
      truth_clusterer->cluster_and_reset(/*argument is if to reset hitsetkey as well*/ false);
      trkid = trkid_new;
      truth_clusterer->is_embedded_track = (truthinfo->isEmbeded(hiter->second->get_trkid()));
      if (Verbosity() > 1000){
        std::cout << " New track " << trkid << " is embed? : " 
          << truth_clusterer->is_embedded_track << std::endl;
      }
      if (truth_clusterer->is_embedded_track) 
      { // build new TrkrTruthTrack
        auto particle = /*(PHG4Particlev3*)*/ truthinfo->GetParticle(trkid);
        int vtxid = particle->get_vtx_id();
        PHG4VtxPoint* vtx = truthinfo->GetVtx(vtxid);
        current_track = new TrkrTruthTrackv1(trkid, particle, vtx) ;
        truthtracks->addTruthTrack(current_track);
        truth_clusterer->set_current_track(current_track);
      }
    }


    // for very high occupancy events, accessing the TrkrHitsets on the node tree 
    // for every drifted electron seems to be very slow
    // Instead, use a temporary map to accumulate the charge from all 
    // drifted electrons, then copy to the node tree later
    
    double eion = hiter->second->get_eion();
    unsigned int n_electrons = gsl_ran_poisson(RandomGenerator.get(), eion * electrons_per_gev);
    //    count_electrons += n_electrons;

    if (Verbosity() > 100)
      std::cout << "  new hit with t0, " << t0 << " g4hitid " << hiter->first
                << " eion " << eion << " n_electrons " << n_electrons
                << " entry z " << hiter->second->get_z(0) << " exit z " 
                << hiter->second->get_z(1) << " avg z" 
                << (hiter->second->get_z(0) + hiter->second->get_z(1)) / 2.0
                << std::endl;

    if (n_electrons == 0) { continue; }

    if (Verbosity() > 100)
    {
      std::cout << std::endl
                << "electron drift: g4hit " << hiter->first << " created electrons: " 
                << n_electrons << " from " << eion * 1000000 << " keV" << std::endl;
      std::cout << " entry x,y,z = " << hiter->second->get_x(0) << "  " 
                << hiter->second->get_y(0) << "  " << hiter->second->get_z(0)
                << " radius " << sqrt(pow(hiter->second->get_x(0), 2) + 
                     pow(hiter->second->get_y(0), 2)) << std::endl;
      std::cout << " exit x,y,z = " << hiter->second->get_x(1) << "  " 
                << hiter->second->get_y(1) << "  " << hiter->second->get_z(1)
                << " radius " << sqrt(pow(hiter->second->get_x(1), 2) + 
                    pow(hiter->second->get_y(1), 2)) << std::endl;
    }

    for (unsigned int i = 0; i < n_electrons; i++)
    {
      // We choose the electron starting position at random from a flat
      // distribution along the path length the parameter t is the fraction of
      // the distance along the path betwen entry and exit points, it has
      // values between 0 and 1
      const double f = gsl_ran_flat(RandomGenerator.get(), 0.0, 1.0);

      const double x_start = hiter->second->get_x(0) + f * (hiter->second->get_x(1) - hiter->second->get_x(0));
      const double y_start = hiter->second->get_y(0) + f * (hiter->second->get_y(1) - hiter->second->get_y(0));
      const double z_start = hiter->second->get_z(0) + f * (hiter->second->get_z(1) - hiter->second->get_z(0));
      const double t_start = hiter->second->get_t(0) + f * (hiter->second->get_t(1) - hiter->second->get_t(0));

      unsigned int side = 0;
      if (z_start > 0) side = 1;

      const double r_sigma = diffusion_trans * sqrt(tpc_length / 2. - std::abs(z_start));
      const double rantrans =
          gsl_ran_gaussian(RandomGenerator.get(), r_sigma) +
          gsl_ran_gaussian(RandomGenerator.get(), added_smear_sigma_trans);

      const double t_path = (tpc_length / 2. - std::abs(z_start)) / drift_velocity;
      const double t_sigma = diffusion_long * sqrt(tpc_length / 2. - std::abs(z_start)) / drift_velocity;
      const double rantime =
          gsl_ran_gaussian(RandomGenerator.get(), t_sigma) +
          gsl_ran_gaussian(RandomGenerator.get(), added_smear_sigma_long) / drift_velocity;
      double t_final = t_start + t_path + rantime;

      if (t_final < min_time || t_final > max_time) continue;

      double z_final;
      if (z_start < 0)
        z_final = -tpc_length / 2. + t_final * drift_velocity;
      else
        z_final = tpc_length / 2. - t_final * drift_velocity;

      const double radstart = std::sqrt(square(x_start) + square(y_start));
      const double phistart = std::atan2(y_start, x_start);
      const double ranphi = gsl_ran_flat(RandomGenerator.get(), -M_PI, M_PI);

      double x_final = x_start + rantrans * std::cos(ranphi);  // Initialize these to be only diffused first, will be overwritten if doing SC distortion
      double y_final = y_start + rantrans * std::sin(ranphi);

      double rad_final = sqrt(square(x_final) + square(y_final));
      double phi_final = atan2(y_final, x_final);

      if (do_ElectronDriftQAHistos)
      {
        z_startmap->Fill(z_start, radstart);                   // map of starting location in Z vs. R
        deltaphinodist->Fill(phistart, rantrans / rad_final);  // delta phi no distortion, just diffusion+smear
        deltarnodist->Fill(radstart, rantrans);                // delta r no distortion, just diffusion+smear
      }

      if (m_distortionMap)
      {
        const double r_distortion = m_distortionMap->get_r_distortion(radstart, phistart, z_start);
        const double phi_distortion = m_distortionMap->get_rphi_distortion(radstart, phistart, z_start) / radstart;
        const double z_distortion = m_distortionMap->get_z_distortion(radstart, phistart, z_start);

        rad_final += r_distortion;
        phi_final += phi_distortion;
        z_final += z_distortion;
        if (z_start < 0)
          t_final = (z_final + tpc_length / 2.0) / drift_velocity;
        else
          t_final = (tpc_length / 2.0 - z_final) / drift_velocity;

        x_final = rad_final * std::cos(phi_final);
        y_final = rad_final * std::sin(phi_final);

        if (do_ElectronDriftQAHistos)
        {
          const double phi_final_nodiff = phistart + phi_distortion;
          const double rad_final_nodiff = radstart + r_distortion;
          deltarnodiff->Fill(radstart, rad_final_nodiff - radstart);    //delta r no diffusion, just distortion
          deltaphinodiff->Fill(phistart, phi_final_nodiff - phistart);  //delta phi no diffusion, just distortion
          deltaphivsRnodiff->Fill(radstart, phi_final_nodiff - phistart);
          deltaRphinodiff->Fill(radstart, rad_final_nodiff * phi_final_nodiff - radstart * phistart);

          // Fill Diagnostic plots, written into ElectronDriftQA.root
          hitmapstart->Fill(x_start, y_start);             // G4Hit starting positions
          hitmapend->Fill(x_final, y_final);               //INcludes diffusion and distortion
          deltar->Fill(radstart, rad_final - radstart);    //total delta r
          deltaphi->Fill(phistart, phi_final - phistart);  // total delta phi
          deltaz->Fill(z_start, z_distortion);             // map of distortion in Z (time)
        }
      }

      // remove electrons outside of our acceptance. Careful though, electrons from just inside 30 cm can contribute in the 1st active layer readout, so leave a little margin
      if (rad_final < min_active_radius - 2.0 || rad_final > max_active_radius + 1.0)
      { continue; }

      if (Verbosity() > 1000)
      {
        std::cout << "electron " << i << " g4hitid " << hiter->first << " f " << f << std::endl;
        std::cout << "radstart " << radstart << " x_start: " << x_start
                  << ", y_start: " << y_start
                  << ",z_start: " << z_start
                  << " t_start " << t_start
                  << " t_path " << t_path
                  << " t_sigma " << t_sigma
                  << " rantime " << rantime
                  << std::endl;

        std::cout << "       rad_final " << rad_final << " x_final " << x_final
          << " y_final " << y_final
          << " z_final " << z_final << " t_final " << t_final 
          << " zdiff " << z_final - z_start << std::endl;
      }

      if (Verbosity() > 0)
      {
        assert(nt);
        nt->Fill(ihit, t_start, t_final, t_sigma, rad_final, z_start, z_final);
      }
      // this fills the cells and updates the hits in temp_hitsetcontainer for this drifted electron hitting the GEM stack
      padplane->MapToPadPlane(truth_clusterer, single_hitsetcontainer.get(),
          temp_hitsetcontainer.get(), hittruthassoc, x_final, y_final, t_final,
          side, hiter, ntpad, nthit);
    }  // end loop over electrons for this g4hit

    TrkrHitSetContainer::ConstRange single_hitset_range = single_hitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
    for (TrkrHitSetContainer::ConstIterator single_hitset_iter = single_hitset_range.first;
         single_hitset_iter != single_hitset_range.second;
         ++single_hitset_iter)
    {
      // we have an itrator to one TrkrHitSet for the Tpc from the single_hitsetcontainer
      TrkrDefs::hitsetkey node_hitsetkey = single_hitset_iter->first;
      const unsigned int layer = TrkrDefs::getLayer(node_hitsetkey);
      const int sector = TpcDefs::getSectorId(node_hitsetkey);
      const int side = TpcDefs::getSide(node_hitsetkey);

      if (Verbosity() > 2)
        std::cout << " hitsetkey " << node_hitsetkey << " layer " << layer << " sector " << sector << " side " << side << std::endl;
      // get all of the hits from the single hitset
      TrkrHitSet::ConstRange single_hit_range = single_hitset_iter->second->getHits();
      for (TrkrHitSet::ConstIterator single_hit_iter = single_hit_range.first;
           single_hit_iter != single_hit_range.second;
           ++single_hit_iter)
      {
        TrkrDefs::hitkey single_hitkey = single_hit_iter->first;

        // Add the hit-g4hit association
        // no need to check for duplicates, since the hit is new
        hittruthassoc->addAssoc(node_hitsetkey, single_hitkey, hiter->first);
        if (Verbosity() > 100)
          std::cout << "        adding assoc for node_hitsetkey " << node_hitsetkey << " single_hitkey " << single_hitkey << " g4hitkey " << hiter->first << std::endl;
      }
    }

    // Dump the temp_hitsetcontainer to the node tree and reset it
    //    - after every "dump_interval" g4hits
    //    - if this is the last g4hit
    if (dump_counter >= dump_interval || count_g4hits == g4hit->size())
    {
      //std::cout << " dump_counter " << dump_counter << " count_g4hits " << count_g4hits << std::endl;

      double eg4hit = 0.0;
      TrkrHitSetContainer::ConstRange temp_hitset_range = temp_hitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
      for (TrkrHitSetContainer::ConstIterator temp_hitset_iter = temp_hitset_range.first;
           temp_hitset_iter != temp_hitset_range.second;
           ++temp_hitset_iter)
      {
        // we have an itrator to one TrkrHitSet for the Tpc from the temp_hitsetcontainer
        TrkrDefs::hitsetkey node_hitsetkey = temp_hitset_iter->first;
        const unsigned int layer = TrkrDefs::getLayer(node_hitsetkey);
        const int sector = TpcDefs::getSectorId(node_hitsetkey);
        const int side = TpcDefs::getSide(node_hitsetkey);
        if (Verbosity() > 100)
          std::cout << "PHG4TpcElectronDrift: temp_hitset with key: " << node_hitsetkey << " in layer " << layer
                    << " with sector " << sector << " side " << side << std::endl;

        // find or add this hitset on the node tree
        TrkrHitSetContainer::Iterator node_hitsetit = hitsetcontainer->findOrAddHitSet(node_hitsetkey);

        // get all of the hits from the temporary hitset
        TrkrHitSet::ConstRange temp_hit_range = temp_hitset_iter->second->getHits();
        for (TrkrHitSet::ConstIterator temp_hit_iter = temp_hit_range.first;
             temp_hit_iter != temp_hit_range.second;
             ++temp_hit_iter)
        {
          TrkrDefs::hitkey temp_hitkey = temp_hit_iter->first;
          TrkrHit *temp_tpchit = temp_hit_iter->second;
          if (Verbosity() > 10 && layer == print_layer)
          {
            std::cout << "      temp_hitkey " << temp_hitkey << " layer " << layer << " pad " << TpcDefs::getPad(temp_hitkey)
                      << " z bin " << TpcDefs::getTBin(temp_hitkey)
                      << "  energy " << temp_tpchit->getEnergy() << " eg4hit " << eg4hit << std::endl;

            eg4hit += temp_tpchit->getEnergy();
	    //            ecollectedhits += temp_tpchit->getEnergy();
            ncollectedhits++;
          }

          // find or add this hit to the node tree
          TrkrHit *node_hit = node_hitsetit->second->getHit(temp_hitkey);
          if (!node_hit)
          {
            // Otherwise, create a new one
            node_hit = new TrkrHitv2();
            node_hitsetit->second->addHitSpecificKey(temp_hitkey, node_hit);
          }

          // Either way, add the energy to it
          node_hit->addEnergy(temp_tpchit->getEnergy());

        }  // end loop over temp hits

        if (Verbosity() > 100 && layer == print_layer)
          std::cout << "  ihit " << ihit << " collected energy = " << eg4hit << std::endl;

      }  // end loop over temp hitsets

      // erase all entries in the temp hitsetcontainer
      temp_hitsetcontainer->Reset();

      // reset the dump counter
      dump_counter = 0;
    }  // end copy of temp hitsetcontainer to node tree hitsetcontainer

    ++ihit;

    single_hitsetcontainer->Reset();

  }  // end loop over g4hits

  truth_clusterer->cluster_and_reset(/*argument is if to reset hitsetkey as well*/ true);

  if (Verbosity() > 2)
  {
    std::cout << "From PHG4TpcElectronDrift: hitsetcontainer printout at end:" << std::endl;
    // We want all hitsets for the Tpc
    TrkrHitSetContainer::ConstRange hitset_range = hitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
    for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first;
         hitset_iter != hitset_range.second;
         ++hitset_iter)
    {
      // we have an itrator to one TrkrHitSet for the Tpc from the trkrHitSetContainer
      TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
      const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
      if (layer != print_layer) continue;
      const int sector = TpcDefs::getSectorId(hitsetkey);
      const int side = TpcDefs::getSide(hitsetkey);

      std::cout << "PHG4TpcElectronDrift: hitset with key: " << hitsetkey << " in layer " << layer << " with sector " << sector << " side " << side << std::endl;

      // get all of the hits from this hitset
      TrkrHitSet *hitset = hitset_iter->second;
      TrkrHitSet::ConstRange hit_range = hitset->getHits();
      for (TrkrHitSet::ConstIterator hit_iter = hit_range.first;
           hit_iter != hit_range.second;
           ++hit_iter)
      {
        TrkrDefs::hitkey hitkey = hit_iter->first;
        TrkrHit *tpchit = hit_iter->second;
        std::cout << "      hitkey " << hitkey << " pad " << TpcDefs::getPad(hitkey) << " z bin " << TpcDefs::getTBin(hitkey)
                  << "  energy " << tpchit->getEnergy() << std::endl;
      }
    }
  }

  if (Verbosity() > 1000)
  {
    std::cout << "From PHG4TpcElectronDrift: hittruthassoc dump:" << std::endl;
    hittruthassoc->identify();

    hittruthassoc->identify();
  }

  ++event_num;  // if doing more than one event, event_num will be incremented.

  if (Verbosity() > 500) 
  {
    std::cout << " TruthTrackContainer results at end of event in PHG4TpcElectronDrift::process_event " << std::endl;
    truthtracks->identify();
  }

  if (Verbosity()>800) {
    truth_clusterer->print(truthtracks);
    truth_clusterer->print_file(truthtracks,"drift_clusters.txt");
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TpcElectronDrift::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 0)
  {
    assert(m_outf);
    assert(nt);
    assert(ntpad);
    assert(nthit);
    assert(ntfinalhit);

    m_outf->cd();
    nt->Write();
    ntpad->Write();
    nthit->Write();
    ntfinalhit->Write();
    m_outf->Close();
  }
  if (do_ElectronDriftQAHistos)
  {
    EDrift_outf.reset(new TFile("ElectronDriftQA.root", "recreate"));
    EDrift_outf->cd();
    deltar->Write();
    deltaphi->Write();
    deltaz->Write();
    deltarnodist->Write();
    deltaphinodist->Write();
    deltarnodiff->Write();
    deltaphinodiff->Write();
    deltaRphinodiff->Write();
    deltaphivsRnodiff->Write();
    hitmapstart->Write();
    hitmapend->Write();
    z_startmap->Write();
    EDrift_outf->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TpcElectronDrift::set_seed(const unsigned int seed)
{
  gsl_rng_set(RandomGenerator.get(), seed);
}

void PHG4TpcElectronDrift::SetDefaultParameters()
{
  // Data on gasses @20 C and 760 Torr from the following source:
  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
  // diffusion and drift velocity for 400kV for NeCF4 50/50 from calculations:
  // http://skipper.physics.sunysb.edu/~prakhar/tpc/HTML_Gases/split.html
  static constexpr double Ne_dEdx = 1.56;   // keV/cm
  static constexpr double CF4_dEdx = 7.00;  // keV/cm
  // double Ne_NPrimary = 12;    // Number/cm
  // double CF4_NPrimary = 51;   // Number/cm
  static constexpr double Ne_NTotal = 43;    // Number/cm
  static constexpr double CF4_NTotal = 100;  // Number/cm
  static constexpr double Tpc_NTot = 0.5 * Ne_NTotal + 0.5 * CF4_NTotal;
  static constexpr double Tpc_dEdx = 0.5 * Ne_dEdx + 0.5 * CF4_dEdx;
  static constexpr double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  set_default_double_param("diffusion_long", 0.012);   // cm/SQRT(cm)
  set_default_double_param("diffusion_trans", 0.004);  // cm/SQRT(cm)
  set_default_double_param("electrons_per_gev", Tpc_ElectronsPerKeV * 1000000.);
  set_default_double_param("min_active_radius", 30.);        // cm
  set_default_double_param("max_active_radius", 78.);        // cm
  set_default_double_param("drift_velocity", 8.0 / 1000.0);  // cm/ns
  set_default_double_param("max_time", 13200.);              //ns
  set_default_double_param("extended_readout_time", 7000.);  //ns

  // These are purely fudge factors, used to increase the resolution to 150 microns and 500 microns, respectively
  // override them from the macro to get a different resolution
  set_default_double_param("added_smear_trans", 0.085);  // cm
  set_default_double_param("added_smear_long", 0.105);   // cm

  return;
}

PHG4TpcElectronDrift::~PHG4TpcElectronDrift() {
  if (truth_clusterer != nullptr) delete truth_clusterer;
}

void PHG4TpcElectronDrift::setTpcDistortion(PHG4TpcDistortion *distortionMap)
{
  m_distortionMap.reset(distortionMap);
}

void PHG4TpcElectronDrift::registerPadPlane(PHG4TpcPadPlane *inpadplane)
{
  if (Verbosity()) std::cout << "Registering padplane " << std::endl;
  padplane.reset(inpadplane);
  padplane->Detector(Detector());
  padplane->UpdateInternalParameters();
  if (Verbosity()) std::cout << "padplane registered and parameters updated" << std::endl;

  return;
}
