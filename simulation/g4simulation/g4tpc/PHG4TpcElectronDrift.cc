// this is the new containers version
// it uses the same MapToPadPlane as the old containers version

#include "PHG4TpcElectronDrift.h"
#include "PHG4TpcPadPlane.h"                            // for PHG4TpcPadPlane
#include "DistortedTrackContainer.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrDefs.h>

#include <tpc/TpcDefs.h>
#include <tpc/TpcHit.h>

#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>
#include <phparameter/PHParameterInterface.h>           // for PHParameterIn...


#include <pdbcalbase/PdbParameterMapContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>                         // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>                           // for PHDataNode
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                               // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                             // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>                                // for PHWHERE

#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TNtuple.h>
#include <TSystem.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>                                // for gsl_rng_alloc

#include <bitset>
#include <cassert>
#include <cmath>                                       // for sqrt, fabs, NAN
#include <cstdlib>                                     // for exit
#include <iostream>
#include <map>                                          // for _Rb_tree_cons...
#include <utility>                                      // for pair

using namespace std;

//_____________________________________________________________
namespace
{
  template<class T> inline constexpr T square(const T& x) { return x*x; }
}


//_____________________________________________________________________
void PHG4TpcElectronDrift::DistortedTrackListClear()
{
  m_container->clearDistortions();
}

//_____________________________________________________________
PHG4TpcElectronDrift::PHG4TpcElectronDrift(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
  // this is used as a buffer for charge collection from a single g4hit
  , temp_hitsetcontainer(new TrkrHitSetContainer())
{
  InitializeParameters();
  RandomGenerator.reset(gsl_rng_alloc(gsl_rng_mt19937));
  set_seed(PHRandomSeed());
  return;
}

//_____________________________________________________________
int PHG4TpcElectronDrift::Init(PHCompositeNode *topNode)
{
  padplane->Init(topNode);
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
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  string paramnodename = "G4CELLPARAM_" + detector;
  string geonodename = "G4CELLPAR_" + detector;
  string tpcgeonodename = "G4GEO_" + detector;
  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
  {
    cout << Name() << " Could not locate G4HIT node " << hitnodename << endl;
    topNode->print();
    gSystem->Exit(1);
    exit(1);
  }

  // new containers
  hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if(!hitsetcontainer)
    {
      PHNodeIterator dstiter(dstNode);
      PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
      if (!DetNode)
  {
    DetNode = new PHCompositeNode("TRKR");
    dstNode->addNode(DetNode);
  }

      hitsetcontainer = new TrkrHitSetContainer();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(hitsetcontainer, "TRKR_HITSET", "PHObject");
      DetNode->addNode(newNode);
    }

  hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");
  if(!hittruthassoc)
    {
      PHNodeIterator dstiter(dstNode);
      PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
      if (!DetNode)
  {
    DetNode = new PHCompositeNode("TRKR");
    dstNode->addNode(DetNode);
  }

      hittruthassoc = new TrkrHitTruthAssoc();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(hittruthassoc, "TRKR_HITTRUTHASSOC", "PHObject");
      DetNode->addNode(newNode);
    }

  seggeonodename = "CYLINDERCELLGEOM_SVTX";  // + detector;
  PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, seggeonodename.c_str());
  if (!seggeo)
  {
    seggeo = new PHG4CylinderCellGeomContainer();
    PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(seggeo, seggeonodename.c_str(), "PHObject");
    runNode->addNode(newNode);
  }

  // Local evaluators
  if( m_enable_distortions )
  {
    // get EVAL node
    iter = PHNodeIterator(dstNode);
    auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
    if( !evalNode )
    {
      // create
      std::cout << "PHG4TpcElectronDrift::InitRun - EVAL node missing - creating" << std::endl;
      evalNode = new PHCompositeNode( "EVAL" );
      dstNode->addNode(evalNode);
    }

    //auto newNode = new PHIODataNode<PHObject>( new Container, "PHG4TpcElectronDrift::Container","PHObject");
    auto newNode = new PHIODataNode<PHObject>( new DistortedTrackContainer, "DistortedTrackContainer","PHObject");
    evalNode->addNode(newNode);
  }

  UpdateParametersWithMacro();
  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(detector);
    runNode->addNode(RunDetNode);
  }
  SaveToNodeTree(RunDetNode, paramnodename);

  // save this to the parNode for use
  PHNodeIterator parIter(parNode);
  PHCompositeNode *ParDetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", detector));
  if (!ParDetNode)
  {
    ParDetNode = new PHCompositeNode(detector);
    parNode->addNode(ParDetNode);
  }
  PutOnParNode(ParDetNode, geonodename);

  // find Tpc Geo
  PHNodeIterator tpcpariter(ParDetNode);
  PHParametersContainer *tpcparams = findNode::getClass<PHParametersContainer>(ParDetNode, tpcgeonodename);
  if (!tpcparams)
  {
    string runparamname = "G4GEOPARAM_" + detector;
    PdbParameterMapContainer *tpcpdbparams = findNode::getClass<PdbParameterMapContainer>(RunDetNode, runparamname);
    if (tpcpdbparams)
    {
      tpcparams = new PHParametersContainer(detector);
      if (Verbosity()) tpcpdbparams->print();
      tpcparams->CreateAndFillFrom(tpcpdbparams, detector);
      ParDetNode->addNode(new PHDataNode<PHParametersContainer>(tpcparams, tpcgeonodename));
    }
    else
    {
      cout << "PHG4TpcElectronDrift::InitRun - failed to find " << runparamname << " in order to initialize " << tpcgeonodename << ". Aborting run ..." << endl;
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
  min_time = 0.0;
  max_time = (tpc_length / 2.) / drift_velocity;
  electrons_per_gev = get_double_param("electrons_per_gev");
  min_active_radius = get_double_param("min_active_radius");
  max_active_radius = get_double_param("max_active_radius");

  // space charge distortions from external file
  if( m_enable_distortions )
  {
    std::cout << "PHG4TpcElectronDrift::InitRun - reading distortions from " << m_distortion_filename << std::endl;
    m_distortion_tfile = TFile::Open( m_distortion_filename.c_str());
    if( !m_distortion_tfile )
    {
      std::cout << "PHG4TpcElectronDrift::InitRun - cannot open " << m_distortion_filename << std::endl;
      exit(1);
    }

    // Open TH3F files only once that contain distortions due to space charge
    hDPint= dynamic_cast<TH3*>(m_distortion_tfile->Get("hIntDistortionP")); assert( hDPint );
    hDRint= dynamic_cast<TH3*>(m_distortion_tfile->Get("hIntDistortionR")); assert( hDRint );
    hDZint= dynamic_cast<TH3*>(m_distortion_tfile->Get("hIntDistortionZ")); assert( hDZint );

    // coordinates
    std::cout << "PHG4TpcElectronDrift::InitRun - coordinates: " << std::bitset<3>(m_coordinates) << std::endl;

    // dump axis limits
    for(const auto& axis:{ hDPint->GetXaxis(), hDPint->GetYaxis(), hDPint->GetZaxis() })
    { std::cout << "PHG4TpcElectronDrift::InitRun - axis: " << axis->GetTitle() << " bins: " << axis->GetNbins() << " limits: " << axis->GetXmin() << " " << axis->GetXmax() << std::endl; }

  }

  auto se = Fun4AllServer::instance();
  dlong = new TH1F("difflong", "longitudinal diffusion", 100, diffusion_long - diffusion_long / 2., diffusion_long + diffusion_long / 2.);
  se->registerHisto(dlong);
  dtrans = new TH1F("difftrans", "transversal diffusion", 100, diffusion_trans - diffusion_trans / 2., diffusion_trans + diffusion_trans / 2.);
  se->registerHisto(dtrans);

  if (Verbosity())
  {
    // eval tree only when verbosity is on
    m_outf = new TFile("nt_out.root", "recreate");
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

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TpcElectronDrift::process_event(PHCompositeNode *topNode)
{
  unsigned int print_layer = 18;

  // container for local evaluations
  //m_container = findNode::getClass<Container>(topNode, "PHG4TpcElectronDrift::Container");
  m_container = findNode::getClass<DistortedTrackContainer>(topNode, "DistortedTrackContainer");



  if( m_container ) m_container->Reset();
  else
    if(Verbosity() > 3) std::cout << PHWHERE << " no DistortedTrackContainer on node tree " << std::endl;

  // g4hits
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << "Could not locate g4 hit node " << hitnodename << endl;
      gSystem->Exit(1);
    }

  PHG4HitContainer::ConstIterator hiter;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();

  double ecollectedhits = 0.0;
  int ncollectedhits = 0;
  double ihit = 0;
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      double t0 = fmax(hiter->second->get_t(0), hiter->second->get_t(1));
      if (t0 > max_time)
	{
	  continue;
	}

      // for very high occupancy events, accessing the TrkrHitsets on the node tree for every drifted electron seems to be very slow
      // Instead, use a temporary map to accumulate the charge from all drifted electrons, then copy to the node tree later

      double eion = hiter->second->get_eion();
      unsigned int n_electrons = gsl_ran_poisson(RandomGenerator.get(), eion * electrons_per_gev);
      if (Verbosity() > 100)
	cout << "  new hit with t0, " << t0 << " g4hitid " << hiter->first
	     << " eion " << eion << " n_electrons " << n_electrons
	     << " entry z " << hiter->second->get_z(0) << " exit z " << hiter->second->get_z(1) << " avg z" << (hiter->second->get_z(0) + hiter->second->get_z(1)) / 2.0
	     << endl;

      if (n_electrons == 0) // gsl_ran_poisson returns unsigned int, cannot be < 0
	{
	  continue;
	}

      if (Verbosity() > 100)
	{
	  cout << endl
	       << "electron drift: ihit " << ihit << " g4hit ID " << hiter->first << " created electrons: " << n_electrons
	       << " from " << eion * 1000000 << " keV" << endl;
	  cout << " entry x,y,z = " << hiter->second->get_x(0) << "  " << hiter->second->get_y(0) << "  " << hiter->second->get_z(0)
	       << " radius " << sqrt(pow(hiter->second->get_x(0), 2) + pow(hiter->second->get_y(0), 2)) << endl;
	  cout << " exit x,y,z = " << hiter->second->get_x(1) << "  " << hiter->second->get_y(1) << "  " << hiter->second->get_z(1)
	       << " radius " << sqrt(pow(hiter->second->get_x(1), 2) + pow(hiter->second->get_y(1), 2)) << endl;
	}

      for (unsigned int i = 0; i < n_electrons; i++)
	{
    // We choose the electron starting position at random from a flat distribution along the path length
    // the parameter t is the fraction of the distance along the path betwen entry and exit points, it has values between 0 and 1
    const double f = gsl_ran_flat(RandomGenerator.get(), 0.0, 1.0);

    const double x_start = hiter->second->get_x(0) + f * (hiter->second->get_x(1) - hiter->second->get_x(0));
    const double y_start = hiter->second->get_y(0) + f * (hiter->second->get_y(1) - hiter->second->get_y(0));
    const double z_start = hiter->second->get_z(0) + f * (hiter->second->get_z(1) - hiter->second->get_z(0));
    const double t_start = hiter->second->get_t(0) + f * (hiter->second->get_t(1) - hiter->second->get_t(0));

    const double radstart = std::sqrt(square(x_start) + square(y_start));
    const double r_sigma = diffusion_trans * sqrt(tpc_length / 2. - fabs(z_start));
    double rantrans = gsl_ran_gaussian(RandomGenerator.get(), r_sigma);
    rantrans += gsl_ran_gaussian(RandomGenerator.get(), added_smear_sigma_trans);

    const double t_path = (tpc_length / 2. - fabs(z_start)) / drift_velocity;
    const double t_sigma = diffusion_long * sqrt(tpc_length / 2. - fabs(z_start)) / drift_velocity;
    double rantime = gsl_ran_gaussian(RandomGenerator.get(), t_sigma);
    rantime += gsl_ran_gaussian(RandomGenerator.get(), added_smear_sigma_long) / drift_velocity;

    // drift time
    const double t_final = t_start + t_path + rantime;
    if (t_final < min_time || t_final > max_time) continue;

    double z_final = 0;
    if (z_start < 0) z_final = -tpc_length / 2. + t_final * drift_velocity;
    else z_final = tpc_length / 2. - t_final * drift_velocity;

    double x_final = 0;
    double y_final = 0;
    double rad_final = 0;
    const double ranphi = gsl_ran_flat(RandomGenerator.get(), -M_PI, M_PI);
    if( m_enable_distortions )
    {
      // get starting azimuth angle. Convert to [0,2pi[ to match histograms
      const double z_abs = std::abs( z_start );
      double phistart = std::atan2(y_start,x_start);
      if( phistart < 0 ) phistart += 2*M_PI;

      // add radial distortion
      const double dr = hDRint->Interpolate(phistart,radstart,z_abs);
      rad_final = (m_coordinates&COORD_R) ? radstart+dr : radstart;

      // add azimutal distortion
      const double phi_final = (m_coordinates&COORD_PHI) ?
        phistart + hDPint->Interpolate(phistart,radstart,z_abs)/radstart :
        phistart;

      // add z distortion
      z_final += (m_coordinates&COORD_Z) ? hDZint->Interpolate(phistart,radstart,z_abs) : 0;

      //std::cout << " radstart  " << radstart << " rad_final "  << rad_final <<  " phistart " << phistart << " phi_final " << phi_final << " zstart " << z_start << " z_final " << z_final << std::endl;

      // convert back to cartesian coordinates, add diffusion
      x_final = rad_final*std::cos(phi_final)+rantrans*cos(ranphi);
      y_final = rad_final*std::sin(phi_final)+rantrans*sin(ranphi);

      // fill
      if( m_container )
      {
        DistortionStruct distortion;
        distortion._r = radstart;
        distortion._phi = phistart;
        distortion._z = z_abs;
        distortion._dr = dr;
        distortion._dphi = phi_final - phistart;
        distortion._dz = 0;
        m_container->addDistortion( distortion );
	//std::cout << " Update DistortionStruct: radstart " << distortion._r << " dr " << distortion._dr << " phistart " << distortion._phi << " dphi " << distortion._dphi << " zstart " << distortion._z << " dz " << distortion._dz << std::endl;  
      }
    } else {
      x_final = x_start + rantrans * cos(ranphi);
      y_final = y_start + rantrans * sin(ranphi);
      rad_final = std::sqrt(square(x_final) + square(y_final));
    }

    // remove electrons outside of our acceptance. Careful though, electrons from just inside 30 cm can contribute in the 1st active layer readout, so leave a little margin
    if (rad_final < min_active_radius - 2.0 || rad_final > max_active_radius + 1.0)
    { continue; }

    if (Verbosity() > 1000)
    {
	      cout << "ihit " << ihit << " electron " << i << " g4hitid " << hiter->first << " f " << f << endl;
	      cout << "radstart " << radstart << " x_start: " << x_start
		   << ", y_start: " << y_start
		   << ",z_start: " << z_start
		   << " t_start " << t_start
		   << " t_path " << t_path
		   << " t_sigma " << t_sigma
		   << " rantime " << rantime
		   << endl;

	      //if( sqrt(x_start*x_start+y_start*y_start) > 68.0 && sqrt(x_start*x_start+y_start*y_start) < 72.0)
	      cout << "       rad_final " << rad_final << " x_final " << x_final << " y_final " << y_final
		   << " z_final " << z_final << " t_final " << t_final << " zdiff " << z_final - z_start << endl;
	    }

	  if (Verbosity() > 0)
	  {
	    assert(nt);
      nt->Fill(ihit, t_start, t_final, t_sigma, rad_final, z_start, z_final);
	  }

	  // this fills the cells and updates the hits in temp_hitsetcontainer for this drifted electron hitting the GEM stack
	  MapToPadPlane(x_final, y_final, z_final, hiter, ntpad, nthit);
	}  // end loop over electrons for this g4hit

      if(Verbosity() > 100)
	cout << "Finished drifting " << n_electrons << " electrons from ihit " << ihit
	     << " now process temp_hitsetcontainer  " << endl;

      // transfer the hits from temp_hitsetcontainer to hitsetcontainer on the node tree
      double eg4hit = 0.0;
      TrkrHitSetContainer::ConstRange temp_hitset_range = temp_hitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
      for (TrkrHitSetContainer::ConstIterator temp_hitset_iter = temp_hitset_range.first;
	   temp_hitset_iter != temp_hitset_range.second;
	   ++temp_hitset_iter)
	{

	  // we have an iterator to one TrkrHitSet for the Tpc from the temp_hitsetcontainer
	  TrkrDefs::hitsetkey node_hitsetkey = temp_hitset_iter->first;
	  const unsigned int layer = TrkrDefs::getLayer(node_hitsetkey);
	  const int sector = TpcDefs::getSectorId(node_hitsetkey);
	  const int side = TpcDefs::getSide(node_hitsetkey);
	  if(Verbosity()>100 && layer)
	    cout << "PHG4TpcElectronDrift: temp_hitset with key: " << node_hitsetkey << " in layer " << layer << " with sector " << sector << " side " << side << endl;

	  // find or add this hitset on the node tree
	  TrkrHitSetContainer::Iterator node_hitsetit = hitsetcontainer->findOrAddHitSet(node_hitsetkey);

	  // get all of the hits from the temporary hitset
	  TrkrHitSet::ConstRange temp_hit_range = temp_hitset_iter->second->getHits();
	  for(TrkrHitSet::ConstIterator temp_hit_iter = temp_hit_range.first;
	      temp_hit_iter != temp_hit_range.second;
	      ++temp_hit_iter)
	    {
	      TrkrDefs::hitkey temp_hitkey = temp_hit_iter->first;
	      TrkrHit *temp_tpchit = temp_hit_iter->second;



	      if(Verbosity() > 100 && layer == print_layer)
		{
		  cout << "      temp_hitkey " << temp_hitkey << " l;ayer " << layer << " pad " << TpcDefs::getPad(temp_hitkey)
		       << " z bin " << TpcDefs::getTBin(temp_hitkey)
		       << "  energy " << temp_tpchit->getEnergy() << " eg4hit " << eg4hit << endl;

		  eg4hit +=  temp_tpchit->getEnergy();
		  ecollectedhits +=  temp_tpchit->getEnergy();
		  ncollectedhits++;
		}

	      // find or add this hit to the node tree
	      TrkrHit *node_hit = node_hitsetit->second->getHit(temp_hitkey);
	      if(!node_hit)
		{
		  // Otherwise, create a new one
		  node_hit = new TpcHit();
		  node_hitsetit->second->addHitSpecificKey(temp_hitkey, node_hit);

		  // Add the hit-g4hit association
		  // no need to check for duplicates, since the hit is new
		  hittruthassoc->addAssoc(node_hitsetkey, temp_hitkey, hiter->first);
		} else {
		// Add the hit-g4hit association
		// TODO: check if duplication can happen
		hittruthassoc->findOrAddAssoc(node_hitsetkey, temp_hitkey, hiter->first);
	      }

	      // Either way, add the energy to it
	      node_hit->addEnergy(temp_tpchit->getEnergy());

	    }  // end loop over temp hits

	  if(Verbosity() > 100 && layer == print_layer)
	    cout << "  ihit " << ihit << " collected energy = " << eg4hit << endl;

	} // end loop over temp hitsets

      // erase all entries in the temp hitsetcontainer
      temp_hitsetcontainer->Reset();

      ihit++;

    } // end loop over g4hits

  if(Verbosity() > 2)
    {
      cout << "From PHG4TpcElectronDrift: hitsetcontainer printout at end:" << endl;
      double eallhits = 0.0;
      int nallhits = 0;

      // We want all hitsets for the Tpc
      TrkrHitSetContainer::ConstRange hitset_range = hitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
      for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first;
	   hitset_iter != hitset_range.second;
	   ++hitset_iter)
	{
	  // we have an itrator to one TrkrHitSet for the Tpc from the trkrHitSetContainer
	  TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
	  const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
	  if(layer != print_layer)  continue;
	  const int sector = TpcDefs::getSectorId(hitsetkey);
	  const int side = TpcDefs::getSide(hitsetkey);

	  cout << "PHG4TpcElectronDrift: hitset with key: " << hitsetkey << " in layer " << layer << " with sector " << sector << " side " << side << endl;

	  // get all of the hits from this hitset
	  TrkrHitSet *hitset = hitset_iter->second;
	  TrkrHitSet::ConstRange hit_range = hitset->getHits();
	  for(TrkrHitSet::ConstIterator hit_iter = hit_range.first;
	      hit_iter != hit_range.second;
	      ++hit_iter)
	    {
	      TrkrDefs::hitkey hitkey = hit_iter->first;
	      TrkrHit *tpchit = hit_iter->second;

	      if(layer == print_layer)
		{
		  nallhits++;
		  eallhits +=  tpchit->getEnergy();
		  cout << "      hitkey " << hitkey << " pad " << TpcDefs::getPad(hitkey) << " z bin " << TpcDefs::getTBin(hitkey)
		       << "  energy " << tpchit->getEnergy() << endl;

		  ntfinalhit->Fill(layer, TpcDefs::getPad(hitkey), TpcDefs::getTBin(hitkey), tpchit->getEnergy());
		}
	    }
	}

      cout << " eallhits = " << eallhits << " nallhits " << nallhits << " for print_layer " << print_layer
	   << " ecollectedhits = " << ecollectedhits << " ncollectedhits " << ncollectedhits << endl;
    }


  if(Verbosity() > 1000)
    {
      cout << "From PHG4TpcElectronDrift: hittruthassoc dump:" << endl;
      hittruthassoc->identify();
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TpcElectronDrift::MapToPadPlane(const double x_gem, const double y_gem, const double t_gem, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit)
{ padplane->MapToPadPlane(temp_hitsetcontainer.get(), hittruthassoc, x_gem, y_gem, t_gem, hiter, ntpad, nthit); }

int PHG4TpcElectronDrift::End(PHCompositeNode *topNode)
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

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TpcElectronDrift::set_seed(const unsigned int seed)
{ gsl_rng_set(RandomGenerator.get(), seed); }

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
  static constexpr double Tpc_NTot = 0.5*Ne_NTotal + 0.5*CF4_NTotal;
  static constexpr double Tpc_dEdx = 0.5*Ne_dEdx + 0.5*CF4_dEdx;
  static constexpr double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  set_default_double_param("diffusion_long", 0.012);   // cm/SQRT(cm)
  set_default_double_param("diffusion_trans", 0.004);  // cm/SQRT(cm)
  set_default_double_param("electrons_per_gev", Tpc_ElectronsPerKeV * 1000000.);
  set_default_double_param("min_active_radius", 30.);        // cm
  set_default_double_param("max_active_radius", 78.);        // cm
  set_default_double_param("drift_velocity", 8.0 / 1000.0);  // cm/ns

  // These are purely fudge factors, used to increase the resolution to 150 microns and 500 microns, respectively
  // override them from the macro to get a different resolution
  set_default_double_param("added_smear_trans", 0.085);  // cm
  set_default_double_param("added_smear_long", 0.105);   // cm

  return;
}

void PHG4TpcElectronDrift::registerPadPlane(PHG4TpcPadPlane *inpadplane)
{
  cout << "Registering padplane " << endl;
  padplane.reset(inpadplane);
  padplane->Detector(Detector());
  padplane->UpdateInternalParameters();
  cout << "padplane registered and parameters updated" << endl;

  return;
}
