// this is the new containers version
// it uses the same MapToPadPlane as the old containers version

#include "PHG4TpcElectronDrift.h"
#include "PHG4TpcPadPlane.h"                            // for PHG4TpcPadPlane

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
#include <TH2.h>
#include <TH3F.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TNtuple.h>
#include <TSystem.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>                                // for gsl_rng_alloc

#include <cassert>
#include <cmath>                                       // for sqrt, fabs, NAN
#include <cstdlib>                                     // for exit
#include <iostream>
#include <map>                                          // for _Rb_tree_cons...
#include <utility>                                      // for pair

using namespace std;

PHG4TpcElectronDrift::PHG4TpcElectronDrift(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , hitsetcontainer(nullptr)
  , temp_hitsetcontainer(new TrkrHitSetContainer())// this is used as a buffer for charge collection from a single g4hit
  , hittruthassoc(nullptr)
  , padplane(nullptr)
  , dlong(nullptr)
  , dtrans(nullptr)
  , m_outf(nullptr)
  , nt(nullptr)
  , nthit(nullptr)
  , ntfinalhit(nullptr)
  , ntpad(nullptr)
  , diffusion_trans(NAN)
  , diffusion_long(NAN)
  , drift_velocity(NAN)
  , electrons_per_gev(NAN)
  , min_active_radius(NAN)
  , max_active_radius(NAN)
  , min_time(NAN)
  , max_time(NAN)
{
  //cout << "Constructor of PHG4TpcElectronDrift" << endl;
  InitializeParameters();
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  set_seed(PHRandomSeed());  // fixed seed is handled in this funtcion
  return;
}

PHG4TpcElectronDrift::~PHG4TpcElectronDrift()
{
  gsl_rng_free(RandomGenerator);
  delete padplane;
  delete temp_hitsetcontainer;
}
  
  int PHG4TpcElectronDrift::Init(PHCompositeNode *topNode)
{
  padplane->Init(topNode);
  event_num = 0;
  if(event_num != 0)
    {
      cout << "somehow event_num is non-zero" << endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}
	      
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
  
  
  TFile *StaticDistFile=new TFile("zSmooth.50kHz.flat_B1.4_E-400.0.ross_phislice_lookup_r16xp36xz40.distortion_map.hist.root");//includes Trees of TH3Fs                                     
  if(StaticDistFile->GetSize() == -1)
    {
      cout << "Distortion file could not be opened!" << endl;
    }
 
  //Open Static Space Charge Maps
  hDPint=(TH3F*)StaticDistFile->Get("hIntDistortionP"); // Open TH3F files only once that contain distortions due to space charge
  hDRint=(TH3F*)StaticDistFile->Get("hIntDistortionR");
  hDZint=(TH3F*)StaticDistFile->Get("hIntDistortionZ");
  hDPdiff=(TH3F*)StaticDistFile->Get("hDistortionP"); // Open TH3F files only once that contain distortions due to space charge
  hDRdiff=(TH3F*)StaticDistFile->Get("hDistortionR");
  hDZdiff=(TH3F*)StaticDistFile->Get("hDistortionZ");
  cout << "no segfault after staticDistFile opening" << endl;
  
  //Open time ordered maps
  TFile *TimeDistFile=new TFile("TimeOrderedDistortions.root");//includes Trees of TH3Fs                                                                                               
  if(TimeDistFile->GetSize() == -1)
    {
      cout << "TimeOrderedDistortion file could not be opened!" << endl;
    }
  TimehDR=0; // Initialize branches 
  TimehDP=0;
  TimehDZ=0;
  TimeInthDR=0;
  TimeInthDP=0;
  TimeInthDZ=0;
  TimeTree = (TTree*)TimeDistFile->Get("TimeDists");
  TimeTree->SetBranchAddress("hDistortionR",&TimehDR);
  TimeTree->SetBranchAddress("hDistortionP",&TimehDP);
  TimeTree->SetBranchAddress("hDistortionZ",&TimehDZ);
  TimeTree->SetBranchAddress("hIntDistortionR",&TimeInthDR);
  TimeTree->SetBranchAddress("hIntDistortionP",&TimeInthDP);
  TimeTree->SetBranchAddress("hIntDistortionZ",&TimeInthDZ);
  cout << "no segfault after TimeDistFile opening" << endl;
  //
  
  Fun4AllServer *se = Fun4AllServer::instance();
  bool do_Centralmem = true;
  if(do_Centralmem)
    { 
      double x_start,y_start,x_final,y_final;
      TFile *CMFile=new TFile("Centralmem.root");//includes TGraph
      if(CMFile->GetSize() == -1)
	{
	  cout << "CM in file could not be opened!" << endl;
	}
      CM_outf = new TFile("CM_PositionTree.root", "recreate");
      if(CM_outf->GetSize() == -1)
	{
	  cout << "CM out file could not be opened!" << endl;
	}
      
      CM=(TGraph*)CMFile->Get("Graph");
      CMTimeDists = new TTree("CMTimeDists","Tree w/ event_num, x_i, y_i, x_f, y_f");
      CMTimeDists->Branch("event_num",&event_num,"event_num/I");
      CMTimeDists->Branch("x_start",&x_start,"x_start/B");
      CMTimeDists->Branch("y_start",&y_start,"y_start/B");
      CMTimeDists->Branch("x_final",&x_final,"x_final/B");
      CMTimeDists->Branch("y_final",&y_final,"y_final/B");
    }
  
  if(true) //Set verbosity later
    {
      //Add all diagnostic histograms
      dlong = new TH1F("difflong", "longitudinal diffusion", 100, diffusion_long - diffusion_long / 2., diffusion_long + diffusion_long / 2.);
      se->registerHisto(dlong);
      dtrans = new TH1F("difftrans", "transversal diffusion", 100, diffusion_trans - diffusion_trans / 2., diffusion_trans + diffusion_trans / 2.);
      se->registerHisto(dtrans);
      EDrift_outf = new TFile("ElectronDriftQA.root", "recreate");
      if(EDrift_outf->GetSize() == -1)
	{
	  cout << "EDriftQA file could not be opened!" << endl;
	}
      hitmapstart = new TH2F("hitmapstart","g4hit starting X-Y locations",1560,-78,78,1560,-78,78);
      se->registerHisto(hitmapstart);
      hitmapend = new TH2F("hitmapend","g4hit final X-Y locations",1560,-78,78,1560,-78,78);
      se->registerHisto(hitmapend);
      z_startmap = new TH2F("z_startmap","g4hit starting Z vs. R locations",2000,-100,100,780,0,78);
      se->registerHisto(z_startmap);
      deltaphi = new TH2F("deltaphi","Total delta phi; phi (rad);#Delta phi (rad)",600,0,2*M_PI,1000,-.01,.01);
      se->registerHisto(deltaphi);   
      deltaphidiff = new TH2F("deltaphidiff","Total delta phi; phi (rad);#Delta phi (rad)",600,0,2*M_PI,1000,-.07,.07);
      se->registerHisto(deltaphidiff); 
      deltaphidifference = new TH2F("deltaphidifference","Total delta phi; phi (rad);#Delta phi (rad)",600,0,2*M_PI,1000,-.01,.01);
      se->registerHisto(deltaphidifference); 
      deltaphidifferencepercent = new TH2F("deltaphidifferencepercent","Total delta phi percent; phi (rad);Delta phi percent",600,0,2*M_PI,1000,-20,20);
      se->registerHisto(deltaphidifferencepercent); 
      deltaphiint = new TH2F("deltaphiint","Total delta phi; phi (rad);#Delta phi (rad)",600,0,2*M_PI,1000,-.07,.07);
      se->registerHisto(deltaphiint); 
      deltatime = new TH1F("deltatime","Total time difference between integrated and differential runtimes per G4hit in us",300,-15,15);
      se->registerHisto(deltatime);
      deltaz = new TH2F("deltaz","Total delta z; z (cm);#Delta z (cm)",1000,0,100,1000,-.5,.5);
      se->registerHisto(deltaz); 
      deltaphinodiff = new TH2F("deltaphinodiff","Total delta phi; phi (rad);#Delta phi (rad)",600,0,2*M_PI,1000,-.01,.01);
      se->registerHisto(deltaphinodiff); 
      deltaphinodist = new TH2F("deltaphinodist","Total delta phi; phi (rad);#Delta phi (rad)",600,0,2*M_PI,1000,-.1,.1);
      se->registerHisto(deltaphinodist); 
      deltar = new TH2F("deltar","Total Delta r; r (cm);#Delta r (cm)",580,20,78,1000,-.01,.01);
      se->registerHisto(deltar);
      deltardiff = new TH2F("deltardiff","Total Delta r; r (cm);#Delta r (cm)",580,20,78,1000,-6,6);
      se->registerHisto(deltardiff);
      deltardifference = new TH2F("deltardifference","Total Delta r; r (cm);#Delta r (cm)",580,20,78,1000,-.005,.005);
      se->registerHisto(deltardifference);
      deltardifferencepercent = new TH2F("deltardifferencepercent","Total Delta r; r (cm);% difference in #Delta r between integrated and differential",580,20,78,1000,-20,20);
      se->registerHisto(deltardifferencepercent);
      deltarint = new TH2F("deltarint","Total Delta r; r (cm);#Delta r (cm)",580,20,78,1000,-6,6);
      se->registerHisto(deltarint);
      deltarnodiff = new TH2F("deltarnodiff","Total Delta r; r (cm);#Delta r (cm)",580,20,78,1000,-.01,.01);
      se->registerHisto(deltarnodiff);
      deltarnodist = new TH2F("deltarnodist","Total Delta r; r (cm);#Delta r (cm)",580,20,78,1000,-6,6);
      se->registerHisto(deltarnodist);
     }
     
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
int PHG4TpcElectronDrift::DistortionIntegral(double radstart, double phistart, double z_start, double* radaddress, double* phiaddress)
{
  double distrp=0;
  double distr=0;
  double newp=0;
  double newr=0;

  Int_t binp = TimehDP->GetXaxis()->FindBin(phistart+M_PI);
  Int_t binr = TimehDR->GetYaxis()->FindBin(radstart);
  z_startmap->Fill(z_start,radstart);

  int nBinZ=TimehDR->GetNbinsZ();// Only need to get number of bins in Z, but nice to see them all as a sanity check                                                                                     
  int nBinR=TimehDR->GetNbinsY();
  int nBinP=TimehDR->GetNbinsX();

  if(Verbosity() > 100){
    cout << "total nbins r = " << nBinR << " total nbins p = " << nBinP << " total nbins z = " << nBinZ << endl;
  }
  newr=radstart;// placeholders for comparing original vs distorted values                                                                                                                                 
  newp=phistart+M_PI;
  double newz=abs(z_start); //z_start not positive definite                       
  for(Int_t binz = TimehDR->GetZaxis()->FindBin(newz); binz<=TimehDR->GetNbinsZ()-1; binz++) // Sum over steps in Z starting at some Z specified from process_event, additional complication due to TH3 having extr\a bins on outside for interpolation      
  {
    if (binz==1)
      {
	binz=2; // Eliminates double counting of first z bin in case findbin puts into bin 1 instead of 2 (since they contain the same values)                                                            
      }
    
    distr = TimehDR->Interpolate(newp,newr,binz);//+TimehDR->Interpolate(newp,newr,binz); //If any of these goes into bin 1, this gives an error due to bin "below" being non-existent, for phi need to implement periodic BCs                        
    distrp = TimehDP->Interpolate(newp,newr,binz);//+TimehDP->Interpolate(newp,newr,binz);

    newr += distr; // Sum over r-distortions                                                                                                                                                              
    newp += distrp/newr; // Sum over phi distortions (distrp is r*phi, need to divide by r to get phi)                                                                                                    

        if(Verbosity() > 100)
     {
	cout << "distortion in rphi is " << distrp << " distortion in r is " << distr << " bin number of phi,r,z is " << binp << " " << binr << " " << binz << endl;
     }
  }
 
  if(Verbosity() > 100){
    cout << "inside DistortionIntegral, starting positions r(20-78cm),p(0-2*pi rad),z(0-100cm) = " << radstart <<" "<< phistart+M_PI <<" " << z_start << endl;
    cout << "bin containing phi = 0 is " <<  TimehDP->GetXaxis()->FindBin(0.0001) << "bin containing phi = 2*pi " <<  TimehDP->GetXaxis()->FindBin(6.28) << endl;
    cout << "final p = " << newp << endl;
    cout << "delta p =" << newp-(phistart+M_PI) << endl;
    cout << "final r = " << newr << endl;
    cout << "delta r =" << newr-radstart << endl;
  }


  *radaddress = newr;
  *phiaddress = newp; // using addresses to allow passing more than one value from inside the function                                                                                                      

  return 0;
}

int PHG4TpcElectronDrift::process_event(PHCompositeNode *topNode)
{

  //Get SC maps corresponding to event_num
  //cout << "event num is " << event_num << endl;
  TimeTree->GetEntry(event_num);
  //
  //Declare stopwatches for time performance characterization
  TStopwatch *diffwatch = new TStopwatch();
  TStopwatch *intwatch = new TStopwatch();
  TStopwatch *totwatch = new TStopwatch();
  intwatch->Stop();
  diffwatch->Stop();     
  //
  
  
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
  if(do_Centralmem)
    {
      hit_begin_end.second = std::next(hit_begin_end.first);
    }
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
    unsigned int n_electrons = gsl_ran_poisson(RandomGenerator, eion * electrons_per_gev);
    if (Verbosity() > 100)
      cout << "  new hit with t0, " << t0 << " g4hitid " << hiter->first
           << " eion " << eion << " n_electrons " << n_electrons
           << " entry z " << hiter->second->get_z(0) << " exit z " << hiter->second->get_z(1) << " avg z" << (hiter->second->get_z(0) + hiter->second->get_z(1)) / 2.0
           << endl;

    if (n_electrons == 0)
      {
	continue;
      }
    
    if (Verbosity() > 100)
      {
        cout << endl
             << "electron drift: g4hit " << hiter->first << " created electrons: " << n_electrons
             << " from " << eion * 1000000 << " keV" << endl;
        cout << " entry x,y,z = " << hiter->second->get_x(0) << "  " << hiter->second->get_y(0) << "  " << hiter->second->get_z(0)
             << " radius " << sqrt(pow(hiter->second->get_x(0), 2) + pow(hiter->second->get_y(0), 2)) << endl;
        cout << " exit x,y,z = " << hiter->second->get_x(1) << "  " << hiter->second->get_y(1) << "  " << hiter->second->get_z(1)
             << " radius " << sqrt(pow(hiter->second->get_x(1), 2) + pow(hiter->second->get_y(1), 2)) << endl;
    }
  if(do_Centralmem)
    {
      n_electrons = 7668;
    }

    for (unsigned int i = 0; i < n_electrons; i++)
    {
      // We choose the electron starting position at random from a flat distribution along the path length
      // the parameter t is the fraction of the distance along the path betwen entry and exit points, it has values between 0 and 1
      double f = gsl_ran_flat(RandomGenerator, 0.0, 1.0);

      x_start = hiter->second->get_x(0) + f * (hiter->second->get_x(1) - hiter->second->get_x(0));
      y_start = hiter->second->get_y(0) + f * (hiter->second->get_y(1) - hiter->second->get_y(0));
      double z_start = hiter->second->get_z(0) + f * (hiter->second->get_z(1) - hiter->second->get_z(0));
      double t_start = hiter->second->get_t(0) + f * (hiter->second->get_t(1) - hiter->second->get_t(0));
      if(do_Centralmem)
	{
	  Double_t ax[7668],ay[7668];
	  CM->GetPoint(i,ax[i],ay[i]);
	  //cout<<i<<"th element of X array: "<<ax[i]<<endl;
	  // cout<<i<<"th element of Y array: "<<ay[i]<<endl;
	  x_start = ax[i]*.1;//[i];                                                   
	  // cout << "x_start is " << x_start << endl;
	  y_start = ay[i]*.1;//[i];                                                                                                                        
	  //cout << "y_start is " << y_start << endl;
	  z_start = 0;
	}
      double r_sigma = diffusion_trans * sqrt(tpc_length / 2. - fabs(z_start));
      double rantrans = gsl_ran_gaussian(RandomGenerator, r_sigma);
      rantrans += gsl_ran_gaussian(RandomGenerator, added_smear_sigma_trans);

      double t_path = (tpc_length / 2. - fabs(z_start)) / drift_velocity;
      double t_sigma = diffusion_long * sqrt(tpc_length / 2. - fabs(z_start)) / drift_velocity;
      double rantime = gsl_ran_gaussian(RandomGenerator, t_sigma);
      rantime += gsl_ran_gaussian(RandomGenerator, added_smear_sigma_long) / drift_velocity;
      double t_final = t_start + t_path + rantime;

      double z_final;
      if (z_start < 0)
        z_final = -tpc_length / 2. + t_final * drift_velocity;
      else
        z_final = tpc_length / 2. - t_final * drift_velocity;

      if ((t_final < min_time || t_final > max_time) && !do_Centralmem)
       {
      cout << "skip this, t_final = " << t_final << " is out of range " << min_time <<  " to " << max_time << endl;
       continue;
       }

      double radstart = sqrt(x_start * x_start + y_start * y_start);
      double phistart = atan2(y_start,x_start);
      double ranphi = gsl_ran_flat(RandomGenerator, -M_PI, M_PI);
      double rad_final,phi_final,x_final,y_final;
      z_startmap->Fill(z_start,radstart); // map of starting location in Z vs. R
      x_final = x_start + rantrans * cos(ranphi); // Initialize these to be only smeared first, will be overwritten if doing diff or int SC distortion                     
      y_final = y_start + rantrans * sin(ranphi); 
      rad_final = sqrt(pow(x_final,2)+pow(y_final,2));
      phi_final = atan2(y_final,x_final);
      if(z_start < 0)
	{
	  z_start = fabs(z_start);// Eventually will require two maps for pos and neg Z 
	}
        
      double r_final_diff = 0;
      double phi_final_diff = 0;
      double r_final_int = 0;
      double phi_final_int = 0;
      bool do_diff_SC_Distortion = true;  
      bool do_Int_SC_Distortion = true; // Decide whether or not to do space charge distortion     
      if (do_Int_SC_Distortion)
	{
	  intwatch->Start(false);
	  rad_final = radstart+TimeInthDR->Interpolate(phistart+M_PI,radstart,z_start);
	  phi_final = phistart+M_PI+(TimeInthDP->Interpolate(phistart+M_PI,radstart,z_start)/radstart);
	  x_final = rad_final*cos(phi_final-M_PI)+rantrans*cos(ranphi);
	  y_final = rad_final*sin(phi_final-M_PI)+rantrans*sin(ranphi);
	  deltarint->Fill(radstart,sqrt(pow(x_final,2)+pow(y_final,2))-sqrt(pow(x_start,2)+pow(y_start,2))); //total delta r
	  deltaphiint->Fill(phistart+M_PI,phistart-atan2(y_final,x_final)); // total delta phi	
	  r_final_int = rad_final;
	  phi_final_int = phi_final;
	  intwatch->Stop();
	}
      if(do_diff_SC_Distortion)
	{
	  diffwatch->Start(false);
	  DistortionIntegral(radstart,phistart,z_start,&rad_final,&phi_final);
	  x_final = rad_final*cos(phi_final-M_PI)+rantrans*cos(ranphi);
	  y_final = rad_final*sin(phi_final-M_PI)+rantrans*sin(ranphi);
	  deltardiff->Fill(radstart,sqrt(pow(x_final,2)+pow(y_final,2))-sqrt(pow(x_start,2)+pow(y_start,2))); //total delta r
	  deltaphidiff->Fill(phistart+M_PI,phistart-atan2(y_final,x_final)); // total delta phi
	  r_final_diff = rad_final;
	  phi_final_diff = phi_final;
	  diffwatch->Stop();	
	}
      //cout << "Values of branch variables before filling are: event_num,x_start,y_start,x_final,y_final " << event_num << " , " << x_start << " , " << y_start << " , " << x_final << " , " << y_final <<  endl;
      CMTimeDists->Fill();
            
      deltardifference->Fill(rad_final,r_final_diff - r_final_int);
      deltardifferencepercent->Fill(rad_final,100*((r_final_diff-radstart) - (r_final_int-radstart))/(r_final_diff-radstart));
      deltaphidifference->Fill(phi_final,phi_final_diff - phi_final_int);
      deltaphidifferencepercent->Fill(phi_final,100*((phi_final_diff-phistart) - (phi_final_int-phistart))/(phi_final_diff-phistart));
      
      // Diagnostic plots, written into ElectronDriftQA.root, some won't make sense if do_SC_Distortion set to false
      hitmapstart->Fill(x_start,y_start); // G4Hit starting positions
      hitmapend->Fill(x_final,y_final);//INcludes diffusion and distortion          
      deltar->Fill(radstart,sqrt(pow(x_final,2)+pow(y_final,2))-sqrt(pow(x_start,2)+pow(y_start,2))); //total delta r
      deltaphi->Fill(phistart+M_PI,phistart-atan2(y_final,x_final)); // total delta phi
      deltaphinodist->Fill(phistart+M_PI,rantrans/rad_final); // delta phi no distortion, just diffusion+smear
      deltarnodist->Fill(radstart,rantrans); // delta r no distortion, just diffusion+smear
      deltaphinodiff->Fill(phistart+M_PI,(phistart+M_PI)-phi_final); //delta phi no diffusion, just distortion
      deltarnodiff->Fill(radstart,(rad_final-radstart));//delta r no diffusion, just distortion      
      deltaz->Fill(z_start,hDZint->Interpolate(phistart+M_PI,radstart,z_start)); // map of distortion in Z (time)
   
      // remove electrons outside of our acceptance. Careful though, electrons from just inside 30 cm can contribute in the 1st active layer readout, so leave a little margin
      if (rad_final < min_active_radius - 2.0 || rad_final > max_active_radius + 1.0)
      {
        continue;
      }

      if (Verbosity() > 1000)
      {
        cout << "electron " << i << " g4hitid " << hiter->first << " f " << f << endl;
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
    ihit++;
   
    deltatime->Fill(diffwatch->RealTime() - intwatch->RealTime());
    // transfer the hits from temp_hitsetcontainer to hitsetcontainer on the node tree
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
	if(Verbosity()>100)   
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
	      }
	    
	    // Either way, add the energy to it
	    node_hit->addEnergy(temp_tpchit->getEnergy());
	    
	    // Add the hit-g4hit association	    
	    hittruthassoc->findOrAddAssoc(node_hitsetkey, temp_hitkey, hiter->first);
	    
	  }  // end loop over temp hits
	if(Verbosity() > 100 && layer == print_layer)
	  cout << "  ihit " << ihit << " collected energy = " << eg4hit << endl;

      } // end loop over temp hitsets

    // erase all entries in the temp hitsetcontainer
    temp_hitsetcontainer->Reset();
    
  } // end loop over g4hits
  
  unsigned int print_layer = 47;  

  if(Verbosity() > 2)
    {
      cout << "From PHG4TpcElectronDrift: hitsetcontainer printout at end:" << endl;
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

	      cout << "      hitkey " << hitkey << " pad " << TpcDefs::getPad(hitkey) << " z bin " << TpcDefs::getTBin(hitkey) 
		   << "  energy " << tpchit->getEnergy() << " adc " << tpchit->getAdc() << endl;
	    }
	}
    }

  if(Verbosity() > 2)
    {
      cout << "From PHG4TpcElectronDrift: hittruthassoc dump:" << endl;
      hittruthassoc->identify();
    }
  if(Verbosity() == 12)
    {
      Int_t nentries = (Int_t)TimeTree->GetEntries();  
      cout << "nentries is " << nentries << endl;
      cout << "deltatime is " << (diffwatch->RealTime() - intwatch->RealTime()) << endl;
      cout << "integrated drifting time is " <<  (intwatch->RealTime()) << endl;
      cout << "differential drifting time is " << (diffwatch->RealTime()) << endl;
      cout << "non-drifting time in PHG4ElectronDrift is " << (totwatch->RealTime()-(diffwatch->RealTime()+intwatch->RealTime())) << endl;
      cout << "total time in PHG4ElectronDrift is " << (totwatch->RealTime()) << endl;
    }
  event_num += 1;
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TpcElectronDrift::MapToPadPlane(const double x_gem, const double y_gem, const double t_gem, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit)
{
  padplane->MapToPadPlane(temp_hitsetcontainer, hittruthassoc, x_gem, y_gem, t_gem, hiter, ntpad, nthit);

  return;
}

int PHG4TpcElectronDrift::End(PHCompositeNode *topNode)
{

  if (true) // choose verbosity later
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
    
    CM_outf->cd();
    CMTimeDists->Write();
    CM_outf->Close();

    EDrift_outf->cd();
    // TimehDR->Write();
    // TimehDP->Write();
    deltar->Write();
    deltardiff->Write();
    deltarint->Write();
    deltatime->Write();
    deltardifference->Write();
    deltardifferencepercent->Write();
    deltaphi->Write();
    deltaphidiff->Write();
    deltaphidifference->Write();
    deltaphidifferencepercent->Write();
    deltaphiint->Write();
    deltaz->Write();
    deltarnodist->Write();
    deltaphinodist->Write();
    deltarnodiff->Write();
    deltaphinodiff->Write();
    hitmapstart->Write();
    hitmapend->Write();
    z_startmap->Write();
    EDrift_outf->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TpcElectronDrift::set_seed(const unsigned int iseed)
{
  seed = iseed;
  gsl_rng_set(RandomGenerator, seed);
}

void PHG4TpcElectronDrift::SetDefaultParameters()
{
  // Data on gasses @20 C and 760 Torr from the following source:
  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
  // diffusion and drift velocity for 400kV for NeCF4 90/10 from calculations:
  // https://www.phenix.bnl.gov/WWW/p/draft/prakhar/tpc/HTML_Gas_Linear/Ne_CF4_90_10.html
  double Ne_dEdx = 1.56;   // keV/cm
  double CF4_dEdx = 7.00;  // keV/cm
  // double Ne_NPrimary = 12;    // Number/cm
  // double CF4_NPrimary = 51;   // Number/cm
  double Ne_NTotal = 43;    // Number/cm
  double CF4_NTotal = 100;  // Number/cm
  double Tpc_NTot = 0.9 * Ne_NTotal + 0.1 * CF4_NTotal;
  double Tpc_dEdx = 0.90 * Ne_dEdx + 0.10 * CF4_dEdx;
  double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  set_default_double_param("diffusion_long", 0.015);   // cm/SQRT(cm)
  set_default_double_param("diffusion_trans", 0.006);  // cm/SQRT(cm)
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
  padplane = inpadplane;
  padplane->Detector(Detector());
  padplane->UpdateInternalParameters();
  cout << "padplane registered and parameters updated" << endl;

  return;
}
