/*!
        \file PHEventDisplay.cxx
        \author Sookhyun Lee
        \brief event display interface,
	       set parameters/switches, call detector modules, control display.
        \version $Revision: 1.2 $
        \date    $Date: 07/26/2016
*/

// STL and BOOST includes
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

// PHOOL and Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTimer.h>
#include <phool/PHCompositeNode.h>
//#include <TMutNode.h>
#include <phool/PHTimeServer.h>
#include <phgeom/PHGeomUtility.h>
#include <phgeom/PHGeomTGeo.h>


// EVE framework includes
#include "TEveManager.h"
#include "TEveBrowser.h"
#include "TEveWindow.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGLAutoRotator.h"
#include "TGLViewer.h"
#include "TFile.h"


// sPHENIX Event Display
#include <PHEveDisplay.h>

// Geant4 truth info display
#include <mG4EveDisplay.h>
#include <mSvtxEveDisplay.h>
#include <mCaloEveDisplay.h>
#include <mJetEveDisplay.h>

#include <TThread.h>
#include <TStyle.h>

#include <PHEventDisplay.h>

using boost::bind;

PHEventDisplay::PHEventDisplay(int w = 1920,
				 int h = 1080,
				 bool _use_fieldmap = true,
				 bool _use_geofile = true,
				 const std::string& _mapname = "sPHEBIX.2d.root",
				 const std::string& _geoname = "geo.root") :
  SubsysReco("PHEventDisplay"),
  _pending_update(false),
  _modules(),
  _mutex(PTHREAD_MUTEX_INITIALIZER),
  _update_thread(),
  _PHEveDisplay(new PHEveDisplay(w, h, _use_fieldmap, _use_geofile, _mapname, _geoname, verbosity)),
  _rot(NULL),
  _status_thread(),
  jet_pt_threshold(5.),
  jet_e_scale(30.),
  calo_e_threshold(0.2),
  is_svtx_on(true),
  is_cemc_on(true),
  is_hcalin_on(true),
  is_hcalout_on(true),
  is_jet_on(true),
  is_truth_on(true),
  use_fieldmap(_use_fieldmap),
  use_geofile(_use_geofile),
  width(w),
  height(h),
  mapname(_mapname),
  geoname(_geoname),
  nevent(0),
  verbosity(0)
{

  TEveManager::Create(kTRUE,"V");
  gEve->GetBrowser()->HideBottomTab();
  TEveWindow::SetMainFrameDefWidth(width);
  TEveWindow::SetMainFrameDefHeight(height);
  TFile::SetCacheFileDir(".");

  gEve->FullRedraw3D(kTRUE);  
  TGLViewer*  v = gEve->GetDefaultGLViewer();
  //  v->ColorSet().Background().SetColor(kMagenta+4);
        v->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
  v->RefreshPadEditor(v);
  v->CurrentCamera().RotateRad(0.0, 1.5707);// theta, phi
  v->DoDraw();

  _PHEveDisplay->set_eve_manager(gEve);

  if (verbosity) std::cout << "PHEventDisplay instantiated."<<std::endl;  
}

PHEventDisplay::~PHEventDisplay()
{
}

int PHEventDisplay::Init(PHCompositeNode *topNode)
{
  SubsysReco::Init(topNode);

  _PHEveDisplay->set_jet_pt_threshold(jet_pt_threshold);
  _PHEveDisplay->set_jet_e_scale(jet_e_scale);
  _PHEveDisplay->set_calo_e_threshold(calo_e_threshold);

  std::cout.precision(5);
  std::cout<< "jet pt threshold: "<< jet_pt_threshold<<std::endl;
  std::cout<< "jet e scale: " << jet_e_scale<< std::endl;
  std::cout<< "calo e threshold: " << calo_e_threshold << std::endl;
  std::cout<< "svtx :" << is_svtx_on << ", cemc : "<< is_cemc_on << ", hcalin : "<< is_hcalin_on << ", hcalout : "<< is_hcalout_on << ", jet : "<<is_jet_on <<", truth : "<< is_truth_on <<", verbosity : " << verbosity <<std::endl;

  return 0;
}
 
int PHEventDisplay::InitRun(PHCompositeNode *topNode)
{

  if (verbosity) std::cout<<"PHEventDisplay - initialize run. "<<std::endl;

  try
    {

    if (verbosity) std::cout<<"PHEventDisplay - load_geometry() begins. "<<std::endl;
    _PHEveDisplay->load_geometry(topNode, gEve);
    if (verbosity) std::cout<<"PHEventDisplay - add_elements() begins."<<std::endl;
    _PHEveDisplay->add_elements(gEve);
    if (verbosity) std::cout<<"PHEventDisplay - config_bfield() begins."<<std::endl;
    _PHEveDisplay->config_bfields();

    if(is_truth_on){
 	register_module<mG4EveDisplay>();
	if(verbosity)
	std::cout<<"PHEventDisplay - mG4EveDisplay module registered."<<std::endl;}
    if(is_svtx_on){
	register_module<mSvtxEveDisplay>();
	if (verbosity) 
	std::cout<<"PHEventDisplay - mSvtxEveDisplay module registered."<<std::endl;}
    if(is_cemc_on || is_hcalin_on || is_hcalout_on){
	register_module<mCaloEveDisplay>();
	if(verbosity)
	std::cout<<"PHEventDisplay - mCaloEveDisplay module registered."<<std::endl;}
    if(is_jet_on){
	register_module<mJetEveDisplay>();
	if(verbosity)
	std::cout<<"PHEventDisplay - mJetEveDisplay module registered."<<std::endl;}
	std::for_each(_modules.begin(),_modules.end(),
                bind(&mPHEveModuleBase::init,_1, topNode));
    if(verbosity) std::cout<<"PHEventDisplay - all modules registered. "<<std::endl;

    }
  catch(std::exception &e)
    {
        std::cout << "Exception caught while initializing sPHENIX event display: " << e.what() << std::endl;
    }


  return 0;
}

int PHEventDisplay::process_event(PHCompositeNode *topNode)
{
  if (verbosity) std::cout<<"PHEventDisplay - setting up to process event." <<std::endl;

  std::for_each(_modules.begin(),
		_modules.end(),
		bind(&mPHEveModuleBase::event,
		_1,
		topNode));
  return 0;
}

int PHEventDisplay::End(PHCompositeNode *topNode)
{
  return 0;
}

void PHEventDisplay::update_scene()
{
  if (verbosity) std::cout << "PHEventDisplay - update_scene() nevent = " <<nevent<<std::endl;
  TThread::Lock(); // Keep the autorotator from segfaulting
  std::for_each(_modules.begin(),_modules.end(),
		bind(&mPHEveModuleBase::clear,_1));
  draw_default();
  TThread::UnLock();
}

void PHEventDisplay::run_evt_in_thread()
{
  if(verbosity) std::cout << "PHEventDisplay - run_evt_in_thread() nevent = "<<nevent<<std::endl;
  if( !pthread_mutex_trylock(&_mutex) ) 
  {
    if( _pending_update ){
	update_scene();
    }
    _pending_update = false;
    _update_thread = boost::make_shared<boost::thread>(bind(&PHEventDisplay::reco_thread, this));
  }
}


void PHEventDisplay::reco_thread()
{
  if (verbosity) std::cout<< "PHEventDisplay - reco_thread() nevent = "<<nevent<<std::endl;
  nevent++;
  Fun4AllServer* se = Fun4AllServer::instance();
  se->run(1);
  if (verbosity>1) std::cout <<"reco_thread() run(1) complete."<<std::endl;
  _pending_update = true;
  if (verbosity>1) std::cout<< "reco_thread() pendinding update "<<std::endl;
  pthread_mutex_unlock(&_mutex);
  if (verbosity>1) std::cout <<"reco_thread() mutex unlock()"<<std::endl;

}

void PHEventDisplay::draw_default()
{
  if (verbosity) std::cout<<"PHEventDisplay - draw_default() begins."<<std::endl;
  if (verbosity>1) std::cout<<"draw_default() FullRedraw3D." <<std::endl;
  gEve->FullRedraw3D(kTRUE);
  if (verbosity>1) std::cout<<"draw_default() TGLViewer." <<std::endl;
  TGLViewer*  v = gEve->GetDefaultGLViewer();
  //  v->ColorSet().Background().SetColor(kMagenta+4);
  v->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
  v->RefreshPadEditor(v);
  v->CurrentCamera().RotateRad(0.0, 1.5707);// theta, phi
  v->DoDraw();

}


void PHEventDisplay::start_rotation()
{
  _rot = gEve->GetDefaultGLViewer()->GetAutoRotator();
  _rot->SetADolly(0.0);
  _rot->SetATheta(0.0);
  _rot->SetWTheta(0.075);
  _rot->SetWPhi(0.1);
  _rot->Start();
}

void
PHEventDisplay::go_fullscreen()
{
  _PHEveDisplay->go_fullscreen(gEve);
}

