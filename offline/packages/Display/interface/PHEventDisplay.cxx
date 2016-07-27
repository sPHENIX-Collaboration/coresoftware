/*!
        \file PHEventDisplay.cxx
        \author Sookhyun Lee
        \brief event display interface,
	       parameters/switches set, detector modules called, control display.
        \version $Revision: 1.1 $
        \date    $Date: 07/26/2016
*/

// STL and BOOST includes
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

// PHOOL and Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHIODataNode.h>
#include <phool/PHTimer.h>
#include <phool/PHCompositeNode.h>
//#include <TMutNode.h>
#include <phool/PHTimeServer.h>

// sPHENIX Event Display
#include <PHEveDisplay.h>

// Geant4 truth info display
#include <mG4EveDisplay.h>
#include <mSvtxEveDisplay.h>
#include <mCaloEveDisplay.h>
#include <mJetEveDisplay.h>

// EVE framework includes
#include "TEveManager.h"
#include <TGLViewer.h>
#include <TGLAutoRotator.h>
#include <TThread.h>

#include <PHEventDisplay.h>

using boost::bind;

PHEventDisplay::PHEventDisplay(int w = 1920,
				 int h = 1080,
				 bool use_fieldmap = true,
				 const std::string& mapname = "sPHEBIX.2d.root",
				 const std::string& geoname = "geo.root",
				 float jet_pt_threshold = 5.,
				 float jet_e_scale = 20.,
				 float calo_e_threshold = 0.2,
				 bool is_svtx_on = true,
				 bool is_calo_on = true,
				 bool is_jet_on = true) :
  SubsysReco("PHEventDisplay"),
  _pending_update(false),
  _modules(),
  _mutex(PTHREAD_MUTEX_INITIALIZER),
  _update_thread(),
  _PHEveDisplay(new PHEveDisplay(w, h, use_fieldmap, mapname, geoname)),
  _eve_manager(NULL),
  _rot(NULL),
  _status_thread(),
  _is_svtx_on(is_svtx_on),
  _is_calo_on(is_calo_on),
  _is_jet_on(is_jet_on)
{
  _eve_manager = _PHEveDisplay->get_eve_instance();
  _PHEveDisplay->set_jet_pt_threshold(jet_pt_threshold);
  _PHEveDisplay->set_jet_e_scale(jet_e_scale);
  _PHEveDisplay->set_calo_e_threshold(calo_e_threshold);
  std::cout << "PHEventDisplay initialized.. "<<std::endl;  
}

PHEventDisplay::~PHEventDisplay()
{
}

int PHEventDisplay::Init(PHCompositeNode *topNode)
{
  SubsysReco::Init(topNode);

  register_module<mG4EveDisplay>();
  std::cout<<"mG4EveDisplay module registered.."<<std::endl;
  if(_is_svtx_on){ 
  register_module<mSvtxEveDisplay>();
  std::cout<<"mSvtxEveDisplay module registered.."<<std::endl;}
  if(_is_calo_on){
  register_module<mCaloEveDisplay>();
  std::cout<<"mCaloEveDisplay module registered.."<<std::endl;}
  if(_is_jet_on){
  register_module<mJetEveDisplay>();
  std::cout<<"mJetEveDisplay module registered.."<<std::endl;}

  std::for_each(_modules.begin(),_modules.end(),
		bind(&mPHEveModuleBase::init,_1, topNode));
  std::cout<<"all modules registered.. "<<std::endl;
  return 0;
}
 
int PHEventDisplay::InitRun(PHCompositeNode *topNode)
{
  std::cout<<"PHEventDisplay - initialize run.. "<<std::endl;
  return 0;
}

int PHEventDisplay::process_event(PHCompositeNode *topNode)
{
  std::cout<<"PHEventDisplay - setting up to process event.." <<std::endl;
  std::for_each(_modules.begin(),
		_modules.end(),
		bind(&mPHEveModuleBase::event,
		_1,
		topNode));
  return 0;
}


int PHEventDisplay::ResetEvent(PHCompositeNode *topNode)
{
  return 0;
}

int PHEventDisplay::End(PHCompositeNode *topNode)
{
  return 0;
}

TEveManager* PHEventDisplay::get_eve_instance()
{
  return _PHEveDisplay->get_eve_instance();
}

void PHEventDisplay::update_scene()
{
  TThread::Lock(); // Keep the autorotator from segfaulting
  std::for_each(_modules.begin(),_modules.end(),
		bind(&mPHEveModuleBase::draw_event,_1));
  _eve_manager->Redraw3D(kFALSE);
  TThread::UnLock();
}

void PHEventDisplay::run_evt_in_thread()
{
  if( !pthread_mutex_trylock(&_mutex) ) 
  {
    if( _pending_update ) update_scene();
    _pending_update = false;
    _update_thread = boost::make_shared<boost::thread>(bind(&PHEventDisplay::reco_thread, this));
  }
}

void PHEventDisplay::reco_thread()
{
  Fun4AllServer* se = Fun4AllServer::instance();
  se->run(1);
  _pending_update = true;
  pthread_mutex_unlock(&_mutex);
}

void PHEventDisplay::start_rotation()
{
  _rot = _eve_manager->GetDefaultGLViewer()->GetAutoRotator();
  _rot->SetADolly(0.0);
  _rot->SetATheta(0.0);
  _rot->SetWTheta(0.075);
  _rot->SetWPhi(0.1);
  _rot->Start();
}

void
PHEventDisplay::go_fullscreen()
{
  _PHEveDisplay->go_fullscreen();
}

