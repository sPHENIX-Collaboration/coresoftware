/*
        \file PHEveDisplay.cxx
        \author Sookhyun Lee
        \brief main display module, 
	       load geometry, configure b-field, draw default.
        \version $Revision: 1.1 $
        \date    $Date: 07/26/2016
*/

// STL and BOOST includes
#include<iostream>
#include<string>
#include<stdexcept>
#include<boost/shared_ptr.hpp>

// EVE class includes
#include "TEveManager.h"
#include "TEveVSDStructs.h"
#include "TEveGeoNode.h"
#include "TEveVector.h"
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "TEvePointSet.h"
#include "TEveWindowManager.h"
#include "TEveWindow.h"
#include "TEveViewer.h"
#include "TEveBrowser.h"

#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"

#include "TGLViewer.h"
#include "TGPack.h"

#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TMath.h"

#include "PHBFieldMap.hh"

#include "PHEveDisplay.h"

PHEveDisplay::PHEveDisplay(int w,
			   int h,
			   bool use_fieldmap,
			   const std::string& mapname,
			   const std::string& geoname) :
  _top_list(NULL),
  _svtx_list(NULL),
  _calo_list(NULL),
  _jet_list(NULL),
  _true_list(NULL),
  cnt_prop(NULL),
  mapped_field(NULL),
  _width(w),
  _height(h),
  _use_fieldmap(use_fieldmap),
  _jet_pt_threshold(5.0),
  _calo_e_threshold(0.2),
  geo_filename(geoname),
  map_filename(mapname)
{
  TEveManager::Create(kTRUE,"V");
  gEve->GetBrowser()->HideBottomTab();
  TEveWindow::SetMainFrameDefWidth(_width);
  TEveWindow::SetMainFrameDefHeight(_height);
  TFile::SetCacheFileDir(".");

  try
    {
	std::cout<<"begin load_geometry()"<<std::endl;
        load_geometry();
	std::cout<<"begin draw_default()"<<std::endl;
        draw_default();
	std::cout<<"begin config_bfield()"<<std::endl;
        config_bfields();
    }
  catch(std::exception &e)
    {
        std::cout << "Exception caught while initializing sPHENIX event display: " << e.what() << std::endl;
    }

}  

PHEveDisplay::~PHEveDisplay()
{
  try
    {
	TEveManager::Terminate();
    }
  catch(std::exception &e)
    {
	std::cout << "Exception caught during deconstruction: " << e.what() << std::endl;
    }

}

void 
PHEveDisplay::load_geometry()
{
  TFile* geom = TFile::Open(geo_filename.c_str());
  if (!geom)
    throw std::runtime_error("Could not open sphenix_geo.root geometry file, aborting.");

  gEve->GetGeometry(geo_filename.c_str());
  //gGeoManager->DefaultColors();
  gStyle->SetPalette(1);
  TGeoVolume* top = gGeoManager->GetTopVolume();
  const int nd = top->GetNdaughters();
  TGeoNode* node[nd];
  TEveGeoTopNode* tnode[nd];
    int det_config = 0;
    //bool is_supp_struc = false;
    if (strcmp(geo_filename.c_str(),"sphenix_mie_geo.root")==0) det_config = 1;
    else if(strcmp(geo_filename.c_str(),"sphenix_maps+tpc_geo.root")==0) det_config = 2;
      for(int i=0 ; i< nd; i++)
      {
        if(det_config==1 && i==18) continue;
        if(det_config==2 && (i==73 || (i>10 && i<71))) continue;
        node[i]=top->GetNode(i);
          node[i]->GetVolume()->SetTransparency(70);// 0: opaque, 100: transparent
	  if(det_config==2)
 	  {
	  //is_supp_struc = (i==4 || i==6 || (i>=8 && i<=10) || i==71 || i==74 || (i>=76&&i<=79));
          //if(is_supp_struc) // set color to gray if support structure
          //node[i]->GetVolume()->SetLineColor(kGray+2);
	    if(i==75 || i==83)
	    { // make hcal transparent
	    TGeoVolume* hcalvol = node[i]->GetVolume();
	    const int nhcal = hcalvol->GetNdaughters();
	    TGeoNode* node_hcal[nhcal];
	      for(int j=0 ; j<nhcal; j++)
	      {
	      node_hcal[j] = hcalvol->GetNode(j);
	      node_hcal[j]->GetVolume()->SetTransparency(90);
	      }
	    }
	  }
        tnode[i] = new TEveGeoTopNode(gGeoManager, node[i]);  
        gEve->AddGlobalElement(tnode[i]);
      }
  geom->Close();
  delete geom;

  _top_list = new TEveElementList("TOP");
  _svtx_list = new TEveElementList("SVTX");
  _calo_list = new TEveElementList("CALO");
  _jet_list = new TEveElementList("JET");
  _true_list = new TEveElementList("TRUE"); 
  

  gEve->AddElement(_top_list);
  gEve->AddElement(_svtx_list,_top_list);
  gEve->AddElement(_calo_list,_top_list);
  gEve->AddElement(_jet_list,_top_list);
  gEve->AddElement(_true_list,_top_list);

} 


void 
PHEveDisplay::draw_default()
{

  gEve->FullRedraw3D(kTRUE);

  TGLViewer *v = gEve->GetDefaultGLViewer();
//  v->ColorSet().Background().SetColor(kMagenta+4);
  v->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
  v->RefreshPadEditor(v);
  v->CurrentCamera().RotateRad(0.0, 1.5707);// theta, phi
  v->DoDraw();

}

void
PHEveDisplay::config_bfields()
{
  if( _use_fieldmap )
    {
	std::cout << "PHEveDisplay::config_bfields:"
		  << "Using mapped magnetic fields for track propagation"
		  << std::endl;
	cnt_prop = new TEveTrackPropagator("cnt_prop", 
					     "Central Field Propagator",
					     new MappedField(map_filename));
	cnt_prop->SetMaxStep(2);
    }
  else
    {
	std::cout << "PHEveDisplay::config_bfields:"
		  << "Analytic Form of Approximate Magnetic Field Not Availiable Yet."
		  << std::endl;
    }

  cnt_prop->SetStepper(TEveTrackPropagator::kHelix);
  cnt_prop->SetMaxR(350);
  cnt_prop->SetMaxZ(450);
}

void
PHEveDisplay::go_fullscreen()
{
  TEveViewer* cur_win = gEve->GetDefaultViewer();
  TEveCompositeFrame* fEveFrame = cur_win->GetEveFrame();
  TEveWindow* return_cont = fEveFrame->GetEveParentAsWindow();

  if (return_cont && ! return_cont->CanMakeNewSlots())
  return_cont = 0;

  TEveCompositeFrameInPack* packframe = dynamic_cast<TEveCompositeFrameInPack*>(fEveFrame);
  if (packframe) {
    TGPack* pack = (TGPack*)(packframe->GetParent());
    pack->HideFrame(fEveFrame);
  }

  TGMainFrame* mf = new TGMainFrame(gClient->GetRoot(),_width,_height);
  gVirtualX->SetMWMHints(mf->GetId(),0,0,0);
  mf->SetCleanup(kLocalCleanup);

  TEveCompositeFrameInMainFrame *slot = new TEveCompositeFrameInMainFrame(mf, 0, mf);
  TEveWindowSlot* ew_slot = TEveWindow::CreateDefaultWindowSlot();
  ew_slot->PopulateEmptyFrame(slot);

  mf->AddFrame(slot, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY));
  slot->MapWindow();

  mf->Layout();
  mf->MapWindow();

  TEveWindow::SwapWindows(ew_slot, cur_win);

  ((TEveCompositeFrameInMainFrame*) fEveFrame)->
    SetOriginalSlotAndContainer(ew_slot, return_cont);

  gEve->GetWindowManager()->HideAllEveDecorations();
  gEve->GetWindowManager()->WindowUndocked(cur_win);

  int offset = -8;

  gVirtualX->MoveResizeWindow(mf->GetId(), 0, offset, _width, _height);
}

void 
PHEveDisplay::update()
{
  gEve->Redraw3D(kFALSE, kTRUE);
}

TEveManager* PHEveDisplay::get_eve_instance() const
{
  return gEve;
}

PHEveDisplay::MappedField::MappedField(const std::string& fname) :
  TEveMagField(),
  _fieldmap(new PHBFieldMap(fname))
{
}

TEveVectorD
PHEveDisplay::MappedField::GetFieldD(Double_t x, Double_t y, Double_t z) const
{
  double loc[3];
  loc[0]=x; loc[1]=y; loc[2]=z; 
  double bvec[3];
   _fieldmap->get_bfield(&loc[0],&bvec[0]);
  TEveVectorD vec(bvec[0],
		  bvec[1],
		  bvec[2]); // unit is Gauss  (1Tesla = 10000 Gauss)

  return vec;
}

