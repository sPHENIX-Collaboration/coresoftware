
/*!
	\file mG4EveDisplay.cxx
	\author Sookhyun Lee
	\brief true tracks and true jets from truthinfo container
	\version $Revision: 1.1 $
 	\date    $Date: 07/26/2016
*/

// STL and BOOST includes
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <boost/bind.hpp>

// PHENIX includes
#include <phool/PHCompositeNode.h>
//#include <PHPoint.h>
#include <phool/getClass.h>

// ROOT and EVE includes
#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TEveTrack.h>
#include <TEvePointSet.h>
#include <TEveElement.h>
#include <TEveJetCone.h>
#include <TEveStraightLineSet.h>
#include <TEveVector.h>
#include <TEveCaloData.h>
#include <TEveCalo.h>
#include <TH2F.h>

// truth info includes
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPointv1.h> // <-> reco  vertex
#include <g4main/PHG4Showerv1.h>   // <-> reco  cluster
#include <g4jets/JetMap.h>         // for jets
#include <g4jets/Jet.h>

#include <PHEveDisplay.h>
#include <mG4EveDisplay.h>


using boost::bind;
using namespace TMath;

mG4EveDisplay::mG4EveDisplay(boost::shared_ptr<PHEveDisplay> dispin) :
  mPHEveModuleBase(),
  _evedisp(dispin),
  _truth(NULL),   // truth info container
  _jetmap(NULL),
  _prop(NULL),
  _true_tracks(NULL),
  _true_jets(NULL),    
  radius(0.3),
  length(300)
{
  _evemanager = _evedisp->get_eve_instance();
  _prop = _evedisp->get_cnt_prop();
  _true_tracks = new TEveTrackList("True tracks");
  _true_jets = new TEveElementList("True jets");
  _evemanager->AddElement(_true_tracks,_evedisp->get_true_list());
  _evemanager->AddElement(_true_jets,_evedisp->get_true_list());
}

mG4EveDisplay::~mG4EveDisplay()
{
}

void
mG4EveDisplay::init(PHCompositeNode* topNode)
{  
}

void 
mG4EveDisplay::init_run(PHCompositeNode* topNode)
{
}

bool
mG4EveDisplay::event(PHCompositeNode* topNode)
{
  std::cout<<"mG4EveDisplay - event.."<<std::endl;
  try
    {
      create_nodes(topNode);
      draw_tracks();      
      draw_jets();
    }
  catch(std::exception& e)
    {
      static bool first(true);
      if( first )
	std::cout << "mG4EveDisplay::event Exception: " << e.what() << std::endl;
      first = false;
    }

  return true;

}

void
mG4EveDisplay::end(PHCompositeNode* topNode)
{
}

void 
mG4EveDisplay::create_nodes(PHCompositeNode* topNode)
{
  _truth = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!_truth) std::cout<<"TruthInfoContainer node not found!! "<<std::endl;
  _jetmap = findNode::getClass<JetMap>(topNode,"AntiKt_Truth_r03"); 
  if (!_jetmap) std::cout<<"JetMap node not found!! "<<std::endl;
  std::cout<<"mG4EveDisplay - nodes created.."<<std::endl;
}

void
mG4EveDisplay::draw_tracks()
{
  std::cout<<"mG4EveDisplay - draw_tracks() begins.. "<<std::endl;

    // true primary particle tracks
    PHG4TruthInfoContainer::ConstRange range = _truth->GetPrimaryParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
	 iter != range.second;
	 ++iter) {
	  PHG4Particle *particle = iter->second;
	  TEveRecTrackT<double>* trk_reco = new TEveRecTrackT<double>();
	  trk_reco->fV.Set(_truth->GetVtx(particle->get_vtx_id())->get_x(),
			   _truth->GetVtx(particle->get_vtx_id())->get_y(),
			   _truth->GetVtx(particle->get_vtx_id())->get_z());
	  trk_reco->fP.Set(particle->get_px(),
			   particle->get_py(),
			   particle->get_pz());
	  int pid = particle->get_pid();
	  int tid = particle->get_track_id();
	  std::cout<<"mG4EveDisplay - particle id : "<<pid<<", track id : "<<tid<<std::endl;
	  if(fabs(pid)==211 || fabs(pid)==321 || fabs(pid)==11 || fabs(pid)==13 || fabs(pid)==15 || fabs(pid)==17){
	  if(pid>0)
	  trk_reco->fSign = 1;
	  else
	  trk_reco->fSign = -1;
//	  _prop->SetMagField(0.);
	  TEveTrack* trk = new TEveTrack(trk_reco, _prop);
	  std::cout<<"mG4EveDisplay - EveTrack instantiated.."<<std::endl;
	  trk->SetSmooth(kTRUE);
	  // Color coding and 
	  if(pid>0)
          trk->SetLineColor(kOrange);
	  else
	  trk->SetLineColor(kBlue);
          trk->SetLineWidth(2.5);
	  std::cout<<"mG4EveDisplay - set cosmetices of track"<<std::endl;	
	  trk->MakeTrack();            
	  std::cout<<"mG4EveDisplay - track made.. "<<std::endl;
	  _true_tracks->AddElement(trk);
	  }
      }

}

void 
mG4EveDisplay::draw_jets()
{
  std::cout<<"mG4EveDisplay - draw_jets() begins.."<<std::endl;
  PHG4VtxPoint* vertex =  _truth->GetPrimaryVtx(_truth->GetPrimaryVertexIndex());

    float vx = vertex->get_x();
    float vy = vertex->get_y();
    float vz = vertex->get_z();

  for (JetMap::Iter iter = _jetmap->begin();
         iter != _jetmap->end();
         ++iter) {
      Jet* truejet = iter->second;

      float id    = truejet->get_id();
     // float ncomp = truejet->size_comp();
      float eta   = truejet->get_eta();
      float phi   = truejet->get_phi();
      float e     = truejet->get_e();
      float e_scale = _evedisp->get_jet_e_scale();
      float pt    = truejet->get_pt();
      float pt_threshold = _evedisp->get_jet_pt_threshold();
      if (pt< pt_threshold) continue;
      std::cout<<"Jet "<< id <<" : eta " <<eta << ", phi "<< phi << ", pt "<< pt <<std::endl;

      TEveStraightLineSet* axis= new TEveStraightLineSet("ConeAxis");
      axis->SetLineColor(kGreen);
      axis->SetLineWidth(2);

      TEveVector coneaxis(  (float) Cos((double)phi/CosH((double)eta)),
                        (float) Sin((double)phi/CosH((double)eta)),
                        (float) TanH((double)eta));
      coneaxis *= e/e_scale*length;
      axis->AddLine(vx, vy, vz, coneaxis.fZ, coneaxis.fY, coneaxis.fZ);
      _evemanager->AddElement(axis,_true_jets);

      TEveJetCone* jetcone = new TEveJetCone("JetCone");
      jetcone->SetPickable(kTRUE);
      jetcone->SetCylinder(250., 250.);
      if(jetcone->AddCone(eta,phi, radius)!= -1)
      _evemanager->AddElement(jetcone,_true_jets);
  }
}

void
mG4EveDisplay::draw_event()
{
  add_elements();
  clear();
}


void 
mG4EveDisplay::clear()
{
  _prop->IncDenyDestroy(); // Keep the track propagator from being destroyed
  _true_tracks->DestroyElements();
  _true_jets->DestroyElements();
}
