/*!
        \file mJetEveDisplay.cxx
        \author Sookhyun Lee
        \brief reconstructed jets
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
#include <fun4all/SubsysReco.h>
#include <phool/getClass.h>

// ROOT and EVE includes
#include <TEveManager.h>
#include <TEveJetCone.h>
#include <TEveStraightLineSet.h>
#include <TEveVector.h>
#include <TEveElement.h>


#include <g4hough/SvtxVertexMap.h>		// for maps & tpc
#include <g4jets/JetMap.h>			// for jets
#include <g4jets/Jet.h>

#include <PHEveDisplay.h>

#include <mJetEveDisplay.h>


using boost::bind;
using namespace TMath;

mJetEveDisplay::mJetEveDisplay(boost::shared_ptr<PHEveDisplay> dispin) :
  mPHEveModuleBase(),
  _evedisp(dispin),
  _vtxmap(NULL),
  _jetmap(NULL),
  _reco_jets(NULL),
  radius(0.3),
  length(300)
{
  _evemanager = _evedisp->get_eve_instance();
  _reco_jets = new TEveElementList("Reco jets");
  _evemanager->AddElement(_reco_jets,_evedisp->get_jet_list());
}

mJetEveDisplay::~mJetEveDisplay()
{
}

void
mJetEveDisplay::init(PHCompositeNode* topNode)
{  
}

void 
mJetEveDisplay::init_run(PHCompositeNode* topNode)
{
}

bool
mJetEveDisplay::event(PHCompositeNode* topNode)
{
  std::cout<<"mJetEveDisplay - event.."<<std::endl;
  try
    {
      create_nodes(topNode);
      draw_jets();      
    }
  catch(std::exception& e)
    {
      static bool first(true);
      if( first )
	std::cout << "mJetEveDisplay::event Exception: " << e.what() << std::endl;
      first = false;
    }

  return true;

}

void
mJetEveDisplay::end(PHCompositeNode* topNode)
{
}

void 
mJetEveDisplay::create_nodes(PHCompositeNode* topNode)
{
  _vtxmap = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");
  if (_vtxmap->empty()) std::cout<<"SvtxVertexMap is empty!!! "<<std::endl;
  _jetmap = findNode::getClass<JetMap>(topNode,"AntiKt_Cluster_r03");
  if (!_jetmap) std::cout<<"JetMap node not found!! "<<std::endl;
//  _jetmap->identify();

}

void
mJetEveDisplay::draw_jets()
{
  std::cout<<"mJetEveDisplay - draw_jet() begins.."<<std::endl;


      SvtxVertex* vertex = (_vtxmap->begin()->second);

      float vx = vertex->get_x();
      float vy = vertex->get_y();
      float vz = vertex->get_z();

  for (JetMap::Iter iter = _jetmap->begin();
         iter != _jetmap->end();
         ++iter) {
      Jet* recojet = iter->second;

      float id    = recojet->get_id();
     // float ncomp = recojet->size_comp();
      float eta   = recojet->get_eta();
      float phi   = recojet->get_phi();
      float e     = recojet->get_e();
      float e_scale = _evedisp->get_jet_e_scale();
      float pt    = recojet->get_pt();
      float pt_threshold = _evedisp->get_jet_pt_threshold();     
      if (pt< pt_threshold) continue;

      std::cout<<"Jet "<< id <<" : eta " <<eta << ", phi "<< phi << ", pt " << pt <<std::endl;

      TEveStraightLineSet* axis= new TEveStraightLineSet("ConeAxis");
      axis->SetLineColor(kGreen);
      axis->SetLineWidth(2);

      TEveVector coneaxis(  (float) Cos((double)phi/CosH((double)eta)), 
			(float) Sin((double)phi/CosH((double)eta)),
			(float) TanH((double)eta));
      coneaxis *= e/e_scale*length;
      axis->AddLine(vx, vy, vz, coneaxis.fZ, coneaxis.fY, coneaxis.fZ);
      _evemanager->AddElement(axis,_reco_jets);

      TEveJetCone* jetcone = new TEveJetCone("JetCone");
      jetcone->SetPickable(kTRUE);
      jetcone->SetCylinder(250., 250.);
      if(jetcone->AddCone(eta,phi, radius)!= -1)	
      _evemanager->AddElement(jetcone,_reco_jets);
  }
}



void
mJetEveDisplay::draw_event()
{
  add_elements();
  clear();    
}

void 
mJetEveDisplay::clear()
{
  _reco_jets->DestroyElements();
}
