
/*!
	\file mCaloEveDisplay.cxx
        \author Sookhyun Lee
        \brief reconstructed energy clusters from cemc/hcalin/hcalout
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
#include <phool/PHIODataNode.h>
//#include <PHPoint.h>
#include <phool/getClass.h>

// ROOT and EVE includes
#include <TEveManager.h>
#include <TEveRGBAPalette.h>
#include <TEveBoxSet.h>
#include <TEveElement.h>
#include <TEveCaloData.h>
#include <TEveCalo.h>
#include <TH2F.h>

// Calo includes
#include <g4cemc/RawTowerContainer.h> 
#include <g4cemc/RawTower.h>
#include <g4cemc/RawTowerGeomContainerv1.h>
#include <g4cemc/RawTowerGeomv3.h>
#include <g4cemc/RawClusterContainer.h>
#include <g4cemc/RawCluster.h>

#include <PHEveDisplay.h>

#include <mCaloEveDisplay.h>


using boost::bind;

mCaloEveDisplay::mCaloEveDisplay(boost::shared_ptr<PHEveDisplay> dispin) :
  mPHEveModuleBase(),
  _evedisp(dispin),
  _cemc_list(NULL),
  _hcalin_list(NULL),
  _hcalout_list(NULL),
  _cemc_boxset(NULL),
  _hcalin_boxset(NULL),
  _hcalout_boxset(NULL),
  _pal(new TEveRGBAPalette(0,130)),
  _cemc_clusters(NULL),
  _hcalin_clusters(NULL),
  _hcalout_clusters(NULL),
  _cemc_towers(NULL),
  _hcalin_towers(NULL),
  _hcalout_towers(NULL),
  _cemc_towergeo(NULL),
  _hcalin_towergeo(NULL),
  _hcalout_towergeo(NULL)
{
  _evemanager = _evedisp->get_eve_instance();

  _cemc_list = new TEveElementList("CEMC");
  _hcalin_list = new TEveElementList("HCALIN");
  _hcalout_list = new TEveElementList("HCALOUT");

  _evemanager->AddElement(_cemc_list,_evedisp->get_calo_list());
  _evemanager->AddElement(_hcalin_list,_evedisp->get_calo_list());
  _evemanager->AddElement(_hcalout_list,_evedisp->get_calo_list());

}

mCaloEveDisplay::~mCaloEveDisplay()
{
}

void
mCaloEveDisplay::init(PHCompositeNode* topNode)
{  
}

void 
mCaloEveDisplay::init_run(PHCompositeNode* topNode)
{
}

bool
mCaloEveDisplay::event(PHCompositeNode* topNode)
{
  std::cout<<"mCaloEveDisplay - event.."<<std::endl;
  try
    {
      create_nodes(topNode);
      draw_clusters(true, true, true);
    }
  catch(std::exception& e)
    {
      static bool first(true);
      if( first )
	std::cout << "mCaloEveDisplay::event Exception: " << e.what() << std::endl;
      first = false;
    }

  return true;

}

void
mCaloEveDisplay::end(PHCompositeNode* topNode)
{
}

void 
mCaloEveDisplay::create_nodes(PHCompositeNode* topNode)
{
  _cemc_clusters = findNode::getClass<RawClusterContainer>(topNode,"CLUSTER_CEMC");
  _cemc_towers = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_CEMC");
  _cemc_towergeo = findNode::getClass<RawTowerGeomContainer>(topNode,"TOWERGEOM_CEMC");
  if(!(_cemc_clusters && _cemc_towers && _cemc_towergeo))
  {
    std::cerr << "Error: Can't find CEMC Clusters & Towers" << std::endl;
  }

  _hcalin_clusters = findNode::getClass<RawClusterContainer>(topNode,"CLUSTER_HCALIN");
  _hcalin_towers = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_HCALIN");
  _hcalin_towergeo = findNode::getClass<RawTowerGeomContainer>(topNode,"TOWERGEOM_HCALIN");
  if(!(_hcalin_clusters && _hcalin_towers && _hcalin_towergeo))
  {
    std::cerr << "Error: Can't find HCALIN Clusters & Towers" << std::endl;
  }
  _hcalout_clusters = findNode::getClass<RawClusterContainer>(topNode,"CLUSTER_HCALOUT");
  _hcalout_towers = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_HCALOUT");
  _hcalout_towergeo = findNode::getClass<RawTowerGeomContainer>(topNode,"TOWERGEOM_HCALOUT");
  if(!(_hcalout_clusters && _hcalout_towers && _hcalout_towergeo))
  {
    std::cerr << "Error: Can't find HCALOUT Cluster & Towers" << std::endl;
  }
  std::cout<<"mCaloEveDisplay - nodes created.."<<std::endl;
}


void
mCaloEveDisplay::draw_clusters(bool is_cemc_on, bool is_hcalin_on, bool is_hcalout_on)
{
  // box set approach
  unsigned int iclust = 0;
  float energy = 0.;
  float e_threshold = _evedisp->get_calo_e_threshold();
  float center_x = 0.;
  float center_y = 0.; 
  float center_z = 0.;
  float size = 0.;
  float scale = 0.;
  RawCluster* cluster;
  RawTower* tower;
  RawTowerGeom* towergeo;

  if(is_cemc_on)
  {
    _cemc_boxset = new TEveBoxSet("CEMC_BOXSET");
    _cemc_boxset->SetPalette(_pal);
    _cemc_boxset->Reset(TEveBoxSet::kBT_AABox,kFALSE,64);

    for (iclust = 0; iclust < _cemc_clusters->size(); iclust++)  
    {
      cluster = _cemc_clusters->getCluster(iclust); 
      RawCluster::TowerConstRange begin_end = cluster->get_towers();
      for (RawCluster::TowerConstIterator iter = begin_end.first;
	iter != begin_end.second;
	++iter)
      {
	tower = _cemc_towers->getTower(iter->first);
	towergeo = _cemc_towergeo->get_tower_geometry(iter->first);

	energy = tower->get_energy();
        if (energy < e_threshold) continue;

        center_x = towergeo->get_center_x();
        center_y = towergeo->get_center_y();
        center_z = towergeo->get_center_z();
        std::cout<<"tower "<<iter->first<<": energy "<< energy
        <<", center x " << center_x
        <<", center y " << center_y
        <<", center z " << center_z <<std::endl;
        size = 4.0;
        scale = TMath::Log(energy)+4;//1MeV :1, 1GeV :4 
        _cemc_boxset->AddBox(center_x-scale*size/2.,
                             center_y-scale*size/2.,
                             center_z-scale*size/2.,
                              scale*size,
                              scale*size,
                              scale*size);
        _cemc_boxset->DigitValue(tower->get_energy()*1000.);
      }         
    }
    _evemanager->AddElement(_cemc_boxset,_cemc_list);
  }
 
  if(is_hcalin_on)
  {
    _hcalin_boxset = new TEveBoxSet("HCALIN_BOXSET");
    _hcalin_boxset->SetPalette(_pal);
    _hcalin_boxset->Reset(TEveBoxSet::kBT_AABox,kFALSE,64);

    for (iclust = 0; iclust < _hcalin_clusters->size(); iclust++)
    {
      cluster = _hcalin_clusters->getCluster(iclust);
      std::cout<<"cluster " << iclust <<" found." <<std::endl;
      RawCluster::TowerConstRange begin_end = cluster->get_towers();
      std::cout<<"tower loop begins. "<<std::endl;
      for (RawCluster::TowerConstIterator iter = begin_end.first;
        iter != begin_end.second;
        ++iter)
      {
        tower = _hcalin_towers->getTower(iter->first);
        towergeo = _hcalin_towergeo->get_tower_geometry(iter->first);

	energy = tower->get_energy();
        if (energy < e_threshold) continue;
	center_x = towergeo->get_center_x();
	center_y = towergeo->get_center_y();
	center_z = towergeo->get_center_z();
        std::cout<<"tower "<<iter->first<<": energy "<< energy
        <<", center x " << center_x
	<<", center y " << center_y
	<<", center z " << center_z <<std::endl;
	size = 4.0;
	scale = TMath::Log(energy)+4;//1MeV :1, 1GeV :4 
        _hcalin_boxset->AddBox(center_x-scale*size/2.,
                               center_y-scale*size/2.,
                               center_z-scale*size/2.,
				scale*size,
				scale*size,
				scale*size);
        _hcalin_boxset->DigitValue(tower->get_energy()*1000.);
      }
    }
    _evemanager->AddElement(_hcalin_boxset,_hcalin_list);
  }

  if(is_hcalout_on)
  {
    _hcalout_boxset = new TEveBoxSet("HCALOUT_BOXSET");
    _hcalout_boxset->SetPalette(_pal);
    _hcalout_boxset->Reset(TEveBoxSet::kBT_AABox,kFALSE,64);

    for (iclust = 0; iclust < _hcalout_clusters->size(); iclust++)
    {
      cluster = _hcalout_clusters->getCluster(iclust);
      std::cout<<"cluster " << iclust <<" found." <<std::endl;
      RawCluster::TowerConstRange begin_end = cluster->get_towers();
      std::cout<<"tower loop begins. "<<std::endl;
      for (RawCluster::TowerConstIterator iter = begin_end.first;
        iter != begin_end.second;
        ++iter)
      {
        tower = _hcalout_towers->getTower(iter->first);
        towergeo = _hcalout_towergeo->get_tower_geometry(iter->first);

        energy = tower->get_energy();
        if (energy < e_threshold) continue;
        center_x = towergeo->get_center_x();
        center_y = towergeo->get_center_y();
        center_z = towergeo->get_center_z();
        std::cout<<"tower "<<iter->first<<": energy "<< energy
        <<", center x " << center_x
        <<", center y " << center_y
        <<", center z " << center_z <<std::endl;
        size = 4.0;
        scale = TMath::Log(energy)+4;//1MeV :1, 1GeV :4 
        _hcalout_boxset->AddBox(center_x-scale*size/2.,
                               center_y-scale*size/2.,
                               center_z-scale*size/2.,
                                scale*size,
                                scale*size,
                                scale*size);
        _hcalout_boxset->DigitValue(tower->get_energy()*1000.);
      }
    }
    _evemanager->AddElement(_hcalout_boxset,_hcalout_list);
  }

}

void
mCaloEveDisplay::draw_event()
{
  add_elements();
  clear();
}

void 
mCaloEveDisplay::clear()
{
  _cemc_list->DestroyElements();
  _cemc_boxset->DestroyElements();
  _hcalin_list->DestroyElements();
  _hcalin_boxset->DestroyElements();
  _hcalout_list->DestroyElements();
  _hcalout_boxset->DestroyElements();
}
