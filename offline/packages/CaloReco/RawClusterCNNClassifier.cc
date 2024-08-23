/*
I trained a crappy CNN to classify EMCal clusters as photon or not photon. S.Li 8.16.2024
*/


#include "RawClusterCNNClassifier.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

// Tower stuff
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>


#include <ffaobjects/EventHeader.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/onnxlib.h>



#include <iostream>
#include <vector>

RawClusterCNNClassifier::RawClusterCNNClassifier(const std::string &name)
    : SubsysReco(name)
{
}

RawClusterCNNClassifier::~RawClusterCNNClassifier() = default;


int RawClusterCNNClassifier::Init(PHCompositeNode */*topNode*/)
{
    //init the onnx model
    std::string modelPath = "/sphenix/u/shuhang98/core_patch/coresoftware/offline/packages/CaloReco/functional_model.onnx";
    onnxmodule = onnxSession(modelPath);
    return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterCNNClassifier::process_event(PHCompositeNode *topNode)
{
    //get the cluster container
    
    std::string clusterNodeName = "CLUSTERINFO_CEMC";
    RawClusterContainer *clusterContainer = findNode::getClass<RawClusterContainer>(topNode, clusterNodeName);
    if (!clusterContainer)
    {
        std::cout << "Could not locate cluster node " << clusterNodeName << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
    }
    //I trained the model with the info from towerinfo container, raw tower should also work
    std::string towerNodeName = "TOWERINFO_CALIB_CEMC";
    TowerInfoContainer *emcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, towerNodeName);
    if (!emcTowerContainer)
    {
        std::cout << "Could not locate tower node " << towerNodeName << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
    RawClusterContainer::ConstIterator clusterIter;
    for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; ++clusterIter) {
    RawCluster *recoCluster = clusterIter->second;
    //reset the prob inplace
    recoCluster->set_prob(-1);
    CLHEP::Hep3Vector vertex(0, 0, 0);
    CLHEP::Hep3Vector E_vec_cluster_Full = RawClusterUtility::GetEVec(*recoCluster, vertex);
    float ET = E_vec_cluster_Full.perp();
    if (ET < minET)
    {
        continue;
    }
    const RawCluster::TowerConstRange tower_begin_end =
            recoCluster->get_towers();

    int maxtowerE = 0;
    int maxtowerieta = -1;
    int maxtoweriphi = -1;

    for (RawCluster::TowerConstIterator tower_iter = tower_begin_end.first;
             tower_iter != tower_begin_end.second; ++tower_iter) {
          RawTowerDefs::keytype tower_key = tower_iter->first;
         
          // get ieta iphi
          int ix = RawTowerDefs::decode_index2(tower_key); // iphi?
          int iy = RawTowerDefs::decode_index1(tower_key); // ieta I  guess?(S.L.)
          RawTowerDefs::CalorimeterId caloid =
              RawTowerDefs::decode_caloid(tower_key);
          // check if cemc, but I guess they shoul all be anyways lol
          if (caloid != RawTowerDefs::CalorimeterId::CEMC) {
            continue;
          }
          // get the towerinfo key
          unsigned int towerinfokey = TowerInfoDefs::encode_emcal(
              iy, ix); // this is the key for the towerinfo container get the towerinfo
          TowerInfo *towerinfo =
              emcTowerContainer->get_tower_at_key(towerinfokey);
          if (!towerinfo) {
            //should not happen
            std::cout << "No towerinfo for tower key " << towerinfokey
                      << std::endl;
            continue;
          }
          float towerE = towerinfo->get_energy();
          if (towerE > maxtowerE) {
            maxtowerE = towerE;
            maxtowerieta = iy;
            maxtoweriphi = ix;
          }

        }
        //find the N by N tower around the max tower
        std::vector<float> input;
        //resize to inputDimx * inputDimy
        int vectorSize = inputDimx * inputDimy;
        input.resize(vectorSize, 0);

        if (maxtowerE > 0)
        {
          int xlength = int((inputDimx - 1) / 2);
          int ylength = int((inputDimy - 1) / 2);
          if(maxtowerieta - ylength < 0 || maxtowerieta + ylength >= 96)
          {
            continue;
          }
          for (int ieta = maxtowerieta - ylength; ieta <= maxtowerieta + ylength; ieta++)
          {
            for (int iphi = maxtoweriphi - xlength; iphi <= maxtoweriphi + xlength; iphi++)
            {
              int mappediphi = iphi;

              if (mappediphi < 0){ mappediphi += 256; }
              if (mappediphi > 255){ mappediphi -= 256; }
              unsigned int towerinfokey = TowerInfoDefs::encode_emcal(ieta, mappediphi);
              TowerInfo *towerinfo = emcTowerContainer->get_tower_at_key(towerinfokey);
              if (!towerinfo)
              {
                // should not happen
                std::cout << "No towerinfo for tower key " << towerinfokey << std::endl;
                std::cout << "ieta: " << ieta << " iphi: " << mappediphi << std::endl;
                continue;
              }
              int index = (ieta - maxtowerieta + ylength) * inputDimx + iphi - maxtoweriphi + xlength;
              input.at(index) = towerinfo->get_energy();
            }
          }
        }
        std::vector<float> prob = onnxInference(onnxmodule, input, 1, inputDimx, inputDimy, inputDimz, outputDim);
        //std::cout << "new prob: " << prob[0] << "ET: " << ET << " original prob: " << recoCluster->get_prob() << std::endl;
        //inplace change for the prob for now
        recoCluster->set_prob(prob[0]);

    }

    return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterCNNClassifier::End(PHCompositeNode */*topNode*/)
{
    delete onnxmodule;
    return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterCNNClassifier::Print(const std::string &what) const
{
    std::cout << "RawClusterCNNClassifier::Print(const std::string &what) const Printing info for " << what << std::endl;
    return;
}