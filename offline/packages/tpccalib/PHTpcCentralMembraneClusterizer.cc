#include "PHTpcCentralMembraneClusterizer.h"

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

/// Tracking includes
#include <trackbase/ActsGeometry.h>
#include <trackbase/CMFlashClusterContainerv1.h>
#include <trackbase/CMFlashClusterv3.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/ActsTransformations.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TVector3.h>

#include <boost/format.hpp>

#include <cmath>
#include <string>

namespace
{

  // stream acts vector3
  [[maybe_unused]] std::ostream &operator<<(std::ostream &out, const TVector3 &v)
  {
    out << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return out;
  }

}  // namespace

//____________________________________________________________________________..
PHTpcCentralMembraneClusterizer::PHTpcCentralMembraneClusterizer(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int PHTpcCentralMembraneClusterizer::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);

  if (_histos)
  {
    m_histogramfile.reset(new TFile(m_histogramfilename.c_str(), "RECREATE"));
    m_histogramfile->cd();
    henergy = new TH1F("henergy", "cluster energy", 200, 0, 2000);
    hxy = new TH2F("hxy", "cluster x:y", 800, -100, 100, 800, -80, 80);
    hz = new TH1F("hz", "cluster z", 220, -2, 2);

    hz_pos = new TH1F("hz_pos", "cluster z>0", 400, 0, 110);
    hz_neg = new TH1F("hz_neg", "cluster z<0", 400, -110, 0);

    hClustE[0] = new TH1F("hRawClusterEnergy", "Cluster Energy Before Merging;E[?]", 200, 0, 2000);
    hClustE[1] = new TH1F("hMatchedClusterEnergy", "Pair Cluster Energy After Merging;E[?]", 200, 0, 2000);
    hClustE[2] = new TH1F("hSoloClusterEnergy", "Lone Cluster Energy After Merging;E[?]", 200, 0, 2000);

    hDist = new TH1F("hDist", "3D distance to nearby clusters on same padrow;dist[cm]", 100, -1, 10);
    hDistRow = new TH2F("hDistRow", "phi distance to nearby clusters vs (lower)row;dist[rad];padrow", 100, -0.001, 0.01, 60, -0.5, 59.5);
    hDist2 = new TH1F("hDist2", "phi distance to nearby clusters on same padrow;dist[rad]", 100, -0.001, 0.01);
    hDistRowAdj = new TH2F("hDistRowAdj", "phi distance to nearby clusters vs (lower)row;dist[rad];(lower) padrow", 100, -0.001, 0.01, 60, -0.5, 59.5);
    hDist2Adj = new TH1F("hDist2Adj", "phi distance to nearby clusters on adjacent padrow;dist[rad]", 100, -0.001, 0.01);
  }

  hrPhi_reco_petalModulo_pos = new TH2F("hrPhi_reco_petalModulo_pos", "Reco R vs #phi Z > 0;#phi;R (cm)", 50, 0.0, M_PI / 9, 400, 20, 78);
  hrPhi_reco_petalModulo_neg = new TH2F("hrPhi_reco_petalModulo_neg", "Reco R vs #phi;#phi Z < 0;R (cm)", 50, 0.0, M_PI / 9, 400, 20, 78);

  for (int i = 0; i < 48; i++)
  {
    hphi_reco_pos[i] = new TH1F(
        (boost::format("hphi_reco_pos_layer%02d") % (7 + i)).str().c_str(),
        (boost::format("Reco #phi for Layer %02d Z > 0;#phi;counts") % (7 + i)).str().c_str(),
        50, 0.0, M_PI / 9);
    hphi_reco_neg[i] = new TH1F(
        (boost::format("hphi_reco_neg_layer%02d") % (7 + i)).str().c_str(),
        (boost::format("Reco #phi for Layer %02d Z < 0;#phi;counts") % (7 + i)).str().c_str(),
        50, 0.0, M_PI / 9);

    if (i < 47)
    {
      hphi_reco_pair_pos[i] = new TH1F(
          (boost::format("hphi_reco_pair_pos_layers%02d_%02d") % (7 + i) % (8 + i)).str().c_str(),
          (boost::format("Reco #phi for Layers %02d and %02d Z > 0;#phi;counts") % (7 + i) % (8 + i)).str().c_str(),
          50, 0.0, M_PI / 9);
      hphi_reco_pair_neg[i] = new TH1F(
          (boost::format("hphi_reco_pair_neg_layers%02d_%02d") % (7 + i) % (8 + i)).str().c_str(),
          (boost::format("Reco #phi for Layers %02d and %02d Z < 0;#phi;counts") % (7 + i) % (8 + i)).str().c_str(),
          50, 0.0, M_PI / 9);
    }
  }

  return ret;
}

//____________________________________________________________________________..
int PHTpcCentralMembraneClusterizer::process_event(PHCompositeNode *topNode)
{
  // local coord conversion below
  auto tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  if (!tgeometry)
  {
    std::cout << PHWHERE << "No Acts geometry on node tree. Can't  continue." << std::endl;
  }

  if (Verbosity() > 0)
  {
    std::cout << std::endl
              << "original size of cluster map: " << _cluster_map->size() << std::endl;
  }

  std::vector<TVector3> pos;       // position vector in cartesian
  std::vector<int> layer;          // cluster layer number
  std::vector<unsigned int> side;  // cluster side
  std::vector<int> i_pair;         // vector for pair matching
  std::vector<float> energy;       // vector for energy values of clusters
  std::vector<bool> isAcrossGap;
  int nTpcClust = 0;

  double mean_z_content_plus = 0.0;
  double mean_z_content_minus = 0.0;

  // first loop over clusters to make mod phi histograms of each layer and each pair of layers
  for (const auto &hitsetkey : _cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    auto clusRange = _cluster_map->getClusters(hitsetkey);
    for (auto clusiter = clusRange.first;
         clusiter != clusRange.second; ++clusiter)
    {
      const auto &[cluskey, cluster] = *clusiter;
      auto glob = tgeometry->getGlobalPosition(cluskey, cluster);
      TVector3 tmp_pos(glob(0), glob(1), glob(2));

      unsigned int cluster_layer = TrkrDefs::getLayer(cluskey);

      double phi = tmp_pos.Phi();

      double phiMod = phi;

      // make sure phiMod is between 0 and pi/9
      while (phiMod < 0.0 || phiMod > M_PI / 9)
      {
        if (phiMod < 0.0)
        {
          phiMod += M_PI / 9;
        }
        else if (phiMod > M_PI / 9)
        {
          phiMod -= M_PI / 9;
        }
      }

      if (tmp_pos.Z() >= 0.0)
      {
        hz_pos->Fill(tmp_pos.Z());

        // mean_z_content_plus += tmp_pos.Z();

        hrPhi_reco_petalModulo_pos->Fill(phiMod, tmp_pos.Perp());
        hphi_reco_pos[cluster_layer - 7]->Fill(phiMod);
        // for layer pairs, if last layer can only go in layer 53-54 pair, if first layer can only go in layer 7-8 pair
        if (cluster_layer < 54)
        {
          hphi_reco_pair_pos[cluster_layer - 7]->Fill(phiMod);
        }
        if (cluster_layer > 7)
        {
          hphi_reco_pair_pos[cluster_layer - 8]->Fill(phiMod);
        }
      }
      else
      {
        hz_neg->Fill(tmp_pos.Z());

        // mean_z_content_minus += tmp_pos.Z();

        hrPhi_reco_petalModulo_neg->Fill(phiMod, tmp_pos.Perp());
        hphi_reco_neg[cluster_layer - 7]->Fill(phiMod);
        // for layer pairs, if last layer can only go in layer 53-54 pair, if first layer can only go in layer 7-8 pair
        if (cluster_layer < 54)
        {
          hphi_reco_pair_neg[cluster_layer - 7]->Fill(phiMod);
        }
        if (cluster_layer > 7)
        {
          hphi_reco_pair_neg[cluster_layer - 8]->Fill(phiMod);
        }
      }
    }
  }

  for (int i = 1; i < hz_pos->GetNbinsX(); i++)
  {
    mean_z_content_plus += hz_pos->GetBinContent(i);
  }
  for (int i = 1; i < hz_neg->GetNbinsX(); i++)
  {
    mean_z_content_minus += hz_neg->GetBinContent(i);
  }

  mean_z_content_plus = mean_z_content_plus / hz_pos->GetNbinsX();
  mean_z_content_minus = mean_z_content_minus / hz_neg->GetNbinsX();

  // use peak in z distributions to identify if there is a CM Flash. Peak should be >20% of entries (needs to be tuned)
  if (hz_pos->GetMaximum() < 5 * mean_z_content_plus || hz_neg->GetMaximum() < 5 * mean_z_content_minus)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // loop over histos for each pair of layers, count number of bins above threshold
  // for a given layer, layer pair with higher average value above threshold will be better match for meta-clustering
  for (int i = 0; i < 47; i++)
  {
    for (int j = 1; j < hphi_reco_pair_pos[i]->GetNbinsX(); j++)
    {
      if (hphi_reco_pair_pos[i]->GetBinContent(j) > (m_metaClusterThreshold + (mean_z_content_plus / 18.0)))
      {
        nPairAbove_pos[i]++;
        pairAboveContent_pos[i] += (hphi_reco_pair_pos[i]->GetBinContent(j) - (m_metaClusterThreshold + (mean_z_content_plus / 18.0)));
      }
      if (hphi_reco_pair_neg[i]->GetBinContent(j) > (m_metaClusterThreshold + (mean_z_content_minus / 18.0)))
      {
        nPairAbove_neg[i]++;
        pairAboveContent_neg[i] += (hphi_reco_pair_neg[i]->GetBinContent(j) - (m_metaClusterThreshold + (mean_z_content_minus / 18.0)));
      }
    }
  }

  // loop over clusters again
  for (const auto &hitsetkey : _cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    auto clusRange = _cluster_map->getClusters(hitsetkey);
    for (auto clusiter = clusRange.first;
         clusiter != clusRange.second; ++clusiter)
    {
      ++m_total_clusters;

      const auto &[cluskey, cluster] = *clusiter;
      auto glob = tgeometry->getGlobalPosition(cluskey, cluster);
      // std::cout << "PHTpcCentralMembraneClusterizer::process_event - key: " << cluskey << "z: " << glob.z() << std::endl;

      float x = glob(0);
      float y = glob(1);
      float z = glob(2);

      if (Verbosity() > 0)
      {
        unsigned int lyr = TrkrDefs::getLayer(cluskey);
        unsigned short cluster_side = TpcDefs::getSide(cluskey);
        std::cout << " z " << z << " side " << cluster_side << " layer " << lyr << " Adc " << cluster->getAdc() << " x " << x << " y " << y << std::endl;
      }

      TVector3 tmp_pos(x, y, z);

      unsigned int lay = TrkrDefs::getLayer(cluskey);

      double phi = tmp_pos.Phi();
      double phiMod = phi;

      while (phiMod < 0.0 || phiMod > M_PI / 9)
      {
        if (phiMod < 0.0)
        {
          phiMod += M_PI / 9;
        }
        else if (phiMod > M_PI / 9)
        {
          phiMod -= M_PI / 9;
        }
      }

      int phiBin = -1;

      bool aboveThreshold = true;

      // Find which phiMod bin cluster is in (with given layer) and if counts below threshold don't use
      if (z >= 0)
      {
        phiBin = hphi_reco_pos[lay - 7]->GetXaxis()->FindBin(phiMod);
        if (hphi_reco_pos[lay - 7]->GetBinContent(phiBin) < m_moduloThreshold)
        {
          aboveThreshold = false;
        }
      }
      else
      {
        phiBin = hphi_reco_neg[lay - 7]->GetXaxis()->FindBin(phiMod);
        if (hphi_reco_neg[lay - 7]->GetBinContent(phiBin) < m_moduloThreshold)
        {
          aboveThreshold = false;
        }
      }

      if (!aboveThreshold)
      {
        continue;
      }

      if (cluster->getAdc() < _min_adc_value)
      {
        continue;
      }
      if (std::abs(z) < _min_z_value)
      {
        continue;
      }

      // require that z be within 1cm of peak in z distributions (separate for each side)
      if ((z >= 0 && std::abs(z - hz_pos->GetBinCenter(hz_pos->GetMaximumBin())) > 1.0) || (z < 0 && std::abs(z - hz_neg->GetBinCenter(hz_neg->GetMaximumBin())) > 1.0))
      {
        continue;
      }

      ++m_accepted_clusters;

      i_pair.push_back(-1);
      energy.push_back(cluster->getAdc());
      nTpcClust++;
      pos.emplace_back(x, y, z);
      layer.push_back((int) (TrkrDefs::getLayer(cluskey)));
      side.push_back(TpcDefs::getSide(cluskey));
      isAcrossGap.push_back(false);
      if (Verbosity() > 0)
      {
        std::cout << ":\t" << x << "\t" << y << "\t" << z << std::endl;
      }
    }
  }

  if (_histos)
  {
    // fill evaluation histograms
    for (int i = 0; i < nTpcClust; ++i)
    {
      for (int j = i + 1; j < nTpcClust; j++)
      {
        // must match clusters that are on the same side
        if (side[i] != side[j])
        {
          continue;
        }

        if (layer[i] == layer[j])
        {
          const TVector3 delta = pos[i] - pos[j];
          const float dphi = std::abs(pos[i].DeltaPhi(pos[j]));
          hDist->Fill(delta.Mag());
          hDist2->Fill(dphi);
          hDistRow->Fill(dphi, layer[i]);
        }

        if (std::abs(layer[i] - layer[j]) == 1)
        {
          // match those centers to the known/expected stripe positions
          const TVector3 delta = pos[i] - pos[j];
          const float dphi = std::abs(pos[i].DeltaPhi(pos[j]));
          hDist2Adj->Fill(dphi);
          hDistRowAdj->Fill(dphi, layer[i]);
        }
      }
    }
  }

  // now for each cluster, find its nearest partner on an adjacent row

  /* TODO: check with Tony. Should be enough to loop for j=i+1 */

  const float maxphidist = 0.003;  // as read off the plots.
  for (int i = 0; i < nTpcClust; i++)
  {
    int layerMatch = -1;

    int nRowsMatch = 2;
    if (layer[i] >= 39)
    {
      nRowsMatch = 4;
    }
    else if (layer[i] >= 23)
    {
      nRowsMatch = 3;
    }

    if (pos[i].Z() >= 0)
    {
      if (layer[i] == 7)
      {
        // if layer is 7, can only get one layer match
        if (nPairAbove_pos[layer[i] - 7] >= nRowsMatch)
        {
          layerMatch = layer[i] + 1;
        }
      }
      else if (layer[i] == 54)
      {
        // if layer is 8, can only get one layer match
        if (nPairAbove_pos[layer[i] - 8] >= nRowsMatch)
        {
          layerMatch = layer[i] - 1;
        }
      }
      else
      {
        // if both pairs of rows have enough bins above threshold, match to pair with higher average
        if (nPairAbove_pos[layer[i] - 7] >= nRowsMatch && nPairAbove_pos[layer[i] - 8] >= nRowsMatch)
        {
          if (pairAboveContent_pos[layer[i] - 7] / nPairAbove_pos[layer[i] - 7] > pairAboveContent_pos[layer[i] - 8] / nPairAbove_pos[layer[i] - 8])
          {
            layerMatch = layer[i] + 1;
          }
          else
          {
            layerMatch = layer[i] - 1;
          }
          // otherwise just use the one that is above threshold
        }
        else if (nPairAbove_pos[layer[i] - 7] >= nRowsMatch)
        {
          layerMatch = layer[i] + 1;
        }
        else if (nPairAbove_pos[layer[i] - 8] >= nRowsMatch)
        {
          layerMatch = layer[i] - 1;
        }
      }
    }
    else
    {
      if (layer[i] == 7)
      {
        if (nPairAbove_neg[layer[i] - 7] >= nRowsMatch)
        {
          layerMatch = layer[i] + 1;
        }
      }
      else if (layer[i] == 54)
      {
        if (nPairAbove_neg[layer[i] - 8] >= nRowsMatch)
        {
          layerMatch = layer[i] - 1;
        }
      }
      else
      {
        if (nPairAbove_neg[layer[i] - 7] >= nRowsMatch && nPairAbove_neg[layer[i] - 8] >= nRowsMatch)
        {
          if (pairAboveContent_neg[layer[i] - 7] / nPairAbove_neg[layer[i] - 7] > pairAboveContent_neg[layer[i] - 8] / nPairAbove_neg[layer[i] - 8])
          {
            layerMatch = layer[i] + 1;
          }
          else
          {
            layerMatch = layer[i] - 1;
          }
        }
        else if (nPairAbove_neg[layer[i] - 7] >= nRowsMatch)
        {
          layerMatch = layer[i] + 1;
        }
        else if (nPairAbove_neg[layer[i] - 8] >= nRowsMatch)
        {
          layerMatch = layer[i] - 1;
        }
      }
    }

    // if the match is default and the layer is on the edge of a module, identify it as being across the gap
    //     if(layerMatch == -1 && (layer[i] == 22 || layer[i] == 23 || layer[i] == 38 || layer[i] == 39 || layer[i] == 7) ) isAcrossGap[i] = true;
    if (layer[i] == 22 || layer[i] == 23 || layer[i] == 38 || layer[i] == 39 || layer[i] == 7)
    {
      isAcrossGap[i] = true;
    }

    float bestphidist = maxphidist;
    for (int j = 0; j < nTpcClust; j++)
    {
      // radial match must be the correct layer identified above
      if (layer[j] != layerMatch)
      {
        continue;
      }

      // redundant to the 'adjacent row' check:  if (i==j) continue; //can't match yourself.
      // must match to an adjacent row.
      if (std::abs(layer[i] - layer[j]) != 1)
      {
        continue;
      }

      // setting agross gap for ones where match was found
      // if( (layer[i] == 22 && layer[j] == 23) || (layer[i] == 23 && layer[j] == 22) || (layer[i] == 38 && layer[j] == 39) || (layer[i] == 39 && layer[j] == 38) ) isAcrossGap[i] = true;

      // must match clusters that are on the same side
      if (side[i] != side[j])
      {
        continue;
      }

      // must match clusters that are close to each other in Z
      if (std::abs(pos[i].Z() - pos[j].Z()) > 0.5)
      {
        continue;
      }

      const float newphidist = std::abs(pos[i].DeltaPhi(pos[j]));
      if (newphidist < bestphidist)
      {
        i_pair[i] = j;
        bestphidist = newphidist;
      }
    }
  }

  // check to see if the cluster pairs each match each other
  std::vector<bool> goodPair;
  bool allGood = true;
  int nGood = 0;

  for (int i = 0; i < nTpcClust; i++)
  {
    int myPair = i_pair[i];
    int itsPair = myPair < 0 ? -1 : i_pair[myPair];
    if (i != itsPair)
    {
      if (Verbosity())
      {
        std::cout << "PHTpcCentralMembraneClusterizer::process_event -"
                  << " i: " << i
                  << " myPair: " << myPair
                  << " itsPair: " << itsPair
                  << std::endl;
      }

      goodPair.push_back(false);
      allGood = false;
    }
    else
    {
      if (i < myPair)
      {
        ++nGood;
      }
      goodPair.push_back(true);
    }
  }

  if (Verbosity())
  {
    if (allGood)
    {
      std::cout << "PHTpcCentralMembraneClusterizer::process_event - all pairs are good" << std::endl;
    }
    else
    {
      std::cout << "PHTpcCentralMembraneClusterizer::process_event - nGood: " << nGood << " out of " << nTpcClust / 2 << std::endl;
    }
  }

  // build the weighted cluster centers
  //==========================
  std::vector<float> aveenergy;
  std::vector<TVector3> avepos;
  std::vector<unsigned int> nclusters;
  std::vector<bool> isREdge;
  std::vector<std::pair<int, int> > pairNums;

  //  int nR2 = 0;
  // int nR3 = 0;

  //  double R2AveE = 0.0;
  // double R3AveE = 0.0;

  std::pair<int, int> tmp_pair;

  for (int i = 0; i < nTpcClust; ++i)
  {
    if (_histos)
    {
      hClustE[0]->Fill(energy[i]);
    }

    if (goodPair[i])
    {
      if (i_pair[i] > i)
      {
        if (_histos)
        {
          hClustE[1]->Fill(energy[i] + energy[i_pair[i]]);
        }

        aveenergy.push_back(energy[i] + energy[i_pair[i]]);

        //	if(

        // The pads measure phi and z accurately
        // They do not measure R! It is taken as the center of the padrow
        // The x and y values are derived from phi and R. Not normally a problem, since tracks cross entire padrow
        // CM flash clusters have limited radial extent, do not necessarily cross padrows completely - then their nominal R is wrong.
        // So:
        //     Get phi from the energy weighted cluster phi values.
        //     Get R from the procedure below

        // Get phi and z centroid
        double avePhi = (pos[i].Phi() * energy[i] + pos[i_pair[i]].Phi() * energy[i_pair[i]]) * (1. / (energy[i] + energy[i_pair[i]]));
        double aveZ = (pos[i].Z() * energy[i] + pos[i_pair[i]].Z() * energy[i_pair[i]]) * (1. / (energy[i] + energy[i_pair[i]]));

        // Single padrow CM flash clusters, R position is not well defined because cluster is smaller than padrow
        // 2-cluster case: Weighting by padrow center radius is not correct because distribution does not fill padrow (needs to be approximately linearly)
        //      Use ratio of component cluster energies to estimate number of sigmas at row boundary
        float efrac = energy[i] / (energy[i] + energy[i_pair[i]]);

        PHG4TpcCylinderGeom *layergeom1 = _geom_container->GetLayerCellGeom(layer[i]);
        double rad1 = layergeom1->get_radius();
        PHG4TpcCylinderGeom *layergeom2 = _geom_container->GetLayerCellGeom(layer[i_pair[i]]);
        double rad2 = layergeom2->get_radius();
        // matching done correctly now, so distance between layers should be directly calculable
        double layer_dr = std::abs(rad1 - rad2);

        double rad_lyr_boundary = rad1 + layer_dr / 2.0;

        // We have to (temporarily) use distortion corrected cluster positions to determine which stripe this came from
        Acts::Vector3 dist_pos(pos[i].X(), pos[i].Y(), pos[i].Z());
        if (_dcc)
        {
          dist_pos = _distortionCorrection.get_corrected_position(dist_pos, _dcc);
        }
        double dist_r = sqrt(dist_pos[0] * dist_pos[0] + dist_pos[1] * dist_pos[1]);
        double cmclus_dr = _cmclus_dr_outer;

        if (dist_r < 41.0)
        {
          // if across boundary, use average
          if (rad2 >= 41.0)
          {
            cmclus_dr = 0.5 * (_cmclus_dr_inner + _cmclus_dr_mid);
          }
          else
          {
            cmclus_dr = _cmclus_dr_inner;
          }
        }
        else if (dist_r >= 41.0 && dist_r < 58.0)
        {
          // if across boundary, use average
          if (rad2 >= 58.0)
          {
            cmclus_dr = 0.5 * (_cmclus_dr_mid + _cmclus_dr_outer);
          }
          else
          {
            cmclus_dr = _cmclus_dr_mid;
          }
        }

        // Use radial width of stripe and efrac to determine where radius at center of distribution must be
        double aveR = rad_lyr_boundary - efrac * cmclus_dr + cmclus_dr / 2.0;

        if (Verbosity() > 0)
        {
          std::cout << " efrac " << efrac << " _cmclus_dr " << cmclus_dr << " rad_lyr_boundary " << rad_lyr_boundary << " aveR " << aveR
                    << " layer i " << layer[i] << " R i " << rad1 << " layer i_pair " << layer[i_pair[i]] << " R i_pair " << rad2 << " layer_dr " << layer_dr << std::endl;
        }

        TVector3 temppos(aveR * cos(avePhi), aveR * sin(avePhi), aveZ);
        avepos.push_back(temppos);
        nclusters.push_back(2);
        if (isAcrossGap[i] && isAcrossGap[i_pair[i]])
        {
          isREdge.push_back(true);
        }
        else
        {
          isREdge.push_back(false);
        }

        tmp_pair.first = i;
        tmp_pair.second = i_pair[i];
        pairNums.push_back(tmp_pair);

        if (Verbosity() > 0)
        {
          std::cout << " layer i " << layer[i] << " energy " << energy[i] << " pos i " << pos[i].X() << "  " << pos[i].Y() << "  " << pos[i].Z()
                    << " layer i_pair " << layer[i_pair[i]] << " energy i_pair " << energy[i_pair[i]]
                    << " pos i_pair " << pos[i_pair[i]].X() << "  " << pos[i_pair[i]].Y() << "  " << pos[i_pair[i]].Z()
                    << " reco pos " << temppos.X() << "  " << temppos.Y() << "  " << temppos.Z()
                    << std::endl;
        }
      }
    }
    else
    {
      if (_histos)
      {
        hClustE[2]->Fill(energy[i]);
      }
      // These single cluster cases have good phi, but do not have a good radius centroid estimate - may want to skip them, record nclusters and identify if across gap
      //      if(layer[i] == 7) isAcrossGap[i] = true;

      aveenergy.push_back(energy[i]);
      avepos.push_back(pos[i]);
      nclusters.push_back(1);
      isREdge.push_back(isAcrossGap[i]);
      tmp_pair.first = i;
      tmp_pair.second = -1;
      pairNums.push_back(tmp_pair);
    }
  }

  // Loop over the vectors and put the clusters on the node tree
  //==============================================
  if (Verbosity() > 1)
  {
    std::cout << " vector size is " << avepos.size() << std::endl;
  }

  for (unsigned int iv = 0; iv < avepos.size(); ++iv)
  {
    auto cmfc = new CMFlashClusterv3();

    cmfc->setX(avepos[iv].X());
    cmfc->setY(avepos[iv].Y());
    cmfc->setZ(avepos[iv].Z());
    cmfc->setAdc(aveenergy[iv]);
    cmfc->setNclusters(nclusters[iv]);
    cmfc->setIsRGap(isREdge[iv]);

    int pair1 = pairNums[iv].first;
    int pair2 = pairNums[iv].second;

    cmfc->setX1(pos[pair1].X());
    cmfc->setY1(pos[pair1].Y());
    cmfc->setZ1(pos[pair1].Z());
    cmfc->setLayer1(layer[pair1]);
    cmfc->setAdc1(energy[pair1]);
    if (pair2 != -1)
    {
      cmfc->setX2(pos[pair2].X());
      cmfc->setY2(pos[pair2].Y());
      cmfc->setZ2(pos[pair2].Z());
      cmfc->setLayer2(layer[pair2]);
      cmfc->setAdc2(energy[pair2]);
    }
    else
    {
      cmfc->setX2(0);
      cmfc->setY2(0);
      cmfc->setZ2(0);
      cmfc->setLayer2(0);
      cmfc->setAdc2(0);
    }

    _corrected_CMcluster_map->addClusterSpecifyKey(iv, cmfc);

    ++m_cm_clusters;

    if (nclusters[iv] == 1)
    {
      ++m_cm_clusters_size1;
    }
    else if (nclusters[iv] == 2)
    {
      ++m_cm_clusters_size2;
    }
  }

  // read back the clusters and make some histograms

  auto clusrange = _corrected_CMcluster_map->getClusters();
  for (auto cmitr = clusrange.first;
       cmitr != clusrange.second;
       ++cmitr)
  {
    auto cmkey = cmitr->first;
    auto cmclus = cmitr->second;

    if (Verbosity() > 0)
    {
      std::cout << "found CM cluster " << cmkey << " with adc " << cmclus->getAdc()
                << " x " << cmclus->getX() << " y " << cmclus->getY() << " z " << cmclus->getZ()
                << " nclusters " << cmclus->getNclusters()
                << std::endl;
    }

    if (_histos)
    {
      henergy->Fill(cmclus->getAdc());
      hxy->Fill(cmclus->getX(), cmclus->getY());
      hz->Fill(cmclus->getZ());
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTpcCentralMembraneClusterizer::End(PHCompositeNode * /*topNode*/)
{
  if (_histos)
  {
    m_histogramfile->cd();

    henergy->Write();
    hxy->Write();
    hz->Write();

    hz_pos->Write();
    hz_neg->Write();

    hClustE[0]->Write();
    hClustE[1]->Write();
    hClustE[2]->Write();
    hDist->Write();
    hDistRow->Write();
    hDist2->Write();
    hDistRowAdj->Write();
    hDist2Adj->Write();

    for (int i = 0; i < 47; i++)
    {
      hphi_reco_pair_pos[i]->Write();
      hphi_reco_pair_neg[i]->Write();
    }

    m_histogramfile->Close();
  }

  // print statistics
  if (Verbosity())
  {
    std::cout
        << "PHTpcCentralMembraneClusterizer::End -"
        << " cluster statistics total: " << m_total_clusters
        << " accepted: " << m_accepted_clusters << " fraction: "
        << 100. * m_accepted_clusters / m_total_clusters << "%"
        << std::endl;

    std::cout
        << "PHTpcCentralMembraneClusterizer::End -"
        << " cm clusters: " << m_cm_clusters
        << std::endl;

    std::cout
        << "PHTpcCentralMembraneClusterizer::End -"
        << " cm clusters size 1: " << m_cm_clusters_size1
        << std::endl;

    std::cout
        << "PHTpcCentralMembraneClusterizer::End -"
        << " cm clusters size 2: " << m_cm_clusters_size2
        << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTpcCentralMembraneClusterizer::GetNodes(PHCompositeNode *topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!_geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // tpc distortion correction
  _dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainer");
  if (_dcc)
  {
    std::cout << "PHTpcCentralMembraneMatcher:   found TPC distortion correction container" << std::endl;
  }

  _corrected_CMcluster_map = findNode::getClass<CMFlashClusterContainer>(topNode, "CORRECTED_CM_CLUSTER");
  if (!_corrected_CMcluster_map)
  {
    std::cout << "Creating node CORRECTED_CM_CLUSTER" << std::endl;
    PHNodeIterator iter(topNode);

    // Looking for the DST node
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    _corrected_CMcluster_map = new CMFlashClusterContainerv1;
    PHIODataNode<PHObject> *TrkrClusterContainerNode = new PHIODataNode<PHObject>(_corrected_CMcluster_map, "CORRECTED_CM_CLUSTER", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
