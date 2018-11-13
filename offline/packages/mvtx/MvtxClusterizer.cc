/**
 * @file mvtx/MvtxClusterizer.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of MvtxClusterizer
 */
#include "MvtxClusterizer.h"

#include "MvtxDefs.h"
#include "MvtxHit.h"

#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv1.h>


#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>

#include <TMatrixF.h>
#include <TVector3.h>

#define BOOST_NO_HASH // Our version of boost.graph is incompatible with GCC-4.3 w/o this flag
#include <boost/bind.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
using namespace boost;

#include <iostream>
#include <stdexcept>
#include <cmath>

using namespace std;

static const float twopi = 2.0 * M_PI;

bool MvtxClusterizer::are_adjacent(const pixel lhs,
                                   const pixel rhs )
{
  
  if ( GetZClustering() )
  {
    // column is first, row is second
    if ( fabs(lhs.first - rhs.first) <= 1 )
    {
      if ( fabs(lhs.second - rhs.second) <= 1 )
      {
        return true;
      }
    }
  }
  else
  {
    if ( fabs(lhs.first - rhs.first) == 0 )
    {
      if ( fabs(lhs.second - rhs.second) <= 1 )
      {
        return true;
      }
    }
  }
  
  return false;
}


MvtxClusterizer::MvtxClusterizer(const string &name) 
 : SubsysReco(name)
 , m_hits(nullptr)
 , m_clusterlist(nullptr)
 , m_makeZClustering(true)
{
  
}

int MvtxClusterizer::InitRun(PHCompositeNode* topNode)
{
  
  //-----------------
  // Add Cluster Node
  //-----------------
  
  PHNodeIterator iter(topNode);
  
  // Looking for the DST node
  PHCompositeNode *dstNode
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) 
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);
  
  // Create the SVX node if required
  PHCompositeNode* svxNode
    = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode", "TRKR"));
  if (!svxNode) 
  {
    svxNode = new PHCompositeNode("TRKR");
    dstNode->addNode(svxNode);
  }
  
  // Create the Cluster node if required
  TrkrClusterContainer *trkrclusters
    = findNode::getClass<TrkrClusterContainer>(dstNode, "TrkrClusterContainer");
  if (!trkrclusters) 
  {
    trkrclusters = new TrkrClusterContainer();
    PHIODataNode<PHObject> *TrkrClusterContainerNode =
      new PHIODataNode<PHObject>(trkrclusters, "TrkrClusterContainer", "PHObject");
    svxNode->addNode(TrkrClusterContainerNode);
  }
  
  //----------------
  // Report Settings
  //----------------
  
  if (Verbosity() > 0)
  {
    cout << "====================== MvtxClusterizer::InitRun() =====================" << endl;
    cout << " Z-dimension Clustering = " << boolalpha << m_makeZClustering << noboolalpha << endl;
    cout << "===========================================================================" << endl;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int MvtxClusterizer::process_event(PHCompositeNode *topNode) {
  
  // get node containing the digitized hits
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TrkrHitSetContainer");
  if (!m_hits)
  {
    cout << PHWHERE << "ERROR: Can't find node TrkrHitSetContainer" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  // get node for clusters
  m_clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TrkrClusterContainer");
  if (!m_clusterlist)
  {
    cout << PHWHERE << " ERROR: Can't find TrkrClusterContainer." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_clusterlist->Reset();
  
  // run clustering
  ClusterMvtx(topNode);
  PrintClusters(topNode);
  
  // done
  return Fun4AllReturnCodes::EVENT_OK;
}

void MvtxClusterizer::ClusterMvtx(PHCompositeNode *topNode) {
  
  if (Verbosity() > 0)
    cout << "Entering MvtxClusterizer::ClusterMvtx " << endl;
  
  //-----------
  // Clustering
  //-----------

  // loop over each MvtxHitSet object (chip)
  TrkrHitSetContainer::ConstRange hitsetrange =
    m_hits->getHitSets(TrkrDefs::TrkrId::mvtxId);
  for ( TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
        hitsetitr != hitsetrange.second;
        ++hitsetitr )
  {
    
    TrkrHitSet* hitset = hitsetitr->second;

    if (Verbosity()>2)
      hitset->identify();
    

    // fill a vector of hits to make things easier
    // D. McGlinchey - this probably isn't necessary. Use iterator directly?
    std::vector<pixel> hitvec;
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    for ( TrkrHitSet::ConstIterator hitr = hitrangei.first;
          hitr != hitrangei.second;
          ++hitr)
    {
      // get the column and row indeces from the hitkey 
      // (we don't actually care about the hit itself for the mvtx)
      unsigned int col = MvtxDefs::getCol(hitr->first);
      unsigned int row = MvtxDefs::getRow(hitr->first);
      hitvec.push_back(make_pair(col, row));
    }
    if ( Verbosity()>2 )
      cout << "hitvec.size(): " << hitvec.size() << endl;
    
    // do the clustering
    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    Graph G;

    // loop over hits in this chip
    for ( unsigned int i = 0; i < hitvec.size(); i++)
    {
      for ( unsigned int j = 0; j < hitvec.size(); j++)
      {
        if ( are_adjacent(hitvec[i], hitvec[j]) )
          add_edge(i, j, G);
      }
    }
    
    
    // Find the connections between the vertices of the graph (vertices are the rawhits,
    // connections are made when they are adjacent to one another)
    vector<int> component(num_vertices(G));
    
    // this is the actual clustering, performed by boost
    connected_components(G, &component[0]);
    
    // Loop over the components(hits) compiling a list of the
    // unique connected groups (ie. clusters).
    set<int> cluster_ids; // unique components
    multimap<int, pixel> clusters;
    for (unsigned int i = 0; i < component.size(); i++)
    {
      cluster_ids.insert( component[i] );
      clusters.insert( make_pair(component[i], hitvec[i]) );
    }

    // loop over the componenets and make clusters
    for (set<int>::iterator clusiter = cluster_ids.begin();
         clusiter != cluster_ids.end();
         clusiter++ )
    {
      
      int clusid = *clusiter;
      pair< multimap<int, pixel>::iterator,
            multimap<int, pixel>::iterator> clusrange = clusters.equal_range(clusid);

      multimap<int, pixel>::iterator mapiter = clusrange.first;

      if (Verbosity() > 2)
        cout << "Filling cluster id " << clusid << endl;

      // make cluster
      TrkrDefs::cluskey ckey = MvtxDefs::genClusKey(hitset->getHitSetKey(),clusid);
      TrkrClusterv1* clus = static_cast<TrkrClusterv1*>((m_clusterlist->findOrAddCluster(ckey))->second);

      // determine the size of the cluster in phi and z
      set<int> phibins;
      set<int> zbins;

      // determine the cluster position...
      double xsum = 0.0;
      double ysum = 0.0;
      double zsum = 0.0;
      unsigned nhits = 0;


      for (mapiter = clusrange.first; mapiter != clusrange.second; mapiter++ )
      {

        // size
        zbins.insert((mapiter->second).first);
        phibins.insert((mapiter->second).second);

        // find the center of the pixel in local coords
        xsum += (mapiter->second).second;
        zsum += (mapiter->second).first + 0.5;

        ++nhits;
      } //mapitr

      double thickness = 50.e-4/28e-4;
      double phisize = phibins.size();
      double zsize = zbins.size();

      double clusx = NAN;
      double clusy = NAN;
      double clusz = NAN;

      clusx = xsum / nhits;
      clusy = ysum / nhits;
      clusz = zsum / nhits;

      clus->setPosition(0, clusx);
      clus->setPosition(1, clusy);
      clus->setPosition(2, clusz);
      clus->setLocal();
      
      clus->setAdc(nhits);

      double invsqrt12 = 1.0 / sqrt(12.0);

      TMatrixF DIM(3, 3);
      DIM[0][0] = pow(0.5 * phisize, 2);

      DIM[0][1] = 0.0;
      DIM[0][2] = 0.0;
      DIM[1][0] = 0.0;
      DIM[1][1] = pow(0.5 * thickness, 2);
      DIM[1][2] = 0.0;
      DIM[2][0] = 0.0;
      DIM[2][1] = 0.0;
      DIM[2][2] = pow(0.5 * zsize, 2);

      TMatrixF ERR(3, 3);
      ERR[0][0] = pow(0.5 * phisize * invsqrt12, 2);
      ERR[0][1] = 0.0;
      ERR[0][2] = 0.0;
      ERR[1][0] = 0.0;
      ERR[1][1] = pow(0.5 * thickness * invsqrt12, 2);
      ERR[1][2] = 0.0;
      ERR[2][0] = 0.0;
      ERR[2][1] = 0.0;
      ERR[2][2] = pow(0.5 * zsize * invsqrt12, 2);


      clus->setSize( 0 , 0 , DIM[0][0] );
      clus->setSize( 0 , 1 , DIM[0][1] );
      clus->setSize( 0 , 2 , DIM[0][2] );
      clus->setSize( 1 , 0 , DIM[1][0] );
      clus->setSize( 1 , 1 , DIM[1][1] );
      clus->setSize( 1 , 2 , DIM[1][2] );
      clus->setSize( 2 , 0 , DIM[2][0] );
      clus->setSize( 2 , 1 , DIM[2][1] );
      clus->setSize( 2 , 2 , DIM[2][2] );

      clus->setError( 0 , 0 , ERR[0][0] );
      clus->setError( 0 , 1 , ERR[0][1] );
      clus->setError( 0 , 2 , ERR[0][2] );
      clus->setError( 1 , 0 , ERR[1][0] );
      clus->setError( 1 , 1 , ERR[1][1] );
      clus->setError( 1 , 2 , ERR[1][2] );
      clus->setError( 2 , 0 , ERR[2][0] );
      clus->setError( 2 , 1 , ERR[2][1] );
      clus->setError( 2 , 2 , ERR[2][2] );


      if ( Verbosity() > 2 )
        clus->identify();

    } // clusitr
  } // hitsetitr

  return;
}


void MvtxClusterizer::PrintClusters(PHCompositeNode * topNode)
{

  if (Verbosity() >= 1) {

    TrkrClusterContainer *clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TrkrClusterContainer");
    if (!clusterlist) return;

    cout << "================= MvtxClusterizer::process_event() ====================" << endl;

    cout << " Found and recorded the following " << clusterlist->size() << " clusters: " << endl;

    clusterlist->identify();

    cout << "===========================================================================" << endl;
  }

  return;
}
