#include "PHG4SvtxBeamSpotReco.h"
#include "SvtxBeamSpot.h"
#include "SvtxVertexMap.h"
#include "SvtxVertex.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

#include <TPrincipal.h>

#include <iostream>

using namespace std;

PHG4SvtxBeamSpotReco::PHG4SvtxBeamSpotReco(const string &name) :
  SubsysReco(name),
  _pca(2),
  _vertexes(NULL),
  _beamspot(NULL),
  _timer(PHTimeServer::get()->insert_new(name)) {
}

int PHG4SvtxBeamSpotReco::InitRun(PHCompositeNode* topNode) {

  // clear for new run
  _pca.Clear();

  // Looking for the PAR node
  PHNodeIterator iter(topNode);
  
  PHCompositeNode *parNode 
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","PAR"));
  if (!parNode) {
    cout << PHWHERE << "PAR Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator pariter(parNode);
  
  // Create the SVX node if required
  PHCompositeNode* svxNode 
    = dynamic_cast<PHCompositeNode*>(pariter.findFirst("PHCompositeNode","SVTX"));
  if (!svxNode) {
    svxNode = new PHCompositeNode("SVTX");
    parNode->addNode(svxNode);
  }
  
  // Create the SvtxBeamSpot node if required
  _beamspot = findNode::getClass<SvtxBeamSpot>(topNode,"SvtxBeamSpot");
  if (!_beamspot) {
    _beamspot = new SvtxBeamSpot();
    PHIODataNode<PHObject> *SvtxBeamSpotNode =
      new PHIODataNode<PHObject>(_beamspot, "SvtxBeamSpot", "PHObject");
    svxNode->addNode(SvtxBeamSpotNode);
  }

  // Pull the reconstructed track information off the node tree...
  _vertexes = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if(!_vertexes) {
    cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  if (Verbosity() > 0) {
    cout << "=================== PHG4SvtxBeamSpotReco::InitRun() =======================" << endl;
    cout << " Storing cumulative beam spot location under PAR/SVTX/SvtxBeamSpot" << endl;
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxBeamSpotReco::process_event(PHCompositeNode *topNode)
{
  _timer.get()->restart();

  if (_vertexes->empty()) return Fun4AllReturnCodes::EVENT_OK;
  
  for (SvtxVertexMap::ConstIter iter = _vertexes->begin();
       iter != _vertexes->end();
       ++iter) {
    const SvtxVertex* vertex = iter->second;
    Double_t data[2] = {vertex->get_x(),vertex->get_y()};
    Double_t* pdata = &data[0];
    _pca.AddRow(pdata);
  }

  // recalculate beam spot x,y
  const TVectorD* MEAN  = _pca.GetMeanValues();
  const TMatrixD* COVAR = _pca.GetCovarianceMatrix();
  
  _beamspot->set_x((*MEAN)[0]);
  _beamspot->set_y((*MEAN)[1]);

  for (unsigned int i = 0; i < 2; ++i) {
    for (unsigned int j = 0; j <= i; ++j) {
      _beamspot->set_error(i,j,(*COVAR)[i][j]);
    }
  }

  if (Verbosity() > 1) _beamspot->identify();
  
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxBeamSpotReco::End(PHCompositeNode* topNode) {

  if (Verbosity() > 0) {
    cout << "=================== PHG4SvtxBeamSpotReco::End() ===========================" << endl;
    _beamspot->identify();
    cout << "===========================================================================" << endl;
  }
  
  _pca.Clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

