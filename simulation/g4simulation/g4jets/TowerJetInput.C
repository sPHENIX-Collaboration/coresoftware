
#include "TowerJetInput.h"

#include "JetInput.h"
#include "Jet.h"
#include "JetV1.h"

#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

// PHENIX Geant4 includes
#include <g4cemc/RawTowerGeom.h>
#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTower.h>

// standard includes
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

TowerJetInput::TowerJetInput(Jet::SRC input)
  : _verbosity(0),
    _input(input) {
}

void TowerJetInput::identify(std::ostream& os) {
  os << "   TowerJetInput: ";
  if      (_input == Jet::CEMC_TOWER)    os << "TOWER_CEMC to Jet::CEMC_TOWER";
  else if (_input == Jet::HCALIN_TOWER)  os << "TOWER_HCALIN to Jet::HCALIN_TOWER";
  else if (_input == Jet::HCALOUT_TOWER) os << "TOWER_HCALOUT to Jet::HCALOUT_TOWER";
  os << endl;
}

std::vector<Jet*> TowerJetInput::get_input(PHCompositeNode *topNode) {
  
  if (_verbosity > 0) cout << "TowerJetInput::process_event -- entered" << endl;

  RawTowerContainer *towers = NULL;
  RawTowerGeom *geom = NULL;
  if (_input == Jet::CEMC_TOWER) {
    towers = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CEMC");
    geom = findNode::getClass<RawTowerGeom>(topNode,"TOWERGEOM_CEMC");
    if (!towers||!geom) {
      return std::vector<Jet*>();
    }
  } else if (_input == Jet::HCALIN_TOWER) {
    towers = findNode::getClass<RawTowerContainer>(topNode,"TOWER_HCALIN");
    geom = findNode::getClass<RawTowerGeom>(topNode,"TOWERGEOM_HCALIN");
    if (!towers||!geom) {
      return std::vector<Jet*>();
    }
  } else if (_input == Jet::HCALOUT_TOWER) {
    towers = findNode::getClass<RawTowerContainer>(topNode,"TOWER_HCALOUT");
    geom = findNode::getClass<RawTowerGeom>(topNode,"TOWERGEOM_HCALOUT");
    if (!towers||!geom) {
      return std::vector<Jet*>();
    }
  } else {
    return std::vector<Jet*>();
  }
  
  std::vector<Jet*> pseudojets;
  RawTowerContainer::ConstRange begin_end = towers->getTowers();
  RawTowerContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter !=  begin_end.second; ++rtiter) {
    RawTower *tower = rtiter->second;

    int bineta = tower->get_bineta();
    int binphi = tower->get_binphi();
    double eta = geom->get_etacenter(bineta);
    double phi = geom->get_phicenter(binphi);

    double pt = tower->get_energy() / cosh(eta);
    double px = pt * cos(phi);
    double py = pt * sin(phi);
    double pz = pt * sinh(eta);

    Jet *jet = new JetV1();
    jet->set_px(px);
    jet->set_py(py);
    jet->set_pz(pz);
    jet->set_e(tower->get_energy());
    jet->insert_comp(_input,(bineta << 16) | binphi);
    pseudojets.push_back(jet);
  }

  if (_verbosity > 0) cout << "TowerJetInput::process_event -- exited" << endl;

  return pseudojets;
}
