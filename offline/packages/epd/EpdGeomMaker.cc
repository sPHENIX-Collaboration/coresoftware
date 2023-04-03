#include "EpdGeomMaker.h"

#include "EpdGeom.h"
#include "EpdGeomV1.h"

#include <string>

#include <fun4all/Fun4AllBase.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>

EpdGeomMaker::EpdGeomMaker(const std::string &name) : SubsysReco(name) {
    this->epd_geom = new EpdGeomV1();
}

EpdGeomMaker::~EpdGeomMaker() {
    delete this->epd_geom;
}

void EpdGeomMaker::CreateNodeTree(PHCompositeNode *topNode) {
    // Find the DST node
    PHNodeIterator node_itr(topNode);
    PHCompositeNode *dst_node = dynamic_cast<PHCompositeNode *>(node_itr.findFirst("PHCompositeNode", "DST"));
    if (!dst_node) {
        std::cout << "PHComposite node created: DST" << std::endl;
        dst_node = new PHCompositeNode("DST");
        topNode->addNode(dst_node);
    }
    PHIODataNode<PHObject> *epd_geom_node = new PHIODataNode<PHObject>(this->epd_geom, "EPD_Geometry", "PHObject");
    dst_node->addNode(epd_geom_node);
}