#include "DumpRawTowerGeom.h"

#include <phool/PHIODataNode.h>

#include <calobase/RawTowerGeom.h>

#include <string>

using namespace std;

typedef PHIODataNode<RawTowerGeom> MyNode_t;

DumpRawTowerGeom::DumpRawTowerGeom(const string &NodeName): DumpObject(NodeName)
{
  return ;
}

int DumpRawTowerGeom::process_Node(PHNode *myNode)
{
  RawTowerGeom *rawtowergeom = NULL;
  MyNode_t *thisNode = static_cast <MyNode_t *> (myNode);
  if (thisNode)
    {
      rawtowergeom = thisNode->getData();
    }
  if (rawtowergeom)
    {
      //TODO: use RawTowerGeomContainer in this dump instead of RawTowerGeom

//      *fout << "radius: " << rawtowergeom->get_radius() << endl;
//      *fout << "thickness: " << rawtowergeom->get_thickness() << endl;
//      *fout << "phibins: " << rawtowergeom->get_phibins() << endl;
//      *fout << "phistep: " << rawtowergeom->get_phistep() << endl;
//      *fout << "phimin: " << rawtowergeom->get_phimin() << endl;
//      *fout << "etabins: " << rawtowergeom->get_etabins() << endl;
//      *fout << "etastep: " << rawtowergeom->get_etastep() << endl;
//      *fout << "etamin: " << rawtowergeom->get_etamin() << endl;
    }
  return 0;
}

