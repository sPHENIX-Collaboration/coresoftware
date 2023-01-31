/*
Creates the EPD Geometry object and adds it to the node tree
*/

#include "EpdGeom.h"

#include <string>

#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>


class EpdGeomMaker : public SubsysReco {
public:
    EpdGeomMaker(const std::string &name = "EpdGeomMaker");
    ~EpdGeomMaker() override;

    void CreateNodeTree(PHCompositeNode *topNode);

private:
    EpdGeom *epd_geom;
};