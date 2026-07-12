#ifndef GLOBALVERTEX_VERTEXDEFS_H
#define GLOBALVERTEX_VERTEXDEFS_H

// Namespace is used to make the Global Vertex and Calo Vertex object not depend
//on each other but allows for CaloVtxReco to pass on which algorithm was used to
// generate the z-vertex

// If you'd like to add an algorithm add it here, then to GlobalVertex::VTXTYPE
// and also GlobalVertexv4::get_position()

namespace VertexDefs
{
  enum CALOALGO
    {
      UNDEFINED=0,
      JETSKEW=1,
      AVGZ=2,
      JETMLP=3,
      CNN=4
    };
};
#endif
