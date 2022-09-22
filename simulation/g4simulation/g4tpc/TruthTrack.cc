
#include <TLorentzVector.h>




  TLorentzVector v1;
  v1.SetPxPyPzE(px1, py1, pz1, E1);
  virtual int get_vtx_id() const { return -9999; }

  PHG4TruthInfoContainer::PHG4VtxPoint* GetVtx(const int vtxid);
