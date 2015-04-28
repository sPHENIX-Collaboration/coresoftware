#ifndef __PHG4SVTXADDCONNECTEDCELLS__
#define __PHG4SVTXADDCONNECTEDCELLS__

#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>
#include <vector>
#include <map>

class PHG4CylinderCell;
class PHG4CylinderCellGeomContainer;
class PHG4CylinderCellContainer;

class PHG4SvtxAddConnectedCells : public SubsysReco {

public:

  PHG4SvtxAddConnectedCells(const char * name = "PHG4SvtxAddConnectedCells");
  ~PHG4SvtxAddConnectedCells(){}
  
  //! module initialization
  int Init(PHCompositeNode *topNode){return 0;}
  
  //! run initialization
  int InitRun(PHCompositeNode *topNode);
  
  //! event processing
  int process_event(PHCompositeNode *topNode);
  
  //! end of process
  int End(PHCompositeNode *topNode){return 0;}


  void set_phi_offset(int offset) {
    connected_phi_offset = offset;
    std::cout << " PHG4SvtxAddConnectedCells: phi bins offset for connected cells set to " << connected_phi_offset << std::endl;
  }
  
  void set_ncells_connected(int layer, int connected) {
    if(layer < 19)
      {
	ncells_connected[layer] = connected;
	std::cout << " PHG4SvtxAddConnectedCells: number of connected cells for layer " << layer << " set to " << connected << std::endl;
      }
    else
      std::cout << " ********* Layer number exceeds maximum layer number of 19, doing nothing!" << std::endl;
  }
    
private:

  static bool lessthan(const PHG4CylinderCell*, 
		       const PHG4CylinderCell*);

  int connected_phi_offset;
  int ncells_connected[20];  // allows for 20 layers max

  // node tree storage pointers
  PHG4CylinderCellContainer* _cells;
  PHG4CylinderCellGeomContainer* _geom_container;

  PHTimeServer::timer _timer;
};

#endif
