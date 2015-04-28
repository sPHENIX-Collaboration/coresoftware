#ifndef PHG4Detector_h
#define PHG4Detector_h

#include <iostream>
#include <string>

class G4UserSteppingAction;
class G4LogicalVolume;
class PHCompositeNode;

//! base class for phenix detector creation
/*! derived classes must implement construct method, which takes the "world" logical volume as argument */
class PHG4Detector
{

  public:

  //! constructor

  PHG4Detector( PHCompositeNode *Node ):
    topNode(Node),
    verbosity(0),
    name("NONAME"),
    overlapcheck(false)
	{}

    PHG4Detector( PHCompositeNode *Node, const std::string &nam ):
      topNode(Node),
      verbosity(0),
      name(nam),
      overlapcheck(false)
	{}

  //! destructor
  virtual ~PHG4Detector( void )
  {}

  //! construct method
  /*!
  construct all logical and physical volumes relevant for given detector and place them
  inside the world logical volume
  */
  virtual void Construct( G4LogicalVolume* world ) = 0;

  virtual void Verbosity(const int v) {verbosity = v;}

  virtual int Verbosity()  const {return verbosity;}

  virtual G4UserSteppingAction* GetSteppingAction() { return 0; }

  virtual std::string GetName() const {return name;}

  virtual void OverlapCheck(const bool chk = true) {overlapcheck = chk;}

  virtual void Print(const std::string &what = "ALL") const 
  {std::cout << name << ": Print method not implemented" << std::endl;}

 protected:
  PHCompositeNode *topNode;
  int verbosity;
  std::string name;
  bool overlapcheck;
};

#endif
