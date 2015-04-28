#ifndef __PHG4VUserSteppingAction_H__
#define __PHG4VUserSteppingAction_H__

#include <G4UserSteppingAction.hh>

#include <PHCompositeNode.h>

class G4Step;

class PHG4UserSteppingAction : public G4UserSteppingAction
{

  public:
  PHG4UserSteppingAction( void ):
    topNode_( 0 )
  {}

  virtual ~PHG4UserSteppingAction()
  {}

  virtual void UserSteppingAction(const G4Step*) = 0;

  //! set top node (from where particle list is retrieved for passing to geant
  virtual void SetTopNode( PHCompositeNode* topNode )
  { topNode_ = topNode; }

  protected:

  //! access top node
  virtual PHCompositeNode* TopNode( void ) const
  { return topNode_; }

  private:

  //! temporary pointer to top node
  PHCompositeNode* topNode_;

};


#endif //__G4PHPHYTHIAREADER_H__
