// $Id: $                                                                                             

/*!
 * \file PHG4CylinderCellSpacalv1.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHG4CYLINDERCELLSPACALV1_H_
#define PHG4CYLINDERCELLSPACALV1_H_

#include <PHG4CylinderCellv1.h>

/*!
 * \brief PHG4CylinderCell_Spacalv1
 */
class PHG4CylinderCell_Spacalv1 : public PHG4CylinderCellv1
{
public:
  PHG4CylinderCell_Spacalv1();
  virtual
  ~PHG4CylinderCell_Spacalv1();

  //! Group hit into each fiber, so we allow fiber-fiber variation studies off-Geant production.
  int
  get_fiber_ID() const
  {
    return fiber_ID;
  }

  //! Group hit into each fiber, so we allow fiber-fiber variation studies off-Geant production.
  void
  set_fiber_ID(int fiberId)
  {
    fiber_ID = fiberId;
  }


protected:

  //! Group hit into each fiber, so we allow fiber-fiber variation studies off-Geant production.
  int fiber_ID;

ClassDef(PHG4CylinderCell_Spacalv1,1)

};

#endif /* PHG4CYLINDERCELLSPACALV1_H_ */
