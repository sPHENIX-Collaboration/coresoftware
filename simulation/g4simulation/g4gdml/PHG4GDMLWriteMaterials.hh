//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: PHG4GDMLWriteMaterials.hh 68053 2013-03-13 14:39:51Z gcosmo $
//
//
// class PHG4GDMLWriteMaterials
//
// Class description:
//
// GDML class for writing materials definitions.

// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#ifndef _PHG4GDMLWRITEMATERIALS_INCLUDED_
#define _PHG4GDMLWRITEMATERIALS_INCLUDED_

#include <Geant4/G4Types.hh>
#include <vector>
#include <Geant4/G4Version.hh>

#include "PHG4GDMLWriteDefine.hh"

class G4Isotope;
class G4Element;
class G4Material;
class G4PhysicsFreeVector;

#if G4VERSION_NUMBER >= 1100
class G4MaterialPropertiesTable;
#else
class G4PhysicsOrderedFreeVector;
#endif

class PHG4GDMLWriteMaterials : public PHG4GDMLWriteDefine
{

 public:

   void AddIsotope(const G4Isotope* const);
   void AddElement(const G4Element* const);
   void AddMaterial(const G4Material* const);

   virtual void MaterialsWrite(xercesc::DOMElement*);

 protected:

   PHG4GDMLWriteMaterials();
   virtual ~PHG4GDMLWriteMaterials();

   void AtomWrite(xercesc::DOMElement*,const G4double&);
   void DWrite(xercesc::DOMElement*,const G4double&);
   void PWrite(xercesc::DOMElement*,const G4double&);
   void TWrite(xercesc::DOMElement*,const G4double&);
   void MEEWrite(xercesc::DOMElement*,const G4double&);
   void IsotopeWrite(const G4Isotope* const);
   void ElementWrite(const G4Element* const);
   void MaterialWrite(const G4Material* const);


#if G4VERSION_NUMBER >= 1100

   void PropertyWrite(xercesc::DOMElement*, const G4Material* const);
   void PropertyVectorWrite(const G4String&,
                            const G4PhysicsFreeVector* const);
   void PropertyConstWrite(const G4String&, const G4double,
                           const G4MaterialPropertiesTable*);
#else

   void PropertyWrite(xercesc::DOMElement*, const G4Material* const);
   void PropertyVectorWrite(const G4String&,
                            const G4PhysicsOrderedFreeVector* const);

#endif


 protected:

   std::vector<const G4Isotope*> isotopeList;
   std::vector<const G4Element*> elementList;
   std::vector<const G4Material*> materialList;
   std::vector<const G4PhysicsFreeVector*> propertyList;
   xercesc::DOMElement* materialsElement;
};

#endif
