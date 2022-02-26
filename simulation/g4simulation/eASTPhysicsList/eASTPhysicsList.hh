////////////////////////////////////////////////////////////////////////////////
//                             
//  eASTPhysicsList.hh 
//  Geant4 physics list for Electron Ion Collider detector simulation
//                                                                  
//  History                                                        
//    Jun.21.2018 : original implementation - Dennis H. Wright (SLAC) 
//    May.06.2021 : migration to eAST - Makoto Asai (SLAC)
//                                                                      
////////////////////////////////////////////////////////////////////////////////

#ifndef eASTPhysicsList_h
#define eASTPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class eASTPhysicsListMessenger;
class G4Region;
class G4ProductionCuts;

#include <map>

class eASTPhysicsList: public G4VModularPhysicsList
{
  public:
    eASTPhysicsList(G4int verb=0);
    ~eASTPhysicsList();

    virtual void ConstructParticle();
    virtual void ConstructProcess();
    virtual void SetCuts();

  private:
    void SetupProcesses();

  private:
    G4bool processesAreRegistered = false;
    eASTPhysicsListMessenger* pMessenger;

    //G4bool addHP = false;      // add Neutron_HP
    G4bool addRDM = false;     // add Radioactive Decay Module
    G4bool addOptical = false; // add optical physics

    G4int stepLimit_opt = -1;  // step limiter option (0:charged, 1:neutral, 2:all, 3:e+/e-)
    std::map<G4Region*,G4double> localStepLimits; // map of region name and limit value
    G4double globalCuts[4];    // for e-, e+ gamma, proton
    std::map<G4Region*,G4ProductionCuts*> localCuts; // map of region name and cuts

  public:
    // Neutron_HP needs further implementation
    //void AddHP(G4bool val = true) { addHP = val; }
    //G4bool IfHP() const { return addHP; }
    void AddRDM(G4bool val = true) { addRDM = val; }
    G4bool IfRDM() const { return addRDM; }
    void AddOptical(G4bool val = true) { addOptical = val; }
    G4bool IfOptical() const { return addOptical; }

    void AddStepLimit(G4int val = 0) { stepLimit_opt = val; }
    G4int IfStepLimit() const { return stepLimit_opt; }
    void SetGlobalStepLimit(G4double);
    G4double GetGlobalStepLimit() const;
    G4Region* SetLocalStepLimit(const G4String&,G4double);
    G4double GetLocalStepLimit(const G4String&) const;
    void SetGlobalCuts(G4double);
    G4double GetGlobalCuts() const { return GetGlobalCut(0); }
    void SetGlobalCut(G4int, G4double);
    G4double GetGlobalCut(G4int i) const { return globalCuts[i]; }
    G4Region* SetLocalCuts(const G4String& reg,G4double val)
    {
      G4Region* regPtr = nullptr;
      for(G4int i=0; i<4; i++)
      {
        regPtr = SetLocalCut(reg,i,val);
        if(!regPtr) return regPtr;
      }
      return regPtr;
    }
    G4double GetLocalCuts(const G4String& reg) const { return GetLocalCut(reg,0); }
    G4Region* SetLocalCut(const G4String&,G4int,G4double);
    G4double GetLocalCut(const G4String&,G4int) const;
    
  private:
    G4Region* FindRegion(const G4String&) const;
};

#endif
