#ifndef HFTRIGGER_H
#define HFTRIGGER_H

//sPHENIX stuff
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <kfparticle_sphenix/KFParticle_particleList.h> 
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phool/getClass.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/IteratorRange.h> 
#include <HepMC/SimpleVector.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//class KFParticle_particleList;

class DecayFinder : public SubsysReco
{
 public:
  DecayFinder();

  explicit DecayFinder(const std::string &name);

  virtual ~DecayFinder(){}

  int Init(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode);

  int End(PHCompositeNode *topNode);

  int parseDecayDescriptor();

  bool findDecay(PHCompositeNode *topNode);

  bool findParticle(std::string particle);

  int checkIfCorrectParticle(HepMC::GenParticle* particle, bool& trackFailedPT, bool& trackFailedETA);

  int deleteElement(int arr[], int n, int x);

  void multiplyVectorByScalarAndSort(std::vector<int> &v, int k);

  void printInfo();

  //User configuration
  void setDecayDescriptor(std::string decayDescriptor) { m_decayDescriptor = decayDescriptor; }

  void triggerOnDecay(bool trigger) { m_triggerOnDecay = trigger;}

  void allowPhotons(bool allow) { m_allowPhotons = allow;}

  void allowPi0(bool allow) { m_allowPi0 = allow;}
 
 private:

  PHHepMCGenEventMap *m_geneventmap = NULL;
  PHHepMCGenEvent *m_genevt = NULL;

  int m_counter = 0;
  int m_nCandFail_pT = 0;
  int m_nCandFail_eta = 0;
  int m_nCandFail_pT_and_eta = 0;
  int m_nCandReconstructable = 0;

  bool m_getChargeConjugate = false;

  std::string m_decayDescriptor;
  bool m_triggerOnDecay = false;
  bool m_allowPi0 = false;
  bool m_allowPhotons = false;

  int m_mother_ID = 0;
  std::vector<int> m_intermediates_ID;
  std::vector<int> m_daughters_ID;

  int m_nTracksFromMother = 0;
  std::vector<int> m_nTracksFromIntermediates;

  std::vector<int> m_motherDecayProducts;
};

#endif  //HFTRIGGER_H
