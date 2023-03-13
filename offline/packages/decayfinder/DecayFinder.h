#ifndef DECAYFINDER_DECAYFINDER_H
#define DECAYFINDER_DECAYFINDER_H

//sPHENIX stuff
#include <fun4all/SubsysReco.h>

#include <g4main/PHG4Particle.h>

#include <cstddef>              // for NULL
#include <string>
#include <vector>
#include <utility>               // for pair

class DecayFinderContainer_v1;
class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHHepMCGenEvent;
class PHHepMCGenEventMap;
namespace HepMC { class GenParticle; }

class DecayFinder : public SubsysReco
{
 public:
  typedef std::vector<std::pair<std::pair<int, int>, int>> Decay;

  DecayFinder();

  explicit DecayFinder(const std::string &name);

  virtual ~DecayFinder() {}

  int Init(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode);

  int End(PHCompositeNode *topNode);

  int parseDecayDescriptor();

  bool findDecay(PHCompositeNode *topNode);

  bool findParticle(const std::string particle);

  void searchHepMCRecord(HepMC::GenParticle* particle, std::vector<int> decayProducts, 
                         bool &breakLoop, bool &hasPhoton, bool &hasPi0, bool &failedPT, bool &failedETA, 
                         std::vector<int> &correctDecayProducts);

  void searchGeant4Record(int barcode, int pid, std::vector<int> decayProducts,
          		 bool &breakLoop, bool &hasPhoton, bool &hasPi0, bool &failedPT, bool &failedETA, 
			 std::vector<int> &correctDecayProducts);

  bool checkIfCorrectHepMCParticle(HepMC::GenParticle *particle, bool &trackFailedPT, bool &trackFailedETA);

  bool checkIfCorrectGeant4Particle(PHG4Particle *particle, bool& hasPhoton, bool& hasPi0, bool& trackFailedPT, bool& trackFailedETA);

  bool compareDecays(std::vector<int> required, std::vector<int> actual);

  int deleteElement(int arr[], int n, int x);

  void multiplyVectorByScalarAndSort(std::vector<int> &v, int k);

  int get_pdgcode(const std::string name);

  int get_charge(const std::string name);

  bool isInRange(float min, float value, float max);

  int createDecayNode(PHCompositeNode *topNode);

  void fillDecayNode(PHCompositeNode* topNode, Decay decay);

  void printInfo();

  void printNode(PHCompositeNode *topNode);

  //User configuration
  /**
   * Use this function to define the decay you want to find in the HepMC record
   * @param[in] decayDescriptor the description of the decay chain, this is a string
   * You define the decay with these rules:
   * @brief You define a particle decaying with "->", the mother on the left, the decay products on the right
   * @brief Set the charge of final state tracks with "^", the particle name on the left and the charge on the right. 
   *        Accepted charges are +, - and 0
   * @brief Use the same rules as above for any intermediatee decays but contain the entire decay within curled 
   *        brackets, "{}"
   * @brief If you also want to find the charge conjugate decay, contain the entire decay descriptor within "[]cc" for 
   *        charge-conjugate. The "cc" is NOT case sensitive
   * @brief The particle names you use must be kept in the TDatabasePDG class from root 
   *        (https://root.cern.ch/doc/master/classTDatabasePDG.html). Print this table to see available particles with 
   *        TDatabasePDG::Instance()->Print()
   * @brief An example of a decay would be: "[B+ -> {D0_bar -> kaon^+ pion^-} pion^+]cc"
   * @note There is an internal list of resonances which, if they appear in the record, will be further analysed. For 
   *       example, the f0(980)->pipi decay is too quick to have a flight distance and so we would only see the pion 
   *       pair in the detector. If you are looking for B_s0 -> J/psi pipi then the decay of the f0 will be studied 
   *       for a pipi final state, basically inclusive decays are handled automatically. If you wish to study the f0 
   *       decay, add it to your decay descriptor and it will automatically be removed from the "skip list"
   */
  void setDecayDescriptor(const std::string &decayDescriptor) { m_decayDescriptor = decayDescriptor; }
  /**
   * @param[in] trigger Set to true to allow further processing of events in which your decay appears, if your decay 
   *            does not appear, all further processing of this event is skipped. This defaults to false so every event 
   *            is proccessed in F4A 
   */
  void triggerOnDecay(bool trigger) { m_triggerOnDecay = trigger; }
  /**
   * @param[in] allow Set to true to allow photons to be associated to your decay
   */
  void allowPhotons(bool allow) { m_allowPhotons = allow; }
  /**
   * @param[in] allow Set to true to allow pi zero to be associated to your decay
   */
  void allowPi0(bool allow) { m_allowPi0 = allow; }
  /**
   * @param[in] allow Set to true to save any of your decays that are found back to the node tree in a DecayFinderContainer
   *           The default name is "decay" and will automatically have "_DecayMap" added to the end
   */
  void saveDST(bool save) { m_save_dst = save; }
  /**
   * @param[in] name Change the default name of the DecayFinderContainer. 
   * @note This name will still have "_DecayMap" added to the end, this cannot be changed
   */
  void setNodeName(const std::string &name) { m_container_name = name; }

  /**
   * @param[in] min The minimum eta threshold for track acceptance
   * @param[in] min The maximum eta threshold for track acceptance
   * @note Set a pseudorapidity threshold range for tracking
   */
  void setEtaRange(float min, float max) { m_eta_low_req = min; m_eta_high_req = max; }

  /**
   * @param[in] pt The minimum pT threshold for track acceptance
   * @note Set a minimum pT threshold for tracking
   */
  void setPTmin(float pt) { m_pt_req = pt; }

 private:
  PHHepMCGenEventMap *m_geneventmap = nullptr;
  PHHepMCGenEvent *m_genevt = nullptr;
  PHG4TruthInfoContainer *m_truthinfo = nullptr;

  double m_eta_high_req = 1.1;
  double m_eta_low_req = -1.1;
  double m_pt_req = 0.2;

  int m_counter = 0;
  int m_nCandFail_pT = 0;
  int m_nCandFail_eta = 0;
  int m_nCandFail_pT_and_eta = 0;
  int m_nCandReconstructable = 0;
  int m_nCandHas_Photon = 0;
  int m_nCandHas_Pi0 = 0;
  int m_nCandHas_Photon_and_Pi0 = 0;
  int m_nCandHas_noPhoton_and_noPi0 = 0;

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

  bool m_save_dst;
  DecayFinderContainer_v1 *m_decayMap = nullptr;
  Decay decayChain;
  std::string m_nodeName;
  std::string m_container_name;
};

#endif  //DECAYFINDER_DECAYFINDER_H
