// make a jet probe and put it onto the node tree
// It will always be a single Jet that the JetReco can get
#ifndef G4_JETS_JETPROBEMAKER__H
#define G4_JETS_JETPROBEMAKER__H

#include <fun4all/SubsysReco.h>
#include <math.h>

class Jet;
class JetContainer;
class TRandom3;
class PHCompositeNode;

class JetProbeMaker : public SubsysReco {
  public: 
    JetProbeMaker(const std::string &name="JetProbeMaker", const int _randseed=0); 
    ~JetProbeMaker() override;

    int process_event(PHCompositeNode* /*topNode*/) override;
    int InitRun(PHCompositeNode *topNode) override;

    void set_abs_eta(float val) { _eta_min = -val; _eta_max=val; }
    void set_eta_min(float val) { _eta_min = val; }
    void set_eta_max(float val) { _eta_max = val; }

    void set_phi_min(float val) { _phi_min=val; }
    void set_phi_max(float val) { _phi_max=val; }

    void set_pt     (float val) { _const_pt = true; _pt = val; _pt_min = val; }
    void set_pt_min (float val) { _pt_min = val; if (_pt != _pt_min) _const_pt = false; }
    void set_pt_max (float val) { _pt = val; if (_pt != _pt_min) _const_pt = false; }

  private:

    TRandom3* _rand;
    float _eta_min  = -0.7;
    float _eta_max  = 0.7;
    float _phi_min  = -M_PI;
    float _phi_max  = M_PI;
    float _pt_min   = 30.;
    float _pt       = 30.;
    bool  _const_pt = true;
    JetContainer* _jets = nullptr;
};

#endif
