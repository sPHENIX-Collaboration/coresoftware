// make a jet probe and put it onto the node tree
// It will always be a single Jet that the JetReco can get
#ifndef G4_JETS_JETPROBEMAKER__H
#define G4_JETS_JETPROBEMAKER__H

#include <fun4all/SubsysReco.h>
#include <gsl/gsl_rng.h>
#include <memory>  // for unique_ptr
#include <math.h>

class Jet;
class JetContainer;
class TRandom3;
class PHCompositeNode;

class JetProbeMaker : public SubsysReco {
  public: 
    JetProbeMaker(const std::string &name="JetProbeMaker");
    ~JetProbeMaker() override { };

    int process_event(PHCompositeNode* /*topNode*/) override;
    int InitRun(PHCompositeNode *topNode) override;

    void set_abs_eta(float val) { _eta_min = -val; _eta_max=val; _eta_range = (_eta_max-_eta_min); }
    void set_eta_min(float val) { _eta_min = val; _eta_range = (_eta_max-_eta_min); }
    void set_eta_max(float val) { _eta_max = val; _eta_range = (_eta_max-_eta_min); }

    void set_phi_min(float val) { _phi_min=val; }
    void set_phi_max(float val) { _phi_max=val; }

    void set_pt     (float val) { _pt_max = val; _pt_min = val; _pt_range = _pt_max-_pt_min; }
    void set_pt_min (float val) { _pt_min = val; _pt_range = _pt_max-_pt_min; }
    void set_pt_max (float val) { _pt_max = val; _pt_range = _pt_max-_pt_min; }

  private:

    class Deleter
    {
      public:
        //! delection operation
        void operator()(gsl_rng *rng) const { gsl_rng_free(rng); }
    };
    std::unique_ptr<gsl_rng, Deleter> m_rng;

    float _eta_min  = -0.7;
    float _eta_max  = 0.7;
    float _eta_range = 1.4;
    float _phi_min  = -M_PI;
    float _phi_max  = M_PI;
    float _pt_min   = 30.;
    float _pt_max   = 30.;
    float _pt_range = 0.;
    JetContainer* _jets = nullptr;
};

#endif
