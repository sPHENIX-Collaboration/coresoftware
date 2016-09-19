/*!
        \file PHEventDisplay.h
        \author Sookhyun Lee
        \brief event display interface,
              set parameters/switches, call detector modules, control display.
        \version $Revision: 1.2 $
        \date    $Date: 07/26/2016
*/

#ifndef __PHEVENTDISPLAY_H__
#define __PHEVENTDISPLAY_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include <vector>

#ifndef __CINT__
#include<boost/shared_ptr.hpp>
#include<boost/thread/thread.hpp>
#else
class shared_ptr;
#endif

#include <pthread.h>

class PHEveDisplay;
class PHCompositeNode;
class TEveManager;
class TGLViewer;
class TGLAutoRotator;
class mPHEveModuleBase;

class PHEventDisplay : public SubsysReco
{
public:
  PHEventDisplay(int w, 
		  int h,
		  bool _use_fieldmap,
		  bool _use_geofile,
		  const std::string& _mapname,
		  const std::string& _geoname);

  ~PHEventDisplay();

  /// Module initialization
  int Init(PHCompositeNode *topNode);
  /// Run initialization
  int InitRun(PHCompositeNode *topNode);
  /// Event processing
  int process_event(PHCompositeNode *topNode);
  /// End of process
  int End(PHCompositeNode *topNode);

  /// Threaded access to Fun4All server
  void run_evt_in_thread();
  void start_rotation();
  void go_fullscreen();
  void set_jet_pt_threshold(float pt){jet_pt_threshold = pt;}
  void set_jet_e_scale(float e_scale){jet_e_scale = e_scale;}
  void set_calo_e_threshold(float e){calo_e_threshold = e;}
  void set_svtx_on(bool svtx_on) {is_svtx_on = svtx_on;}
  void set_cemc_on(bool cemc_on) {is_cemc_on = cemc_on;}
  void set_hcalin_on(bool hcalin_on){is_hcalin_on = hcalin_on;}
  void set_hcalout_on(bool hcalout_on){is_hcalout_on = hcalout_on;}
  void set_jet_on(bool jet_on) {is_jet_on = jet_on;}
  void set_truth_on(bool truth_on) {is_truth_on = truth_on;}
  void set_verbosity(int verb) {verbosity = verb;}


private:
  void reco_thread();
  void draw_default();  
  void update_scene();

#ifndef __CINT__
  typedef boost::shared_ptr<mPHEveModuleBase> pBase;

  template<typename T> boost::shared_ptr<T>
    register_module()
    {
	pBase ptr = pBase(new T(_PHEveDisplay)); // Store in a base type shared_ptr
	_modules.push_back(ptr);
	return boost::dynamic_pointer_cast<T>(ptr); // Cast back to derived for return
    }
#endif

  bool _pending_update;

  //! Reconstruction mutex
#ifndef __CINT__
  std::vector<boost::shared_ptr<mPHEveModuleBase> > _modules;
  pthread_mutex_t _mutex;
  boost::shared_ptr<boost::thread> _update_thread;
  boost::shared_ptr<PHEveDisplay> _PHEveDisplay;
#endif

  TGLAutoRotator* _rot;

#ifndef __CINT__
  boost::shared_ptr<boost::thread> _status_thread;
#endif
  float jet_pt_threshold;
  float jet_e_scale;
  float calo_e_threshold;
  bool is_svtx_on;
  bool is_cemc_on;
  bool is_hcalin_on;
  bool is_hcalout_on;
  bool is_jet_on;
  bool is_truth_on;
  bool use_fieldmap;
  bool use_geofile;
  int width;
  int height;
  std::string mapname;
  std::string geoname;
  int nevent;
  int verbosity;
};

#endif // __PHEVENTDISPLAY_H__
