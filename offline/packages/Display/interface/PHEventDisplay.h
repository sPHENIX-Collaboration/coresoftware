/*!
        \file PHEventDisplay.h
        \author Sookhyun Lee
        \brief event display interface,
               parameters/switches set, detector modules called, control display.
        \version $Revision: 1.1 $
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
class TGLAutoRotator;
class mPHEveModuleBase;

class PHEventDisplay : public SubsysReco
{
public:
  PHEventDisplay(int w, 
		  int h,
		  bool use_fieldmap,
		  const std::string& mapname,
		  const std::string& geoname,
		  float jet_pt_threshold,
		  float jet_e_scale,
		  float calo_e_threshold,
		  bool  is_svtx_on,
		  bool  is_calo_on,
		  bool  is_jet_on);

  ~PHEventDisplay();

  /// Module initialization
  int Init(PHCompositeNode *topNode);
  /// Run initialization
  int InitRun(PHCompositeNode *topNode);
  /// Event processing
  int process_event(PHCompositeNode *topNode);
  /// Event reset
  int ResetEvent(PHCompositeNode *topNode);
  /// End of process
  int End(PHCompositeNode *topNode);

  /// Return a pointer to the underlying TEveManager
  TEveManager* get_eve_instance();
  /// Threaded access to Fun4All server
  void run_evt_in_thread();

  void start_rotation();
  void go_fullscreen();

private:
  void reco_thread();
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

  TEveManager* _eve_manager;
  TGLAutoRotator* _rot;

#ifndef __CINT__
  boost::shared_ptr<boost::thread> _status_thread;
#endif
  bool _is_svtx_on;
  bool _is_calo_on;
  bool _is_jet_on;
};

#endif // __PHONLINEDISPLAY_H__
