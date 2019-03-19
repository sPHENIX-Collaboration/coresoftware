/*!
 *  \file		Track.h
 *  \brief		Data structure and output of the fitting.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef PHGENFIT_TRACK_H
#define PHGENFIT_TRACK_K

#include <GenFit/Track.h>

#include <trackbase/TrkrDefs.h>

//STL
#include <map>
#include <memory>
#include <vector>

//GenFit

namespace genfit
{
class AbsTrackRep;
class StateOnPlane;
class Track;

}  // namespace genfit

namespace PHGenFit
{
class Measurement;

class Track
{
 public:
  //! Default ctor
  Track(genfit::AbsTrackRep* rep, TVector3 seed_pos, TVector3 seed_mom, TMatrixDSym seed_cov, const int v = 0);

  //! Copy constructor
  Track(const PHGenFit::Track& t);

  //! Default dtor
  ~Track();

  //! Add measurement
  int addMeasurement(PHGenFit::Measurement* measurement);
  int addMeasurements(std::vector<PHGenFit::Measurement*>& measurements);

  int deleteLastMeasurement();

  //!
  int updateOneMeasurementKalman(
      const std::vector<PHGenFit::Measurement*>& measurements,
      std::map<double, std::shared_ptr<PHGenFit::Track> >& incr_chi2s_new_tracks,
      const int base_tp_idx = -1,
      const int direction = 1,
      const float blowup_factor = 1.,
      const bool use_fitted_state = false) const;

  /*!
	 * track_point 0 is the first one, and -1 is the last one
	 */
  double extrapolateToPlane(genfit::MeasuredStateOnPlane& state, TVector3 O, TVector3 n, const int tr_point_id = 0) const;
  //!
  double extrapolateToLine(genfit::MeasuredStateOnPlane& state, TVector3 line_point, TVector3 line_direction, const int tr_point_id = 0) const;
  //!
  double extrapolateToCylinder(genfit::MeasuredStateOnPlane& state, double radius, TVector3 line_point, TVector3 line_direction, const int tr_point_id = 0, const int direction = 1) const;
  //!
  double extrapolateToPoint(genfit::MeasuredStateOnPlane& state, TVector3 P, const int tr_point_id = 0) const;

  //!
  genfit::MeasuredStateOnPlane* extrapolateToPlane(TVector3 O, TVector3 n, const int tr_point_id = 0) const;
  //!
  genfit::MeasuredStateOnPlane* extrapolateToLine(TVector3 line_point, TVector3 line_direction, const int tr_point_id = 0) const;
  //!
  genfit::MeasuredStateOnPlane* extrapolateToCylinder(double radius, TVector3 line_point, TVector3 line_direction, const int tr_point_id = 0, const int direction = 1) const;
  //!
  genfit::MeasuredStateOnPlane* extrapolateToPoint(TVector3 P, const int tr_point_id = 0) const;
  //!
  genfit::Track* getGenFitTrack() { return _track; }

  genfit::Track* getGenFitTrack() const { return _track; }

  double get_chi2() const;

  double get_ndf() const;

  double get_charge() const;

  TVector3 get_mom() const;

  // old tracking
  const std::vector<unsigned int>& get_cluster_IDs() const
  {
    return _clusterIDs;
  }
  void set_cluster_IDs(const std::vector<unsigned int>& clusterIDs)
  {
    _clusterIDs = clusterIDs;
  }

  // new tracking
  const std::vector<TrkrDefs::cluskey>& get_cluster_keys() const
  {
    return _clusterkeys;
  }
  void set_cluster_keys(const std::vector<TrkrDefs::cluskey>& clusterkeys)
  {
    _clusterkeys = clusterkeys;
  }


  int get_verbosity() const
  {
    return verbosity;
  }

  void set_verbosity(int verbosity)
  {
    this->verbosity = verbosity;
  }

  //SMART(genfit::Track) getGenFitTrack() {return _track;}

 private:
#if defined(__CINT__) && ! defined(__CLING__)
  Track operator=(Track &trk) {}
#else
  Track operator=(Track &trk) = delete; 
#endif

  int verbosity;

  genfit::Track* _track;
  //std::vector<PHGenFit::Measurement*> _measurements;
  std::vector<unsigned int> _clusterIDs;
  std::vector<TrkrDefs::cluskey> _clusterkeys;

  //SMART(genfit::Track) _track;
};
}  // namespace PHGenFit

#endif
