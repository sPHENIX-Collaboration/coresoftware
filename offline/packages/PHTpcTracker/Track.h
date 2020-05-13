/*!
 *  \file		Track.h
 *  \brief		Data structure and output of the fitting.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef PHGENFIT2_TRACK_H
#define PHGENFIT2_TRACK_H

#include <trackbase/TrkrDefs.h>

#include <TMatrixDSymfwd.h>
#include <TVector3.h>

//STL
#include <map>
#include <memory>
#include <vector>

namespace PHGenFit { class Measurement; }

namespace genfit
{
  class AbsTrackRep;
  class MeasuredStateOnPlane;
  class Track;
}  // namespace genfit

namespace PHGenFit2
{
  class Track
  {
   public:
    //! Default ctor
    Track(genfit::AbsTrackRep* rep, TVector3 seed_pos, TVector3 seed_mom, TMatrixDSym seed_cov, const int v = 0);

    //! Copy constructor
    Track(const PHGenFit2::Track& t);

    //! Default dtor
    ~Track();

    //! Add measurement
    int addMeasurement(PHGenFit::Measurement* measurement, int id = -1);
    int addMeasurements(std::vector<PHGenFit::Measurement*>& measurements, int id = -1);

    int deleteLastMeasurement();

    //!
    int updateOneMeasurementKalman(
        const std::vector<PHGenFit::Measurement*>& measurements,
        std::map<double, std::shared_ptr<PHGenFit2::Track> >& incr_chi2s_new_tracks,
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

    bool get_track_info(TVector3& pos, TVector3& mom, double& charge, int& nhits, double& length);

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

    void set_vertex_id(const unsigned int vert_id)
    {
      _vertex_id = vert_id;
    }

    unsigned int get_vertex_id() const
    {
      //std::cout << " Track: returning vertex_id = " << _vertex_id << std::endl;
      return _vertex_id;
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
#if defined(__CINT__) && !defined(__CLING__)
    Track operator=(Track& trk)
    {
    }
#else
    Track operator=(Track& trk) = delete;
#endif

    int verbosity;

    genfit::Track* _track;
    //std::vector<PHGenFit::Measurement*> _measurements;
    std::vector<unsigned int> _clusterIDs;
    std::vector<TrkrDefs::cluskey> _clusterkeys;
    unsigned int _vertex_id;

    //SMART(genfit::Track) _track;
  };
}  // namespace PHGenFit2

#endif
