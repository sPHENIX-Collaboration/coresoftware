#ifndef TPC_TPCCLUSTERZCROSSINGCORRECTION_H
#define TPC_TPCCLUSTERZCROSSINGCORRECTION_H

/*!
 * \file TpcClusterZCrossingCorrection.h
 * \brief applies correction to TPC cluster Z for bunch crossing time offset
 * \author Tony Frawley, March 2022
 */

class TpcClusterZCrossingCorrection
{
  public:

  //!@name accessors
  //@{

  //! drift velocity (cm/ns)
  static float get_vdrift() { return _vdrift; }

  //! time between crossing (ns)
  static float get_time_between_crossings() { return _time_between_crossings; }

  //! apply correction on a given z
  static float correctZ(float zinit, unsigned int side, short int crossing);
  //@}

  //!@name modifiers
  //@{
  //! drift velocity (cm/ns)
  static void set_vdrift( float value ) { _vdrift = value; }

  //! time between crossing (ns)
  static void set_time_between_crossings( float value ) { _time_between_crossings = value; }
  //@}

  // TODO: move to private
  //!@name parameters
  //@{
  //! drift velocity (cm/ns)
  static float _vdrift;

  private:

  //! time between crossing (ns)
  static float _time_between_crossings;

  //@}

};

#endif
