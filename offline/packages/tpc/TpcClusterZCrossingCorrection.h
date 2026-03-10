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
  static double get_vdrift() { return _vdrift; }

  //! time between crossing (ns)
  static double get_time_between_crossings() { return _time_between_crossings; }

  //! apply correction on a given z
  static double correctZ(double zinit, unsigned int side, short int crossing);
  //@}

  //!@name modifiers
  //@{
  //! drift velocity (cm/ns)
  static void set_vdrift( double value ) { _vdrift = value; }

  //! time between crossing (ns)
  static void set_time_between_crossings( double value ) { _time_between_crossings = value; }
  //@}

  // TODO: move to private
  //!@name parameters
  //@{
  //! drift velocity (cm/ns)
  static double _vdrift;

  private:

  //! time between crossing (ns)
  static double _time_between_crossings;

  //@}

};

#endif
