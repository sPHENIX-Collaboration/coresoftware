// $Id: PHG4PhPyReader.h,v 1.3 2014/08/27 16:56:03 mccumber Exp $                                                                                             

/*!
 * \file PHG4PhPyReader.h
 * \brief PHPythia Input to PHG
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.3 $
 * \date $Date: 2014/08/27 16:56:03 $
 */
#ifndef PHG4PHPYREADER_H_
#define PHG4PHPYREADER_H_

#include <fun4all/SubsysReco.h>

#include <vector>

class PHCompositeNode;
class PHPythiaHeader;
class PHPythiaContainer;
class TMCParticle;

//! PHPythia Input to PHG4
class PHG4PhPyReader : public SubsysReco
{
public:
  PHG4PhPyReader(const std::string &name = "PHG4PhPyReader");
  virtual
  ~PHG4PhPyReader();
  int
  Init(PHCompositeNode* topNode);
  int
  process_event(PHCompositeNode* topNode);
  int
  End(PHCompositeNode* topNode);

  void set_input_node(std::string name) {
    _input_node = name;
  }

  double
  get_event_vx() const
  {
    return _event_vx;
  }

  double
  get_event_vy() const
  {
    return _event_vy;
  }

  double
  get_event_vz() const
  {
    return _event_vz;
  }

  std::vector<double>
  get_vtx_offset() const
  {
    return _vtx_offset;
  }

  std::vector<double>
  get_vtx_sigma() const
  {
    return _vtx_sigma;
  }

  void
  set_vtx_offset(const double x, const double y, const double z)
  {
    _vtx_offset[0] = x;
    _vtx_offset[1] = y;
    _vtx_offset[2] = z;
  }

  void
  set_vtx_sigma(const double x, const double y, const double z)
  {
    _vtx_sigma[0] = x;
    _vtx_sigma[1] = y;
    _vtx_sigma[2] = z;
  }

  //! KS range for alive particle
  /*! Ref: PYTHIA 6.4
   * The actual status code used to distinguish between
different classes of entries is given in the KS column; codes in the range 1-10 correspond
to remaining entries, and those above 10 to those that have fragmented or decayed.
   */
  enum {
    kMIN_ALIVE_KS = 0,
    kMAX_ALIVE_KS = 10
  };

protected:
  std::string _input_node;

  PHPythiaContainer* phpythia;

  //! vertex offset
  std::vector<double> _vtx_offset;

  //! vertex smearing
  std::vector<double> _vtx_sigma;

  //! event - by - event vertex x
  double _event_vx;

  //! event - by - event vertex y
  double _event_vy;

  //! event - by - event vertex z
  double _event_vz;
};

#endif /* PHG4PHPYREADER_H_ */
