// $Id: Prototype3DSTReader.h,v 1.7 2015/02/27 23:42:23 jinhuang Exp $

/*!
 * \file Prototype3DSTReader.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.7 $
 * \date $Date: 2015/02/27 23:42:23 $
 */

#ifndef Prototype3DSTReader_H_
#define Prototype3DSTReader_H_

#include "RawTower_Prototype3.h"
#include "RawTower_Temperature.h"

#include <fun4all/SubsysReco.h>

#include <TClonesArray.h>

#include <iostream>
#include <memory>
#include <string>
#include <vector>

class TTree;


/*!
 * \brief Prototype3DSTReader save information from DST to an evaluator, which could include hit. particle, vertex, towers and jet (to be activated)
 */
class Prototype3DSTReader : public SubsysReco
{
public:
  Prototype3DSTReader(const std::string &filename);
  virtual
  ~Prototype3DSTReader();

  //! full initialization
  int
  Init(PHCompositeNode *);

  //! event processing method
  int
  process_event(PHCompositeNode *);

  //! end of run method
  int
  End(PHCompositeNode *);

  void
  AddTower(const std::string &name)
  {
    _tower_postfix.push_back(name);
  }

  void
  AddTowerTemperature(const std::string &name)
  {
    _towertemp_postfix.push_back(name);
  }

  void
  AddRunInfo(const std::string &name)
  {
    _runinfo_list.push_back(name);
  }
  void
  AddEventInfo(const std::string &name)
  {
    _eventinfo_list.push_back(name);
  }

  //! zero suppression for all calorimeters
  double
  get_tower_zero_sup()
  {
    return _tower_zero_sup;
  }

  //! zero suppression for all calorimeters
  void
  set_tower_zero_sup(double b)
  {
    _tower_zero_sup = b;
  }

protected:

//  std::vector<std::string> _node_postfix;
  std::vector<std::string> _tower_postfix;
  //! tower temperature
  std::vector<std::string> _towertemp_postfix;
//  std::vector<std::string> _jet_postfix;
//  std::vector<std::string> _node_name;
  std::vector<std::string> _runinfo_list;
  std::vector<std::string> _eventinfo_list;

  int nblocks;

#ifndef __CINT__

  typedef std::shared_ptr<TClonesArray> arr_ptr;

  struct record
  {
    unsigned int _cnt;
    std::string _name;
    arr_ptr _arr;
    TClonesArray * _arr_ptr;
    double _dvalue;

    enum enu_type
    {
      typ_hit, typ_part, typ_vertex, typ_tower, typ_jets, typ_runinfo, typ_eventinfo, typ_towertemp
    };
    enu_type _type;
  };
  typedef std::vector<record> records_t;
  records_t _records;

  typedef RawTower_Prototype3 RawTower_type;

  typedef RawTower_Temperature RawTowerT_type;
#endif

  int _event;

  std::string _out_file_name;

//  TFile * _file;
  TTree * _T;

  //! zero suppression for all calorimeters
  double _tower_zero_sup;

  void
  build_tree();
};

#endif /* Prototype3DSTReader_H_ */
