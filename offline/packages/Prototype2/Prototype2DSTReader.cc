// $Id: Prototype2DSTReader.cc,v 1.11 2015/01/06 02:52:07 jinhuang Exp $

/*!
 * \file Prototype2DSTReader.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.11 $
 * \date $Date: 2015/01/06 02:52:07 $
 */

#include "Prototype2DSTReader.h"

#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllReturnCodes.h>
//#include <PHGeometry.h>

#include <phool/getClass.h>
#include <phool/getClass.h>

#include <calobase/RawTowerContainer.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>

#include <TTree.h>

#include <map>
#include <set>
#include <cassert>

#include<sstream>


using namespace std;

Prototype2DSTReader::Prototype2DSTReader(const string &filename) :
    SubsysReco("Prototype2DSTReader"), nblocks(0), _event(0), //
    _out_file_name(filename), /*_file(NULL), */_T(NULL), //
    _tower_zero_sup(-10000000)
{

}

Prototype2DSTReader::~Prototype2DSTReader()
{
  cout << "Prototype2DSTReader::destructor - Clean ups" << endl;

  if (_T)
    {
      _T->ResetBranchAddresses();
    }

  _records.clear();
}

int
Prototype2DSTReader::Init(PHCompositeNode*)
{

  const static int arr_size = 100;

  if (_tower_postfix.size())
    {
      cout
          << "Prototype2DSTReader::Init - zero suppression for calorimeter towers = "
          << _tower_zero_sup << " GeV" << endl;
    }
  for (vector<string>::const_iterator it = _runinfo_list.begin();
      it != _runinfo_list.end(); ++it)
    {
      const string & nodenam = *it;

      record rec;
      rec._cnt = 0;
      rec._name = nodenam;
      rec._arr = NULL;
      rec._arr_ptr = NULL;
      rec._dvalue = 0;
      rec._type = record::typ_runinfo;

      _records.push_back(rec);

      nblocks++;
    }

  for (vector<string>::const_iterator it = _tower_postfix.begin();
      it != _tower_postfix.end(); ++it)
    {
      const char * class_name = RawTower_type::Class()->GetName();

      const string & nodenam = *it;

      string hname = Form("TOWER_%s", nodenam.c_str());
//      _node_name.push_back(hname);
      cout << "Prototype2DSTReader::Init - saving raw tower info from node: "
          << hname << " - " << class_name << endl;

      record rec;
      rec._cnt = 0;
      rec._name = hname;
      rec._arr = make_shared<TClonesArray>(class_name, arr_size);
      rec._arr_ptr = rec._arr.get();
      rec._dvalue = 0;
      rec._type = record::typ_tower;

      _records.push_back(rec);

      nblocks++;
    }

  for (vector<string>::const_iterator it = _towertemp_postfix.begin();
      it != _towertemp_postfix.end(); ++it)
    {
      const string & nodenam = *it;
      string hname = Form("TOWER_TEMPERATURE_%s", nodenam.c_str());

      cout << "Prototype2DSTReader::Init - saving average tower temperature info from node: "
          << hname<< endl;

      record rec;
      rec._cnt = 0;
      rec._name = hname;
      rec._arr = NULL;
      rec._arr_ptr = NULL;
      rec._dvalue = 0;
      rec._type = record::typ_towertemp;

      _records.push_back(rec);

      nblocks++;
    }
  cout << "Prototype2DSTReader::Init - requested " << nblocks << " nodes"
      << endl;

  build_tree();

  return 0;
}

void
Prototype2DSTReader::build_tree()
{
  cout << "Prototype2DSTReader::build_tree - output to " << _out_file_name
      << endl;

  static const int BUFFER_SIZE = 32000;

  // open TFile
  PHTFileServer::get().open(_out_file_name, "RECREATE");

  _T = new TTree("T", "Prototype2DSTReader");

  nblocks = 0;
  for (records_t::iterator it = _records.begin(); it != _records.end(); ++it)
    {
      record & rec = *it;

      cout << "Prototype2DSTReader::build_tree - Add " << rec._name << endl;

      if (rec._type == record::typ_runinfo)
        {

          const string name_cnt = rec._name;
          const string name_cnt_desc = name_cnt + "/D";
          _T->Branch(name_cnt.c_str(), &(rec._dvalue), name_cnt_desc.c_str(),
              BUFFER_SIZE);
        }
      else if (rec._type == record::typ_tower)
        {

          const string name_cnt = "n_" + rec._name;
          const string name_cnt_desc = name_cnt + "/I";
          _T->Branch(name_cnt.c_str(), &(rec._cnt), name_cnt_desc.c_str(),
              BUFFER_SIZE);
          _T->Branch(rec._name.c_str(), &(rec._arr_ptr), BUFFER_SIZE, 99);
        }
      else if (rec._type == record::typ_towertemp)
        {

          const string name_cnt = rec._name + "_AVG";
          const string name_cnt_desc = name_cnt + "/D";
          _T->Branch(name_cnt.c_str(), &(rec._dvalue), name_cnt_desc.c_str(),
              BUFFER_SIZE);
        }

      nblocks++;
    }

  cout << "Prototype2DSTReader::build_tree - added " << nblocks << " nodes"
      << endl;

  _T->SetAutoSave(16000);
}

int
Prototype2DSTReader::process_event(PHCompositeNode* topNode)
{

//  const double significand = _event / TMath::Power(10, (int) (log10(_event)));
//
//  if (fmod(significand, 1.0) == 0 && significand <= 10)
//    cout << "Prototype2DSTReader::process_event - " << _event << endl;
  _event++;

  for (records_t::iterator it = _records.begin(); it != _records.end(); ++it)
    {
      record & rec = *it;

      rec._cnt = 0;

      if (rec._type == record::typ_hit)
        {
          assert(0);
        } //      if (rec._type == record::typ_hit)
      else if (rec._type == record::typ_tower)
        {
          assert(rec._arr.get() == rec._arr_ptr);
          assert(rec._arr.get());
          rec._arr->Clear();

          if (Verbosity() >= 2)
            cout << "Prototype2DSTReader::process_event - processing tower "
                << rec._name << endl;

          RawTowerContainer *hits = findNode::getClass<RawTowerContainer>(
              topNode, rec._name);
          if (!hits)
            {
              if (_event < 2)
                cout
                    << "Prototype2DSTReader::process_event - Error - can not find node "
                    << rec._name << endl;

            }
          else
            {
              RawTowerContainer::ConstRange hit_range = hits->getTowers();

              if (Verbosity() >= 2)
                cout << "Prototype2DSTReader::process_event - processing "
                    << rec._name << " and received " << hits->size()
                    << " tower hits" << endl;

              for (RawTowerContainer::ConstIterator hit_iter = hit_range.first;
                  hit_iter != hit_range.second; hit_iter++)
                {
                  RawTower * hit_raw = hit_iter->second;

                  RawTower_type * hit = dynamic_cast<RawTower_type *>(hit_raw);
//                  RawTower * hit = hit_iter->second;

                  assert(hit);

                  if (hit->get_energy() < _tower_zero_sup)
                    {

                      if (Verbosity() >= 2)
                        cout
                            << "Prototype2DSTReader::process_event - suppress low energy tower hit "
                            << rec._name << " @ ("
//                            << hit->get_thetaMin()
//                            << ", " << hit->get_phiMin()
                            << "), Energy = " << hit->get_energy() << endl;

                      continue;
                    }

                  new ((*(rec._arr.get()))[rec._cnt]) RawTower_type();

                  if (Verbosity() >= 2)
                    cout
                        << "Prototype2DSTReader::process_event - processing Tower hit "
                        << rec._name << " @ ("
//                        << hit->get_thetaMin() << ", "
//                        << hit->get_phiMin()
                        << "), Energy = " << hit->get_energy() << " - "
                        << rec._arr.get()->At(rec._cnt)->ClassName() << endl;

//                  rec._arr->Print();

                  RawTower_type * new_hit =
                      dynamic_cast<RawTower_type *>(rec._arr.get()->At(rec._cnt));
                  assert(new_hit);

                  *new_hit = (*hit);

                  rec._cnt++;
                }
            } // if (!hits)
        } //      if (rec._type == record::typ_hit)
      else if (rec._type == record::typ_towertemp)
        {
          if (Verbosity() >= 2)
            cout
                << "Prototype2DSTReader::process_event - processing tower temperature "
                << rec._name << endl;

          RawTowerContainer *hits = findNode::getClass<RawTowerContainer>(
              topNode, rec._name);
          if (!hits)
            {
              if (_event < 2)
                cout
                    << "Prototype2DSTReader::process_event - Error - can not find node "
                    << rec._name << endl;

            }
          else
            {
              RawTowerContainer::ConstRange hit_range = hits->getTowers();

              if (Verbosity() >= 2)
                cout << "Prototype2DSTReader::process_event - processing "
                    << rec._name << " and received " << hits->size()
                    << " tower hits" << endl;

              rec._cnt = 0;

              for (RawTowerContainer::ConstIterator hit_iter = hit_range.first;
                  hit_iter != hit_range.second; hit_iter++)
                {
                  RawTower * hit_raw = hit_iter->second;

                  RawTowerT_type * hit = dynamic_cast<RawTowerT_type *>(hit_raw);
//                  RawTower * hit = hit_iter->second;

                  assert(hit);

                  rec._dvalue += hit->get_temperature_from_entry(
                      hit->get_nr_entries() - 1);

                  ++rec._cnt;
                }

              rec._dvalue /= rec._cnt;
            } // if (!hits)
        } //      if (rec._type == record::typ_hit)
      else if (rec._type == record::typ_runinfo)
        {

          PdbParameterMap *info = findNode::getClass<PdbParameterMap>(topNode,
              "RUN_INFO");

          assert(info);

          PHParameters run_info_copy("RunInfo");
          run_info_copy.FillFrom(info);

          rec._dvalue = run_info_copy.get_double_param(rec._name);

        } //      if (rec._type == record::typ_hit)
      else if (rec._type == record::typ_part)
        {
          assert(0);
        } //      if (rec._type == record::typ_part)
      else if (rec._type == record::typ_vertex)
        {
          assert(0);
        } //          if (_load_all_particle)

    } //  for (records_t::iterator it = _records.begin(); it != _records.end(); ++it)

  if (_T)
    _T->Fill();

  return 0;
} //  for (records_t::iterator it = _records.begin(); it != _records.end(); ++it)

int
Prototype2DSTReader::End(PHCompositeNode * /*topNode*/)
{
  cout << "Prototype2DSTReader::End - Clean ups" << endl;

  if (_T)
    {
      PHTFileServer::get().cd(_out_file_name);
      _T->Write();
      _T->ResetBranchAddresses();
    }

  _records.clear();

  return 0;
}
