// $Id: PHG4DSTReader.cc,v 1.11 2015/01/06 02:52:07 jinhuang Exp $

/*!
 * \file PHG4DSTReader.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision: 1.11 $
 * \date $Date: 2015/01/06 02:52:07 $
 */

#include "PHG4DSTReader.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerv1.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitEval.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4VtxPointv1.h>

#include <g4jets/Jet.h>
#include <g4jets/JetMap.h>
#include <g4jets/Jetv1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>

#include <TClass.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TString.h>
#include <TTree.h>

#include <cassert>
#include <iostream>
#include <set>
#include <sstream>
#include <utility>

using namespace std;

typedef PHG4Particlev2 part_type;
typedef PHG4HitEval hit_type;
typedef PHG4VtxPointv1 vertex_type;
typedef RawTowerv1 RawTower_type;
typedef Jetv1 PHPyJet_type;

PHG4DSTReader::PHG4DSTReader(const string &filename)
  : SubsysReco("PHG4DSTReader")
  , nblocks(0)
  , _event(0)
  ,  //
  _out_file_name(filename)
  , /*_file(nullptr), */ _T(nullptr)
  ,  //
  _save_particle(true)
  , _load_all_particle(false)
  , _load_active_particle(
        true)
  , _save_vertex(true)
  , _tower_zero_sup(.0)
{
  // TODO Auto-generated constructor stub
}

PHG4DSTReader::~PHG4DSTReader()
{
  if (Verbosity() > 1)
  {
    cout << "PHG4DSTReader::destructor - Clean ups" << endl;
  }
  if (_T)
  {
    _T->ResetBranchAddresses();
  }

  _records.clear();
}

int PHG4DSTReader::Init(PHCompositeNode *)
{
  const int arr_size = 100;

  //  BOOST_FOREACH(string nodenam, _node_postfix)
  for (vector<string>::const_iterator it = _node_postfix.begin();
       it != _node_postfix.end(); ++it)
  {
    const char *class_name = hit_type::Class()->GetName();

    const string &nodenam = *it;

    string hname = Form("G4HIT_%s", nodenam.c_str());
    //      _node_name.push_back(hname);
    if (Verbosity() > 0)
    {
      cout << "PHG4DSTReader::Init - saving hits from node: " << hname << " - "
           << class_name << endl;
    }

    record rec;
    rec._cnt = 0;
    rec._name = hname;
    rec._arr = boost::make_shared<TClonesArray>(class_name, arr_size);
    rec._arr_ptr = rec._arr.get();
    rec._type = record::typ_hit;

    _records.push_back(rec);

    nblocks++;
  }

  if (_tower_postfix.size() && Verbosity() > 0)
  {
    cout << "PHG4DSTReader::Init - zero suppression for calorimeter towers = "
         << _tower_zero_sup << " GeV" << endl;
  }
  for (vector<string>::const_iterator it = _tower_postfix.begin();
       it != _tower_postfix.end(); ++it)
  {
    const char *class_name = RawTower_type::Class()->GetName();

    const string &nodenam = *it;

    string hname = Form("TOWER_%s", nodenam.c_str());
    //      _node_name.push_back(hname);
    if (Verbosity() > 0)
    {
      cout << "PHG4DSTReader::Init - saving raw tower info from node: " << hname
           << " - " << class_name << endl;
    }
    record rec;
    rec._cnt = 0;
    rec._name = hname;
    rec._arr = boost::make_shared<TClonesArray>(class_name, arr_size);
    rec._arr_ptr = rec._arr.get();
    rec._type = record::typ_tower;

    _records.push_back(rec);

    nblocks++;
  }

  for (vector<string>::const_iterator it = _jet_postfix.begin();
       it != _jet_postfix.end(); ++it)
  {
    const char *class_name = PHPyJet_type::Class()->GetName();

    const string &hname = *it;
    if (Verbosity() > 0)
    {
      cout << "PHG4DSTReader::Init - saving jets from node: " << hname << " - "
           << class_name << endl;
    }
    record rec;
    rec._cnt = 0;
    rec._name = hname;
    rec._arr = boost::make_shared<TClonesArray>(class_name, arr_size);
    rec._arr_ptr = rec._arr.get();
    rec._type = record::typ_jets;

    _records.push_back(rec);

    nblocks++;
  }

  if (_save_particle)
  {
    //save particles

    const char *class_name = part_type::Class()->GetName();

    if (Verbosity() > 0)
    {
      cout << "PHG4DSTReader::Init - saving Particles node:" << class_name
           << endl;
    }
    record rec;
    rec._cnt = 0;
    rec._name = "PHG4Particle";
    rec._arr = boost::make_shared<TClonesArray>(class_name, arr_size);
    rec._arr_ptr = rec._arr.get();
    rec._type = record::typ_part;

    _records.push_back(rec);

    nblocks++;
  }

  if (_save_vertex)
  {
    //save particles

    const char *class_name = vertex_type::Class()->GetName();

    if (Verbosity() > 0)
    {
      cout << "PHG4DSTReader::Init - saving vertex for particles under node:"
           << class_name << endl;
    }
    record rec;
    rec._cnt = 0;
    rec._name = "PHG4VtxPoint";
    rec._arr = boost::make_shared<TClonesArray>(class_name, arr_size);
    rec._arr_ptr = rec._arr.get();
    rec._type = record::typ_vertex;

    _records.push_back(rec);

    nblocks++;
  }
  if (Verbosity() > 0)
  {
    cout << "PHG4DSTReader::Init - requested " << nblocks << " nodes" << endl;
  }
  build_tree();

  return 0;
}

void PHG4DSTReader::build_tree()
{
  if (Verbosity() > 0)
  {
    cout << "PHG4DSTReader::build_tree - output to " << _out_file_name << endl;
  }
  static const int BUFFER_SIZE = 32000;

  // open TFile
  PHTFileServer::get().open(_out_file_name, "RECREATE");

  _T = new TTree("T", "PHG4DSTReader");

  nblocks = 0;
  for (records_t::iterator it = _records.begin(); it != _records.end(); ++it)
  {
    record &rec = *it;

    if (Verbosity() > 0)
    {
      cout << "PHG4DSTReader::build_tree - Add " << rec._name << endl;
    }
    const string name_cnt = "n_" + rec._name;
    const string name_cnt_desc = name_cnt + "/I";
    _T->Branch(name_cnt.c_str(), &(rec._cnt), name_cnt_desc.c_str(),
               BUFFER_SIZE);
    _T->Branch(rec._name.c_str(), &(rec._arr_ptr), BUFFER_SIZE, 99);

    nblocks++;
  }

  if (Verbosity() > 0)
  {
    cout << "PHG4DSTReader::build_tree - added " << nblocks << " nodes" << endl;
  }
  _T->SetAutoSave(16000);
}

int PHG4DSTReader::process_event(PHCompositeNode *topNode)
{
  //  const double significand = _event / TMath::Power(10, (int) (log10(_event)));
  //
  //  if (fmod(significand, 1.0) == 0 && significand <= 10)
  //    cout << "PHG4DSTReader::process_event - " << _event << endl;
  _event++;

  //clean ups
  _particle_set.clear();
  _vertex_set.clear();

  PHG4TruthInfoContainer *truthInfoList = findNode::getClass<
      PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!truthInfoList)
  {
    if (_event < 2)
      cout
          << "PHG4DSTReader::process_event - Error - can not find node G4TruthInfo. Quit processing!"
          << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  for (records_t::iterator it = _records.begin(); it != _records.end(); ++it)
  {
    record &rec = *it;

    rec._cnt = 0;
    assert(rec._arr.get() == rec._arr_ptr);
    assert(rec._arr.get());
    rec._arr->Clear();

    if (rec._type == record::typ_hit)
    {
      if (Verbosity() >= 2)
        cout << "PHG4DSTReader::process_event - processing " << rec._name
             << endl;

      PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode,
                                                                    rec._name);
      if (!hits)
      {
        if (_event < 2)
          cout
              << "PHG4DSTReader::process_event - Error - can not find node "
              << rec._name << endl;
      }
      else
      {
        PHG4HitContainer::ConstRange hit_range = hits->getHits();

        if (Verbosity() >= 2)
          cout << "PHG4DSTReader::process_event - processing "
               << rec._name << " and received " << hits->size() << " hits"
               << endl;

        for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first;
             hit_iter != hit_range.second; hit_iter++)
        {
          PHG4Hit *hit = hit_iter->second;

          //                  hit_type * hit = dynamic_cast<hit_type *>(hit_raw);
          //
          assert(hit);

          new ((*(rec._arr.get()))[rec._cnt]) hit_type(hit);

          hit_type *new_hit =
              dynamic_cast<hit_type *>(rec._arr.get()->At(rec._cnt));
          assert(new_hit);

          //                  for (int i = 0; i < 2; i++)
          //                    {
          //                      new_hit->set_x(i, hit->get_x(i));
          //                      new_hit->set_y(i, hit->get_y(i));
          //                      new_hit->set_z(i, hit->get_z(i));
          //                      new_hit->set_t(i, hit->get_t(i));
          //                    }
          //                  new_hit->set_edep(hit->get_edep());
          //                  new_hit->set_layer(hit->get_layer());
          //                  new_hit->set_hit_id(hit->get_hit_id());
          //                  new_hit->set_trkid(hit->get_trkid());

          //                  *new_hit = (*hit);

          if (_load_active_particle)
            _particle_set.insert(hit->get_trkid());

          if (Verbosity() >= 2)
            cout << "PHG4DSTReader::process_event - processing "
                 << rec._name << " and hit " << hit->get_hit_id()
                 << " with track id " << hit->get_trkid() << endl;

          rec._cnt++;
        }
      }  // if (!hits)
    }    //      if (rec._type == record::typ_hit)
    else if (rec._type == record::typ_tower)
    {
      if (Verbosity() >= 2)
        cout << "PHG4DSTReader::process_event - processing tower "
             << rec._name << endl;

      RawTowerContainer *hits = findNode::getClass<RawTowerContainer>(
          topNode, rec._name);
      if (!hits)
      {
        if (_event < 2)
          cout
              << "PHG4DSTReader::process_event - Error - can not find node "
              << rec._name << endl;
      }
      else
      {
        RawTowerContainer::ConstRange hit_range = hits->getTowers();

        if (Verbosity() >= 2)
          cout << "PHG4DSTReader::process_event - processing "
               << rec._name << " and received " << hits->size()
               << " tower hits" << endl;

        for (RawTowerContainer::ConstIterator hit_iter = hit_range.first;
             hit_iter != hit_range.second; hit_iter++)
        {
          RawTower *hit_raw = hit_iter->second;

          RawTower_type *hit = dynamic_cast<RawTower_type *>(hit_raw);
          //                  RawTower * hit = hit_iter->second;

          assert(hit);

          if (hit->get_energy() < _tower_zero_sup)
          {
            if (Verbosity() >= 2)
              cout
                  << "PHG4DSTReader::process_event - suppress low energy tower hit "
                  << rec._name << " @ ("
                  //                            << hit->get_thetaMin()
                  //                            << ", " << hit->get_phiMin()
                  << "), Energy = " << hit->get_energy() << endl;

            continue;
          }

          new ((*(rec._arr.get()))[rec._cnt]) RawTower_type();

          if (Verbosity() >= 2)
            cout
                << "PHG4DSTReader::process_event - processing Tower hit "
                << rec._name << " @ ("
                //                        << hit->get_thetaMin() << ", "
                //                        << hit->get_phiMin()
                << "), Energy = " << hit->get_energy() << " - "
                << rec._arr.get()->At(rec._cnt)->ClassName() << endl;

          //                  rec._arr->Print();

          RawTower_type *new_hit =
              dynamic_cast<RawTower_type *>(rec._arr.get()->At(rec._cnt));
          assert(new_hit);

          //                  //v1 copy
          //                  new_hit->get_bineta(hit->get_bineta());
          //                  new_hit->get_binphi(hit->get_binphi());
          //
          //                  //v2 copy
          //                  new_hit->set_thetaMin(hit->get_thetaMin());
          //                  new_hit->set_thetaSize(hit->get_thetaSize());
          //                  new_hit->set_phiMin(hit->get_phiMin());
          //                  new_hit->set_phiSize(hit->get_phiSize());
          //                  new_hit->set_zMin(hit->get_zMin());
          //                  new_hit->set_zSize(hit->get_zSize());

          *new_hit = (*hit);

          rec._cnt++;
        }
      }  // if (!hits)
    }    //      if (rec._type == record::typ_hit)
    else if (rec._type == record::typ_jets)
    {
      //          std::cout
      //              << "PHG4DSTReader::AddJet - Error - temp. disabled until jet added back to sPHENIX software"
      //              << std::endl;
      //
      if (Verbosity() >= 2)
        cout << "PHG4DSTReader::process_event - processing jets "
             << rec._name << endl;

      JetMap *hits = findNode::getClass<JetMap>(topNode, rec._name);
      if (!hits)
      {
        if (_event < 2)
          cout
              << "PHG4DSTReader::process_event - Error - can not find node "
              << rec._name << endl;
      }
      else
      {
        if (Verbosity() >= 2)
          cout << "PHG4DSTReader::process_event - processing "
               << rec._name << " and received " << hits->size() << " jets"
               << endl;

        // for every recojet
        for (JetMap::Iter iter = hits->begin(); iter != hits->end();
             ++iter)
        {
          Jet *hit_raw = iter->second;

          if (Verbosity() >= 2)
            cout << "PHG4DSTReader::process_event - processing jet "
                 << rec._name << " @ (" << hit_raw->get_eta() << ", "
                 << hit_raw->get_phi() << "), pT = " << hit_raw->get_pt()
                 << " - with raw type " << hit_raw->ClassName() << endl;

          PHPyJet_type *hit = dynamic_cast<PHPyJet_type *>(hit_raw);

          assert(hit);

          new ((*(rec._arr.get()))[rec._cnt]) PHPyJet_type();

          PHPyJet_type *new_hit =
              dynamic_cast<PHPyJet_type *>(rec._arr.get()->At(rec._cnt));
          assert(new_hit);

          *new_hit = (*hit);

          rec._cnt++;
        }
      }  // if (!hits)
    }    //      if (rec._type == record::typ_hit)
    else if (rec._type == record::typ_part)
    {
      if (_load_all_particle)
      {
        static bool once = true;
        if (once)
        {
          cout
              << "PHG4DSTReader::process_event - will load all particle from G4TruthInfo"
              << endl;

          once = false;
        }

        for (auto particle_iter = truthInfoList->GetMap().begin();
             particle_iter != truthInfoList->GetMap().end();
             ++particle_iter)
        {
          _particle_set.insert(particle_iter->first);
        }

      }  //          if (_load_all_particle)
      else
      {
        static bool once = true;
        if (once)
        {
          if (Verbosity() > 0)
          {
            cout << "PHG4DSTReader::process_event - will load primary particle from G4TruthInfo" << endl;
          }
          once = false;
        }

        PHG4TruthInfoContainer::ConstRange primary_range =
            truthInfoList->GetPrimaryParticleRange();

        for (PHG4TruthInfoContainer::ConstIterator particle_iter =
                 primary_range.first;
             particle_iter != primary_range.second;
             ++particle_iter)
        {
          _particle_set.insert(particle_iter->first);

        }  //          if (_load_all_particle) else
      }
      for (PartSet_t::const_iterator i = _particle_set.begin();
           i != _particle_set.end(); ++i)
      {
        auto particle = truthInfoList->GetParticle(*i);
        if (!particle)
        {
          cout
              << "PHG4DSTReader::process_event - ERROR - can not find particle/track ID "
              << *i << " in G4TruthInfo" << endl;

          continue;
        }

        add_particle(rec, particle);
      }  // for(PartSet_t::const_iterator i = _particle_set.begin();i!=_particle_set.end();i++)

    }  //      if (rec._type == record::typ_part)
    else if (rec._type == record::typ_vertex)
    {
      static bool once = true;
      if (once)
      {
        if (Verbosity() > 0)
        {
          cout << "PHG4DSTReader::process_event - will load vertex from G4TruthInfo" << endl;
        }
        once = false;
      }

      for (PartSet_t::const_iterator i = _vertex_set.begin();
           i != _vertex_set.end(); ++i)
      {
        PHG4VtxPoint *v = truthInfoList->GetVtx(*i);
        if (!v)
        {
          cout
              << "PHG4DSTReader::process_event - ERROR - can not find vertex ID "
              << *i << " in G4TruthInfo" << endl;

          continue;
        }

        new ((*(rec._arr.get()))[rec._cnt]) vertex_type();

        if (Verbosity() >= 2)
          cout << "PHG4DSTReader::process_event - saving vertex id "
               << *i << endl;

        vertex_type *new_v =
            static_cast<vertex_type *>(rec._arr.get()->At(rec._cnt));
        assert(new_v);

        new_v->set_x(v->get_x());
        new_v->set_y(v->get_y());
        new_v->set_z(v->get_z());
        new_v->set_t(v->get_t());
        new_v->set_id(v->get_id());

        rec._cnt++;
      }  // else if (rec._type == record::typ_vertex)

    }  //          if (_load_all_particle)

  }  //  for (records_t::iterator it = _records.begin(); it != _records.end(); ++it)

  if (_T)
    _T->Fill();

  return 0;
}  //  for (records_t::iterator it = _records.begin(); it != _records.end(); ++it)

void PHG4DSTReader::add_particle(PHG4DSTReader::record &rec, PHG4Particle *part)
{
  assert(part);

  //              if (Verbosity() >= 2)
  //                cout << "PHG4DSTReader::process_event - Particle type is "
  //                    << p_raw->GetName() << " - " << p_raw->ClassName()
  //                    << " get_track_id = " << p_raw->get_track_id() << endl;

  //              part_type * part = dynamic_cast<part_type *>(p_raw);

  //              assert(part);

  new ((*(rec._arr.get()))[rec._cnt]) part_type();

  part_type *new_part = static_cast<part_type *>(rec._arr.get()->At(rec._cnt));
  assert(new_part);

  new_part->set_track_id(part->get_track_id());
  new_part->set_vtx_id(part->get_vtx_id());
  new_part->set_parent_id(part->get_parent_id());
  new_part->set_primary_id(part->get_primary_id());
  new_part->set_name(part->get_name());
  new_part->set_pid(part->get_pid());
  new_part->set_px(part->get_px());
  new_part->set_py(part->get_py());
  new_part->set_pz(part->get_pz());
  new_part->set_e(part->get_e());

  _vertex_set.insert(part->get_vtx_id());

  rec._cnt++;
}

int PHG4DSTReader::End(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 1)
  {
    cout << "PHG4DSTReader::End - Clean ups" << endl;
  }

  if (_T)
  {
    PHTFileServer::get().cd(_out_file_name);
    _T->Write();
    _T->ResetBranchAddresses();
  }

  _records.clear();

  return 0;
}
