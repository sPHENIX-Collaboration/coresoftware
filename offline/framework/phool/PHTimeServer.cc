// $Id: PHTimeServer.C,v 1.5 2014/01/12 04:14:40 pinkenbu Exp $

/*!
\file		PHTimeServer.cxx
\brief	 PHTimer server for accessing external information
\author	Hugo Pereira
\version $Revision: 1.5 $
\date		$Date: 2014/01/12 04:14:40 $
*/

#include "PHTimeServer.h"

#include <cstdio>
#include <stdexcept>

using namespace std;

//_________________________________________________________
PHTimeServer::timer PHTimeServer::insert_new(const string& key)
{
  string tmp_key(key);

  int version = 0;
  while ((_timers.find(tmp_key)) != _timers.end())
  {
    version++;
    ostringstream o;
    o << key << "_" << version;
    tmp_key = o.str();
  }

  // create a new timer
  _timers.insert(pair<string, timer>(tmp_key, timer(tmp_key, _timer_id)));
  _timer_id++;

  // returns timer
  return _timers.find(tmp_key)->second;
}
//_________________________________________________________
PHTimeServer::timer PHTimeServer::insert_new_single_shot(const string& key)
{
  string tmp_key(key);

  int version = 0;
  while ((_single_shot_timers.find(tmp_key)) != _single_shot_timers.end())
  {
    version++;
    ostringstream o;
    o << key << "_" << version;
    tmp_key = o.str();
  }

  // create a new timer
  _single_shot_timers.insert(pair<string, timer>(tmp_key, timer(tmp_key, _single_shot_timer_id)));
  _single_shot_timer_id++;

  // returns timer
  return _single_shot_timers.find(tmp_key)->second;
}

//_________________________________________________________
PHTimeServer::timer PHTimeServer::get_timer(const string& key)
{
  // check for existing timer matching key
  time_iterator _iter = _timers.find(key);
  if (_iter != _timers.end())
    return _iter->second;
  else
  {
    ostringstream what;
    what << "unknown timer \"" << key << "\" requested.";
    throw invalid_argument(what.str());
  }
}

//_________________________________________________________
PHTimeServer::timer PHTimeServer::get_single_shot_timer(const string& key)
{
  // check for existing timer matching key
  time_iterator _iter = _single_shot_timers.find(key);
  if (_iter != _single_shot_timers.end())
    return _iter->second;
  else
  {
    ostringstream what;
    what << "unknown timer \"" << key << "\" requested.";
    throw invalid_argument(what.str());
  }
}

//_________________________________________________________
void PHTimeServer::print(ostream& out) const
{
  PHTimer::PRINT(out, "Mutoo PHTimeServer");

  // run over normal timers
  for (const_time_iterator iter = _timers.begin(); iter != _timers.end(); ++iter)
  {
    char str[512];
    sprintf(str, "%-20s [%2i] - %-6g ms (%s)-.",
            iter->second.get()->get_name().c_str(),
            iter->second.get_uid(),
            iter->second.get()->elapsed(),
            (char*) ((iter->second.get()->get_state() == PHTimer::RUN) ? " (running)" : " (stopped)"));
    out << str << endl;
  }

  // run over single_shot timers
  PHTimer::PRINT(out, "Mutoo PHTimeServer - single_shots");
  for (const_time_iterator iter = _single_shot_timers.begin(); iter != _single_shot_timers.end(); ++iter)
  {
    char str[512];
    sprintf(str, "single_shot - %-20s [%2i] - %-6g ms (%s)-.",
            iter->second.get()->get_name().c_str(),
            iter->second.get_uid(),
            iter->second.get()->elapsed(),
            (char*) ((iter->second.get()->get_state() == PHTimer::RUN) ? " (running)" : " (stopped)"));
    out << str << endl;
  }

  PHTimer::PRINT(out, "**");

  return;
}

//_________________________________________________________
void PHTimeServer::print_stat(ostream& out) const
{
  // print nothing if no timer was registered
  if (_timers.empty() && _single_shot_timers.empty()) return;

  // header
  PHTimer::PRINT(out, "Mutoo PHTimeServer statistics");

  // normal timers
  for (const_time_iterator iter = _timers.begin(); iter != _timers.end(); ++iter)
    if (iter->second.get()->get_ncycle())
    {
      char str[512];
      sprintf(str, "%-20s [%2i] - Accumulated time: %-6g ms. cycles: %-10u. Time per cycle: %-6g ms",
              iter->second.get()->get_name().c_str(),
              iter->second.get_uid(),
              iter->second.get()->get_accumulated_time(),
              iter->second.get()->get_ncycle(),
              iter->second.get()->get_time_per_cycle());
      out << str << endl;
    }

  // single shot timers
  PHTimer::PRINT(out, "Mutoo PHTimeServer single_shots statistics");
  for (const_time_iterator iter = _single_shot_timers.begin(); iter != _single_shot_timers.end(); ++iter)
    if (iter->second.get()->get_ncycle())
    {
      char str[512];
      sprintf(str, "single_shot - %-20s [%2i] - accumulated: %-6g ms.",
              iter->second.get()->get_name().c_str(),
              iter->second.get_uid(),
              iter->second.get()->get_accumulated_time());
      out << str;

      // check timer _was_ single shot
      if (iter->second.get()->get_ncycle() != 1)
        out << " WARNING: single_shot started more than once.";

      out << endl;
    }

  PHTimer::PRINT(out, "**");

  return;
}
