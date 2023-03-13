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

//_________________________________________________________
PHTimeServer::timer PHTimeServer::insert_new(const std::string& key)
{
  std::string tmp_key(key);

  int version = 0;
  while ((_timers.find(tmp_key)) != _timers.end())
  {
    version++;
    std::ostringstream o;
    o << key << "_" << version;
    tmp_key = o.str();
  }

  // create a new timer
  _timers.insert(std::pair<std::string, timer>(tmp_key, timer(tmp_key, _timer_id)));
  _timer_id++;

  // returns timer
  return _timers.find(tmp_key)->second;
}
//_________________________________________________________
PHTimeServer::timer PHTimeServer::insert_new_single_shot(const std::string& key)
{
  std::string tmp_key(key);

  int version = 0;
  while ((_single_shot_timers.find(tmp_key)) != _single_shot_timers.end())
  {
    version++;
    std::ostringstream o;
    o << key << "_" << version;
    tmp_key = o.str();
  }

  // create a new timer
  _single_shot_timers.insert(std::pair<std::string, timer>(tmp_key, timer(tmp_key, _single_shot_timer_id)));
  _single_shot_timer_id++;

  // returns timer
  return _single_shot_timers.find(tmp_key)->second;
}

//_________________________________________________________
PHTimeServer::timer PHTimeServer::get_timer(const std::string& key)
{
  // check for existing timer matching key
  time_iterator _iter = _timers.find(key);
  if (_iter != _timers.end())
    return _iter->second;
  else
  {
    std::ostringstream what;
    what << "unknown timer \"" << key << "\" requested.";
    throw std::invalid_argument(what.str());
  }
}

//_________________________________________________________
PHTimeServer::timer PHTimeServer::get_single_shot_timer(const std::string& key)
{
  // check for existing timer matching key
  time_iterator _iter = _single_shot_timers.find(key);
  if (_iter != _single_shot_timers.end())
    return _iter->second;
  else
  {
    std::ostringstream what;
    what << "unknown timer \"" << key << "\" requested.";
    throw std::invalid_argument(what.str());
  }
}

//_________________________________________________________
void PHTimeServer::print(std::ostream& out) const
{
  PHTimer::PRINT(out, "Mutoo PHTimeServer");

  // run over normal timers
  for (const auto & _timer : _timers)
  {
    char str[512];
    sprintf(str, "%-20s [%2i] - %-6g ms (%s)-.",
            _timer.second.get()->get_name().c_str(),
            _timer.second.get_uid(),
            _timer.second.get()->elapsed(),
            (char*) ((_timer.second.get()->get_state() == PHTimer::RUN) ? " (running)" : " (stopped)"));
    out << str << std::endl;
  }

  // run over single_shot timers
  PHTimer::PRINT(out, "Mutoo PHTimeServer - single_shots");
  for (const auto & _single_shot_timer : _single_shot_timers)
  {
    char str[512];
    sprintf(str, "single_shot - %-20s [%2i] - %-6g ms (%s)-.",
            _single_shot_timer.second.get()->get_name().c_str(),
            _single_shot_timer.second.get_uid(),
            _single_shot_timer.second.get()->elapsed(),
            (char*) ((_single_shot_timer.second.get()->get_state() == PHTimer::RUN) ? " (running)" : " (stopped)"));
    out << str << std::endl;
  }

  PHTimer::PRINT(out, "**");

  return;
}

//_________________________________________________________
void PHTimeServer::print_stat(std::ostream& out) const
{
  // print nothing if no timer was registered
  if (_timers.empty() && _single_shot_timers.empty()) return;

  // header
  PHTimer::PRINT(out, "Mutoo PHTimeServer statistics");

  // normal timers
  for (const auto & _timer : _timers)
    if (_timer.second.get()->get_ncycle())
    {
      char str[512];
      sprintf(str, "%-20s [%2i] - Accumulated time: %-6g ms. cycles: %-10u. Time per cycle: %-6g ms",
              _timer.second.get()->get_name().c_str(),
              _timer.second.get_uid(),
              _timer.second.get()->get_accumulated_time(),
              _timer.second.get()->get_ncycle(),
              _timer.second.get()->get_time_per_cycle());
      out << str << std::endl;
    }

  // single shot timers
  PHTimer::PRINT(out, "Mutoo PHTimeServer single_shots statistics");
  for (const auto & _single_shot_timer : _single_shot_timers)
    if (_single_shot_timer.second.get()->get_ncycle())
    {
      char str[512];
      sprintf(str, "single_shot - %-20s [%2i] - accumulated: %-6g ms.",
              _single_shot_timer.second.get()->get_name().c_str(),
              _single_shot_timer.second.get_uid(),
              _single_shot_timer.second.get()->get_accumulated_time());
      out << str;

      // check timer _was_ single shot
      if (_single_shot_timer.second.get()->get_ncycle() != 1)
        out << " WARNING: single_shot started more than once.";

      out << std::endl;
    }

  PHTimer::PRINT(out, "**");

  return;
}
