// $Id: PHTimeServer.h,v 1.5 2014/01/12 16:15:05 pinkenbu Exp $
#ifndef __PHTIMESERVER_H__
#define __PHTIMESERVER_H__

/*!
\file		PHTimeServer.h
\brief	 PHTimer server for accessing external information
\author	Hugo Pereira
\version $Revision: 1.5 $
\date		$Date: 2014/01/12 16:15:05 $
*/

#include "PHTimer.h"

#include <iostream>
#include <memory>
#include <map>
#include <string>

#ifndef __CINT__
#include <boost/smart_ptr.hpp>
#endif
/*!
\class	 PHTimeServer
\brief	 PHTimer server for accessing external information
*/

class PHTimeServer
{
 public:
  //! wrapper around PHTimer, for storage in a map
  class timer
  {
   public:
    timer(const std::string& key, unsigned short uid)
      : _timer(new PHTimer(key))
      , _uid(uid)
    {
    }

    PHTimer* get(void)
    {
      return _timer.get();
    }

    const PHTimer* get(void) const
    {
      return _timer.get();
    }

    unsigned short get_uid(void) const
    {
      return _uid;
    }

   private:
#ifndef __CINT__
    std::shared_ptr<PHTimer> _timer;
#endif
    unsigned short _uid;
  };

  //! singleton accessor
  static PHTimeServer* get(void)
  {
    static PHTimeServer* _server = new PHTimeServer();
    return _server;
  }

  //! destructor
  virtual ~PHTimeServer()
  {
  }

  //! insert new timer in map.
  timer insert_new(const std::string&);

  //! insert new single_shot timer in map.
  timer insert_new_single_shot(const std::string&);

  //! retrieve existing timer. throw exception if not found
  timer get_timer(const std::string&);

  //! retrieve existing timer. throw exception if not found
  timer get_single_shot_timer(const std::string&);

  //! dump all registered timer value.
  void print(std::ostream& out = std::cout) const;

  //! dump all registered timer statistics.
  void print_stat(std::ostream& out = std::cout) const;

 protected:
  //! constructor
  PHTimeServer()
    : _timer_id(0)
    , _single_shot_timer_id(0)
  {
  }

 private:
  //! map
  typedef std::map<std::string, timer> time_map;

  //! iterator over the map.
  typedef std::map<std::string, timer>::iterator time_iterator;

  //! iterator over the map.
  typedef std::map<std::string, timer>::const_iterator const_time_iterator;

  //! list of timers
  time_map _timers;

  //! list of single shot timers
  time_map _single_shot_timers;

  //! running timer unique id
  unsigned short _timer_id;

  //! running single shot timer unique id
  unsigned short _single_shot_timer_id;

 public:
  //! light iterator over PHTimer map
  class iterator
  {
   public:
    //! get PHTimer associated to current iterator position, advance iterator
    PHTimeServer::timer* next()
    {
      if (_iter == _map.end()) return 0;
      PHTimeServer::timer* out(&_iter->second);
      _iter++;
      return out;
    }

    //! get PHTimer associated to current iterator position
    PHTimeServer::timer* current()
    {
      if (_iter == _map.end()) return 0;
      return &_iter->second;
    }

   protected:
    //! creator
    iterator(PHTimeServer::time_map map)
      : _map(map)
      , _iter(_map.begin())
    {
    }

   private:
    //! map of PHTimer
    PHTimeServer::time_map _map;

    //! iterator over the map
    PHTimeServer::time_iterator _iter;

    friend class PHTimeServer;
  };

  //! return iterator over the map, located at begin
  iterator range(void) { return iterator(_timers); }
};
#endif
