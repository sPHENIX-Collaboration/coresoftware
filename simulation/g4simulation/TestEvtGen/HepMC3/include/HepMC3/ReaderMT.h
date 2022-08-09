// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READERMT_H
#define HEPMC3_READERMT_H
///
/// @file  ReaderMT.h
/// @brief Definition of class \b ReaderMT
///
/// @class HepMC3::ReaderMT
/// @brief Multithreader GenEvent I/O parsing
///
/// @ingroup IO
///
#include <set>
#include <string>
#include <fstream>
#include <istream>
#include <iterator>
#include <thread>
#include "HepMC3/Reader.h"
#include "HepMC3/GenEvent.h"
namespace HepMC3 {
template <class T, size_t m_number_of_threads>  class ReaderMT : public Reader
{
private:
    bool m_go_try_cache; //!< Flag to trigger using the cached event
    std::vector< std::shared_ptr<T> > m_readers; //!< Vector of all active readers
    std::vector< std::pair<GenEvent, bool> > m_events; //!< Vector of events
    std::vector< std::thread > m_threads;
    static void read_function(std::pair<GenEvent,bool>& e, std::shared_ptr<T> r)
    {
        e.second = r->read_event(e.first);
        r->skip(m_number_of_threads-1);
        if (r->failed()) r->close();
    }
public:
    ReaderMT(const std::string& filename): m_go_try_cache(true) {
        m_events.reserve(m_number_of_threads);
        m_readers.reserve(m_number_of_threads);
        m_threads.reserve(m_number_of_threads);
        for (size_t i = 0; i < m_number_of_threads; ++i) {
            m_readers.push_back(std::make_shared<T>(filename));
            m_readers.back()->skip(m_number_of_threads-1-i);
        }
    }
    ~ReaderMT() {
        m_readers.clear();
        m_events.clear();
        m_threads.clear();
    }
    bool skip(const int) override  {
        return false;///Not implemented
    }
    bool read_event(GenEvent& evt)  override {
        if ( !m_events.empty() ) {
            evt = m_events.back().first;
            m_events.pop_back();
            return true;
        }
        m_events.clear();
        m_threads.clear();
        m_go_try_cache = true;
        m_threads.reserve(m_number_of_threads);
        m_events.reserve(m_number_of_threads);
        for (size_t i = 0; i < m_number_of_threads; ++i) {
            m_events.push_back(std::pair<GenEvent, bool>(GenEvent(Units::GEV,Units::MM), true));
            m_threads.push_back(std::thread(read_function, std::ref(m_events.at(i)), m_readers.at(i)));
        }
        for (auto& th : m_threads) {
            th.join();
        }
        m_threads.clear();

        m_events.erase(std::remove_if(m_events.begin(), m_events.end(),[](std::pair<GenEvent, bool>& x) {
            return !x.second;
        }), m_events.end());

        if (m_events.empty()) {
            m_go_try_cache = false;
            return false;
        }
        evt = m_events.back().first;
        m_events.pop_back();
        return true;
    }
    bool failed()  override {
        for (auto& reader: m_readers)    if (reader && !reader->failed()) return false;
        if ( !m_events.empty() ) return false;
        if ( m_go_try_cache ) return false;
        return true;
    }
    void close()   override {
        for (auto& reader: m_readers) if (reader) reader->close();
    }
};
}
#endif
