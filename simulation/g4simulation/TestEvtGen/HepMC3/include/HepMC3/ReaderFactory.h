// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READERFACTORY_H
#define HEPMC3_READERFACTORY_H

#include <memory>
#include <string>
#include <sys/stat.h>
#include <string.h>

#include "HepMC3/ReaderAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
#include "HepMC3/ReaderHEPEVT.h"
#include "HepMC3/ReaderLHEF.h"
#include "HepMC3/ReaderPlugin.h"
#include "HepMC3/CompressedIO.h"

namespace HepMC3 {


std::shared_ptr<Reader> deduce_reader(std::istream &stream);

std::shared_ptr<Reader> deduce_reader(std::shared_ptr<std::istream> stream);

/**
 * @brief This function deduces the type of input file based on the name/URL
 * and its content, and will return an appropriate Reader object
 *
 * @todo Too big for inlining: move into a .cc implementation file?
 * @todo Need a DEBUG verbosity flag
 */
std::shared_ptr<Reader> deduce_reader(const std::string &filename)
{
    std::string libHepMC3rootIO = "libHepMC3rootIO.so.3";
#if defined(__darwin__) || defined(__APPLE__)
    libHepMC3rootIO = "libHepMC3rootIO.dylib";
#endif
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32)) && !defined(__CYGWIN__)
    libHepMC3rootIO = "HepMC3rootIO.dll";
#endif

    std::string libHepMC3protobufIO = "libHepMC3protobufIO.so.3";
#if defined(__darwin__) || defined(__APPLE__)
    libHepMC3protobufIO = "libHepMC3protobufIO.dylib";
#endif
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32)) && !defined(__CYGWIN__)
    libHepMC3protobufIO = "HepMC3protobufIO.dll";
#endif

    bool remote = false;
    bool pipe = false;
    if (filename.find("http://") != std::string::npos)    remote = true;
    if (filename.find("https://") != std::string::npos)   remote = true;
    if (filename.find("root://") != std::string::npos)    remote = true;
    if (filename.find("gsidcap://") != std::string::npos) remote = true;

    std::vector<std::string> head;
    if (!remote)
    {
        struct stat   buffer;
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32)) && !defined(__CYGWIN__)
        if (!(stat (filename.c_str(), &buffer) == 0))
#else
        if (!(stat (filename.c_str(), &buffer) == 0 && (S_ISFIFO(buffer.st_mode) || S_ISREG(buffer.st_mode) || S_ISLNK(buffer.st_mode))))
#endif
        {
            HEPMC3_ERROR("deduce_reader: file " << filename << " does not exist or is not a regular file/FIFO/link");
            return std::shared_ptr<Reader> (nullptr);
        }

        std::ifstream* file= new std::ifstream(filename);
        if (!file)
        {
            HEPMC3_ERROR("deduce_reader could not open file for testing HepMC version: " << filename);
            return std::shared_ptr<Reader>(nullptr);
        }
        if (!file->is_open()) {
            HEPMC3_ERROR("deduce_reader could not open file for testing HepMC version: " << filename);
            file->close();
            return std::shared_ptr<Reader>(nullptr);
        }

#if defined(__linux__) || defined(__darwin__)|| defined(__APPLE__) || defined(__FreeBSD__) || defined(__sun)
        pipe = S_ISFIFO(buffer.st_mode);
        if (pipe) {
            HEPMC3_DEBUG(0, "deduce_reader: the file " << filename << " is a pipe");
            return deduce_reader(*file);
        }
#endif

        std::string line;
        size_t nonempty = 0;
        while (std::getline(*file, line) && nonempty < 3) {
            if (line.empty()) continue;
            nonempty++;
            head.push_back(line);
        }
        file->close();
        delete file;
    }
    // Assure there are at least two elements in the vector:
    head.push_back("");
    head.push_back("");
    HEPMC3_DEBUG(0, "deduce_reader: Attempt ReaderRootTree for " << filename);
    if ( strncmp(head.at(0).c_str(), "root", 4) == 0 || remote)
        return   std::make_shared<ReaderPlugin>(filename,libHepMC3rootIO,std::string("newReaderRootTreefile"));
    if (!remote)
    {
        HEPMC3_DEBUG(0, "deduce_reader: Attempt ProtobufIO for " << filename);
        if ( strncmp(head.at(0).c_str(),"hmpb",4) == 0 )
            return std::make_shared<ReaderPlugin>(filename,libHepMC3protobufIO,std::string("newReaderprotobuffile"));
#if HEPMC3_USE_COMPRESSION
        HEPMC3_DEBUG(0, "Attempt ReaderGZ for " << filename);
        char buf[6];
        snprintf(buf,6,"%s",head.at(0).c_str());
        Compression det = detect_compression_type(buf, buf + 6);
        if ( det != Compression::plaintext ) {
            HEPMC3_DEBUG(0, "Detected supported compression " << std::to_string(det));
            return deduce_reader(std::shared_ptr< std::istream >(new ifstream(filename.c_str())));
        }
#endif
        HEPMC3_DEBUG(0, "Attempt ReaderAscii for " << filename);
        if ( strncmp(head.at(0).c_str(),"HepMC::Version",14) == 0 && strncmp(head.at(1).c_str(), "HepMC::Asciiv3", 14) == 0 )
            return std::shared_ptr<Reader>((Reader*) ( new ReaderAscii(filename)));
        HEPMC3_DEBUG(0, "Attempt ReaderAsciiHepMC2 for " << filename);
        if ( strncmp(head.at(0).c_str(),"HepMC::Version",14) == 0 && strncmp(head.at(1).c_str(), "HepMC::IO_GenEvent", 18) == 0 )
            return std::shared_ptr<Reader>((Reader*) ( new ReaderAsciiHepMC2(filename)));
        HEPMC3_DEBUG(0, "Attempt ReaderLHEF for " << filename);
        if ( strncmp(head.at(0).c_str(), "<LesHouchesEvents", 17) == 0)
            return std::shared_ptr<Reader>((Reader*) ( new ReaderLHEF(filename)));
        HEPMC3_DEBUG(0, "Attempt ReaderHEPEVT for " << filename);
        std::stringstream st_e(head.at(0).c_str());
        char attr = ' ';
        bool HEPEVT = true;
        int m_i,m_p;
        while (true)
        {
            if (!(st_e >> attr)) {
                HEPEVT=false;
                break;
            }
            if (attr == ' ') continue;
            if (attr != 'E') {
                HEPEVT = false;
                break;
            }
            HEPEVT=static_cast<bool>(st_e >> m_i >> m_p);
            break;
        }
        if (HEPEVT) return std::shared_ptr<Reader>((Reader*) ( new ReaderHEPEVT(filename)));
    }
    HEPMC3_DEBUG(0, "deduce_reader: all attempts failed for " << filename);
    return std::shared_ptr<Reader>(nullptr);
}


/** @brief This function will deduce the type of input stream based on its content and will return appropriate Reader*/
std::shared_ptr<Reader> deduce_reader(std::istream &stream)
{
    std::vector<std::string> head;
    head.push_back("");
    size_t back = 0;
    size_t backnonempty = 0;
    while ( (back < 200 && backnonempty < 100) && stream) {
        char c = stream.get();
        back++;
        if (c == '\n') {
            if (head.back().length() != 0) head.push_back("");
        } else {
            head.back() += c;
            backnonempty++;
        }
    }
    if (!stream)
    {
        HEPMC3_WARNING("Input stream is too short or invalid.");
        return std::shared_ptr<Reader>(nullptr);
    }

    for (size_t i = 0; i < back; i++)  stream.unget();

    if ( strncmp(head.at(0).c_str(),"hmpb",4) == 0 )
    {
        std::string libHepMC3protobufIO = "libHepMC3protobufIO.so.3";
#if defined(__darwin__) || defined(__APPLE__)
        libHepMC3protobufIO = "libHepMC3protobufIO.dylib";
#endif
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32)) && !defined(__CYGWIN__)
        libHepMC3protobufIO = "HepMC3protobufIO.dll";
#endif

        return std::make_shared<ReaderPlugin>(stream,libHepMC3protobufIO,std::string("newReaderprotobufstream"));
    }

    if ( strncmp(head.at(0).c_str(), "HepMC::Version", 14) == 0 && strncmp(head.at(1).c_str(), "HepMC::Asciiv3", 14) == 0 )
    {
        HEPMC3_DEBUG(0, "Attempt ReaderAscii");
        return std::shared_ptr<Reader>((Reader*) ( new ReaderAscii(stream)));
    }

    if ( strncmp(head.at(0).c_str(), "HepMC::Version", 14) == 0 && strncmp(head.at(1).c_str(), "HepMC::IO_GenEvent", 18) == 0 )
    {
        HEPMC3_DEBUG(0, "Attempt ReaderAsciiHepMC2");
        return std::shared_ptr<Reader>((Reader*) ( new ReaderAsciiHepMC2(stream)));
    }

    if ( strncmp(head.at(0).c_str(), "<LesHouchesEvents", 17) == 0)
    {
        HEPMC3_DEBUG(0, "Attempt ReaderLHEF");
        return std::shared_ptr<Reader>((Reader*) ( new ReaderLHEF(stream)));
    }
    HEPMC3_DEBUG(0, "Attempt ReaderHEPEVT");
    std::stringstream st_e(head.at(0).c_str());
    char attr = ' ';
    bool HEPEVT = true;
    int m_i, m_p;
    while (true)
    {
        if (!(st_e >> attr)) {
            HEPEVT = false;
            break;
        }
        if (attr == ' ') continue;
        if (attr != 'E') {
            HEPEVT = false;
            break;
        }
        HEPEVT = static_cast<bool>(st_e >> m_i >> m_p);
        break;
    }
    if (HEPEVT) return std::shared_ptr<Reader>((Reader*) ( new ReaderHEPEVT(stream)));
    HEPMC3_DEBUG(0, "deduce_reader: all attempts failed");
    return std::shared_ptr<Reader>(nullptr);
}
/** @brief This function will deduce the type of input stream based on its content and will return appropriate Reader*/
std::shared_ptr<Reader> deduce_reader(std::shared_ptr<std::istream> stream)
{
    std::vector<std::string> head;
    head.push_back("");
    size_t back = 0;
    size_t backnonempty = 0;
    while ( (back < 200 && backnonempty < 100) && stream) {
        char c = stream->get();
        back++;
        if (c == '\n') {
            if (head.back().length() != 0) head.push_back("");
        } else {
            head.back() += c;
            backnonempty++;
        }
    }
    if (!stream)
    {
        HEPMC3_WARNING("Input stream is too short or invalid.");
        return std::shared_ptr<Reader>(nullptr);
    }

    for (size_t i = 0; i < back; i++)  stream->unget();

    if ( strncmp(head.at(0).c_str(),"hmpb",4) == 0 )
    {
        std::string libHepMC3protobufIO = "libHepMC3protobufIO.so.3";
#if defined(__darwin__) || defined(__APPLE__)
        libHepMC3protobufIO = "libHepMC3protobufIO.dylib";
#endif
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32)) && !defined(__CYGWIN__)
        libHepMC3protobufIO = "HepMC3protobufIO.dll";
#endif

        return std::make_shared<ReaderPlugin>(*stream,libHepMC3protobufIO,std::string("newReaderprotobufstream"));
    }

    if ( strncmp(head.at(0).c_str(), "HepMC::Version", 14) == 0 && strncmp(head.at(1).c_str(), "HepMC::Asciiv3", 14) == 0 )
    {
        HEPMC3_DEBUG(0, "Attempt ReaderAscii");
        return std::shared_ptr<Reader>((Reader*) ( new ReaderAscii(stream)));
    }

    if ( strncmp(head.at(0).c_str(), "HepMC::Version",14) == 0 && strncmp(head.at(1).c_str(), "HepMC::IO_GenEvent", 18) == 0 )
    {
        HEPMC3_DEBUG(0, "Attempt ReaderAsciiHepMC2");
        return std::shared_ptr<Reader>((Reader*) ( new ReaderAsciiHepMC2(stream)));
    }

    if ( strncmp(head.at(0).c_str(), "<LesHouchesEvents", 17) == 0)
    {
        HEPMC3_DEBUG(0, "Attempt ReaderLHEF");
        return std::shared_ptr<Reader>((Reader*) ( new ReaderLHEF(stream)));
    }
    HEPMC3_DEBUG(0, "Attempt ReaderHEPEVT");
    std::stringstream st_e(head.at(0).c_str());
    char attr = ' ';
    bool HEPEVT = true;
    int m_i,m_p;
    while (true)
    {
        if (!(st_e >> attr)) {
            HEPEVT = false;
            break;
        }
        if (attr == ' ') continue;
        if (attr != 'E') {
            HEPEVT = false;
            break;
        }
        HEPEVT = static_cast<bool>(st_e >> m_i >> m_p);
        break;
    }
    if (HEPEVT) return std::shared_ptr<Reader>((Reader*) ( new ReaderHEPEVT(stream)));
    HEPMC3_DEBUG(0, "deduce_reader: all attempts failed");
    return std::shared_ptr<Reader>(nullptr);
}
}
#endif
