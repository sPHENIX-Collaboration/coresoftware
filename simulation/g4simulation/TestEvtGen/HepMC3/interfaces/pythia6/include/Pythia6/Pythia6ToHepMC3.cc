// -*- C++ -*-
//
// This file is part of HepMC3
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef Pythia6_Pythia6ToHepMC3_H
#define Pythia6_Pythia6ToHepMC3_H
#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32)) && !defined(__CYGWIN__)
#define  hepmc3_delete_writer_ HEPMC3_DELETE_WRITER
#define  hepmc3_convert_event_ HEPMC3_CONVERT_EVENT
#define  hepmc3_clear_event_ HEPMC3_CLEAR_EVENT
#define  hepmc3_write_event_ HEPMC3_WRITE_EVENT
#define  hepmc3_set_cross_section_ HEPMC3_SET_CROSS_SECTION
#define  hepmc3_set_pdf_info_ HEPMC3_SET_PDF_INFO
#define  hepmc3_set_event_number_ HEPMC3_SET_EVENT_NUMBER
#define  hepmc3_set_hepevt_address_ HEPMC3_SET_HEPEVT_ADDRESS
#define  hepmc3_set_attribute_int_ HEPMC3_SET_ATTRIBUTE_INT
#define  hepmc3_set_attribute_double_ HEPMC3_SET_ATTRIBUTE_DOUBLE
#define  hepmc3_new_writer_ HEPMC3_NEW_WRITER
#define  hepmc3_new_weight_ HEPMC3_NEW_WEIGHT
#define  hepmc3_set_weight_by_index_ HEPMC3_SET_WEIGHT_BY_INDEX
#define  hepmc3_set_weight_by_name_ HEPMC3_SET_WEIGHT_BY_NAME
#endif
#ifdef  DUMMYPYTHIA6TOHEPMC3
extern "C" {

    int hepmc3_delete_writer_(const int & position)
    {
        return -1;
    }
    int hepmc3_convert_event_(const int & position)
    {
        return -1;
    }
    int hepmc3_write_event_(const int & position)
    {
        return -1;
    }
    int hepmc3_clear_event_(const int & position)
    {
        return -1;
    }
    int hepmc3_set_cross_section_(const int & position, const double& x, const double& xe, const int& n1, const int& n2)
    {
        return -1;
    }

    int hepmc3_set_pdf_info_(const int & position, const int& parton_id1, const int& parton_id2, const double& x1, const double& x2,
                             const double& scale_in, const double& xf1, const double& xf2,
                             const int& pdf_id1, const int& pdf_id2)
    {
        return -1;
    }
    int hepmc3_set_hepevt_address_(int* a)
    {
        return -1;
    }
    int hepmc3_set_attribute_int_(const int & position, const int & attval, const char* attname)
    {
        return -1;
    }
    int hepmc3_set_attribute_double_(const int & position, const double & attval, const char* attname)
    {
        return -1;
    }
    int hepmc3_new_writer_(const int & position, const int & mode, const char* ffilename)
    {
        return  -1;
    }
    int hepmc3_new_weight_(const int & position, const char* name)
    {
        return -1;
    }
    int hepmc3_set_weight_by_index_(const int & position, const double& val, const int & pos)
    {
        return -1;
    }
    int hepmc3_set_weight_by_name_(const int & position, const double& val, const char* name)
    {
        return -1;
    }
}


#else
#include "HepMC3/HEPEVT_Wrapper_Template.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/Writer.h"
#include "HepMC3/WriterHEPEVT.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/WriterPlugin.h"
#include "HepMC3/Print.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
using namespace HepMC3;
#ifndef PYTHIA6HEPEVTSIZE
#define PYTHIA6HEPEVTSIZE 10000
#endif
/** The conversion interface, templated version */
HEPEVT_Wrapper_Template<PYTHIA6HEPEVTSIZE> hepmc3_gInterface;
/** Storage for the output objects (Writers)*/
std::map<int, std::pair<std::shared_ptr<Writer>, GenEvent*> > hepmc3_gWriters;
/** Storage for the GenRunInfo objects associated with the outputs */
std::map<int, std::shared_ptr<GenRunInfo> >  hepmc3_gGenRunInfos;
/** Interface to acces the enets from C++, e.g. Rivet */
GenEvent* hepmc3_gWriters_get_event(const int & position)
{
    if (hepmc3_gWriters.count(position) == 0) {
        printf("Warning in %s: Writer at position %i does not exist\n", __FUNCTION__, position);
        return NULL;
    }
    return    hepmc3_gWriters[position].second;
}
/** Interfaces for C/Fortran */
extern "C" {

    int hepmc3_delete_writer_(const int & position)
    {
        if (hepmc3_gWriters.count(position) == 0) {
            printf("Warning in %s: Writer at position %i does not exist\n", __FUNCTION__, position);
            return 1;
        }
        hepmc3_gWriters[position].first->close();
        hepmc3_gWriters.erase(hepmc3_gWriters.find(position));
        return 0;

    }
    int hepmc3_convert_event_(const int & position)
    {
        if (hepmc3_gWriters.count(position) == 0) {
            printf("Warning in %s: Writer at position %i does not exist\n", __FUNCTION__, position);
            return 1;
        }
        if (!hepmc3_gInterface.m_hepevtptr)
        {
            printf("Error in %s: HEPEVT block does not exist\n", __FUNCTION__);
            return 1;
        }
        hepmc3_gWriters[position].second = new GenEvent(Units::GEV, Units::MM);
        for( int i = 1; i <= hepmc3_gInterface.number_entries(); i++ )
            if (hepmc3_gInterface.m_hepevtptr->jmohep[i-1][1]<hepmc3_gInterface.m_hepevtptr->jmohep[i-1][0])  hepmc3_gInterface.m_hepevtptr->jmohep[i-1][1] = hepmc3_gInterface.m_hepevtptr->jmohep[i-1][0];
        hepmc3_gInterface.HEPEVT_to_GenEvent(hepmc3_gWriters[position].second);
        if (hepmc3_gGenRunInfos.count(position)  ==  0) hepmc3_gGenRunInfos[position] = std::make_shared<GenRunInfo>();
        hepmc3_gWriters[position].second->set_run_info(hepmc3_gGenRunInfos[position]);
        hepmc3_gWriters[position].second->weights() = std::vector<double>(hepmc3_gGenRunInfos[position]->weight_names().size(), 1.0);
        return 0;
    }
    int hepmc3_write_event_(const int & position)
    {
        if (hepmc3_gWriters.count(position) == 0) {
            printf("Warning in %s: Writer at position %i does not exist\n", __FUNCTION__, position);
            return 1;
        }
        hepmc3_gWriters[position].first->write_event(*(hepmc3_gWriters[position].second));
        return 0;
    }
    int hepmc3_clear_event_(const int & position)
    {
        if (hepmc3_gWriters.count(position) == 0) {
            printf("Warning in %s: Writer at position %i does not exist\n", __FUNCTION__, position);
            return 1;
        }
        hepmc3_gWriters[position].second->clear();
        return 0;
    }
    int hepmc3_set_cross_section_(const int & position, const double& x, const double& xe, const int& n1, const int& n2)
    {
        if (hepmc3_gWriters.count(position) == 0) {
            printf("Warning in %s: Writer at position %i does not exist\n", __FUNCTION__, position);
            return 1;
        }
        GenCrossSectionPtr cs = std::make_shared< GenCrossSection>();
        cs->set_cross_section(x, xe, n1, n2);
        hepmc3_gWriters[position].second->set_cross_section(cs);
        return 0;
    }

    int hepmc3_set_pdf_info_(const int & position, const int& parton_id1, const int& parton_id2, const double& x1, const double& x2,
                             const double& scale_in, const double& xf1, const double& xf2,
                             const int& pdf_id1, const int& pdf_id2)
    {
        if (hepmc3_gWriters.count(position) == 0) {
            printf("Warning in %s: Writer at position %i does not exist\n", __FUNCTION__, position);
            return 1;
        }
        GenPdfInfoPtr pdf=std::make_shared< GenPdfInfo>();
        pdf->set(parton_id1, parton_id2, x1, x2, scale_in, xf1, xf2, pdf_id1, pdf_id2);
        hepmc3_gWriters[position].second->set_pdf_info(pdf);
        return 0;
    }
    int hepmc3_set_hepevt_address_(int* a)
    {
      printf("Info in %s: setting /hepevt/ block adress\n", __FUNCTION__);
      hepmc3_gInterface.set_hepevt_address((char*)a);
      return 0;
    }
    int hepmc3_set_attribute_int_(const int & position, const int & attval, const char* attname)
    {
        if (hepmc3_gWriters.count(position) == 0) {
            printf("Warning in %s: Writer at position %i does not exist\n", __FUNCTION__, position);
            return 1;
        }
        hepmc3_gWriters[position].second->add_attribute(attname, std::make_shared<IntAttribute>(attval));
        return 0;
    }
    int hepmc3_set_attribute_double_(const int & position, const double & attval, const char* attname)
    {
        if (hepmc3_gWriters.count(position) == 0) {
            printf("Warning in %s: Writer at position %i does not exist\n", __FUNCTION__, position);
            return 1;
        }
        hepmc3_gWriters[position].second->add_attribute(attname, std::make_shared<DoubleAttribute>(attval));
        return 0;
    }

    int hepmc3_new_writer_(const int & position, const int & mode, const char* ffilename)
    {
        std::string libHepMC3rootIO="libHepMC3rootIO.so";
#ifdef __darwin__
        libHepMC3rootIO="libHepMC3rootIO.dylib";
#endif
#ifdef WIN32
        libHepMC3rootIO="HepMC3rootIO.dll";
#endif
        std::string filename=std::string(ffilename);
        int r_position=position;
        if (r_position == 0)
        {
            if (hepmc3_gWriters.size() == 0) r_position = 1;
            if (hepmc3_gWriters.size()  !=  0) r_position = hepmc3_gWriters.rend()->first+1;
        }
        if (hepmc3_gWriters.count(r_position) != 0) {
            printf("Error in %s: Writer at position %i already exists\n", __FUNCTION__, r_position);
            exit(1);
        }
        if (hepmc3_gGenRunInfos.count(r_position) != 0) {
            printf("Warning in %s: RunInfo at position %i already exists\n", __FUNCTION__, r_position);
        }
        else
        {
            hepmc3_gGenRunInfos[r_position]=std::make_shared<GenRunInfo>();
        }

        switch (mode)
        {
        case 1:
            hepmc3_gWriters[r_position] = std::pair<std::shared_ptr<Writer>, GenEvent*>(std::make_shared<WriterAscii>(filename.c_str(), hepmc3_gGenRunInfos[position]), new GenEvent(hepmc3_gGenRunInfos[position], Units::GEV, Units::MM));
            break;
        case 2:
            hepmc3_gWriters[r_position] = std::pair<std::shared_ptr<Writer>, GenEvent*>(std::make_shared<WriterAsciiHepMC2>(filename.c_str(), hepmc3_gGenRunInfos[position]), new GenEvent(hepmc3_gGenRunInfos[position], Units::GEV, Units::MM));
            break;
        case 3:
            hepmc3_gWriters[r_position] = std::pair<std::shared_ptr<Writer>, GenEvent*>(std::make_shared<WriterHEPEVT>(filename.c_str()), new GenEvent(hepmc3_gGenRunInfos[position], Units::GEV, Units::MM));
            break;
        case 4:
            hepmc3_gWriters[r_position] = std::pair<std::shared_ptr<Writer>, GenEvent*>(std::make_shared<WriterPlugin>(filename.c_str(), libHepMC3rootIO, std::string("newWriterRootfile"), hepmc3_gGenRunInfos[position]), new GenEvent(hepmc3_gGenRunInfos[position], Units::GEV, Units::MM));
            break;
        case 5:
            hepmc3_gWriters[r_position] = std::pair<std::shared_ptr<Writer>, GenEvent*>(std::make_shared<WriterPlugin>(filename.c_str(), libHepMC3rootIO, std::string("newWriterRootTreefile"), hepmc3_gGenRunInfos[position]), new GenEvent(hepmc3_gGenRunInfos[position], Units::GEV, Units::MM));
            break;
        default:
            printf("Error in %s:Output format %d is unknown or not supported.\n", __FUNCTION__, mode);
            exit(2);
            break;
        }
        return  r_position;
    }
    int hepmc3_new_weight_(const int & position, const char* name)
    {
        if (hepmc3_gGenRunInfos.count(position) == 0) {
            printf("Warning in %s: RunInfo at position %i does not exist\n", __FUNCTION__, position);
            return 1;
        }
        if (hepmc3_gGenRunInfos[position]->weight_index(std::string(name)) >= 0) return 0;
        std::vector<std::string> weight_names = hepmc3_gGenRunInfos[position]->weight_names();
        weight_names.push_back(std::string(name));
        hepmc3_gWriters[position].second->weights().push_back(1.0);
        hepmc3_gGenRunInfos[position]->set_weight_names(weight_names);
        return 0;
    }
    int hepmc3_set_weight_by_index_(const int & position, const double& val, const int & index)
    {
        if (hepmc3_gWriters.count(position) == 0) {
            printf("Warning in %s: Writer at position %i does not exist\n", __FUNCTION__, position);
            return 1;
        }
        if (hepmc3_gWriters[position].second->weights().size() < index) {
            printf("Warning in %s: Event has no weight with index %i\n", __FUNCTION__, index);
            return 2;
        }
        hepmc3_gWriters[position].second->weights()[index] = val;
        return 0;
    }
    int hepmc3_set_weight_by_name_(const int & position, const double& val, const char* name)
    {
        if (hepmc3_gWriters.count(position) == 0) {
            printf("Warning in %s: Writer at position %i does not exist\n", __FUNCTION__, position);
            return 1;
        }
        hepmc3_new_weight_(position, name);
        hepmc3_gWriters[position].second->weight(std::string(name)) = val;
        return 0;
    }
}
#endif
#endif
