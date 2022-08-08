// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2021 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_HEPEVT_WRAPPER_RUNTIME_H
#define HEPMC3_HEPEVT_WRAPPER_RUNTIME_H
#include <iostream>
#include <cstdio>
#include <set>
#include <map>
#include <cstring> // memset
#include <algorithm> //min max for VS2017
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/HEPEVT_Helpers.h"

/**
 *  @file HEPEVT_Wrapper_Runtime.h
 *  @brief Definition of \b class HEPEVT_Wrapper_Runtime
 *
 *  @class HepMC3::HEPEVT_Wrapper_Runtime
 *  @brief An interface to HEPEVT common block implemented to deal with varying block size in runtime
 */
namespace HepMC3
{

class HEPEVT_Wrapper_Runtime
{
//
// Functions
//
public:
    /** @brief Default constructor */
    HEPEVT_Wrapper_Runtime() {m_max_particles=0; m_hepevtptr=nullptr;};
    /** @brief Default destructor */
    ~HEPEVT_Wrapper_Runtime() {};
    /** @brief Print information from HEPEVT common block */
    void print_hepevt( std::ostream& ostr = std::cout ) const;
    /** @brief Print particle information */
    void print_hepevt_particle( int index, std::ostream& ostr = std::cout ) const;
    /** @brief Set all entries in HEPEVT to zero */
    void zero_everything();
    /** @brief Convert GenEvent to HEPEVT*/
    bool GenEvent_to_HEPEVT( const GenEvent* evt ) { return GenEvent_to_HEPEVT_nonstatic(evt, this);};
    /** @brief Convert HEPEVT to GenEvent*/
    bool HEPEVT_to_GenEvent( GenEvent* evt ) const { return HEPEVT_to_GenEvent_nonstatic(evt, this);};
    /** @brief Tries to fix list of daughters */
    bool fix_daughters();
private:
    /** @brief  Fortran common block HEPEVT */
    std::shared_ptr<struct HEPEVT_Pointers<double> >  m_hepevtptr;
    /** @brief Block size */
    int m_max_particles;
    /** @brief  Internalstorage storage. Optional.*/
    std::vector<char> m_internal_storage;
//
// Accessors
//
public:
    void   allocate_internal_storage(); //!< Allocates m_internal_storage storage in smart pointer to hold HEPEVT of fixed size
    void   copy_to_internal_storage( char *c, int N ); //!< Copies the content of foreign common block into the internal storage
    void   set_max_number_entries( unsigned int size ) { m_max_particles = size; }//!< Set block size
    void   set_hepevt_address(char *c); //!< Set Fortran block address
    int    max_number_entries()  const     { return m_max_particles;                              } //!< Block size
    int    event_number()     const        { return *(m_hepevtptr->nevhep);             } //!< Get event number
    int    number_entries()  const         { return *(m_hepevtptr->nhep);               } //!< Get number of entries
    int    status(const int index )   const     { return m_hepevtptr->isthep[index-1];    } //!< Get status code
    int    id(const int index )     const       { return m_hepevtptr->idhep[index-1];     } //!< Get PDG particle id
    int    first_parent(const int index ) const { return m_hepevtptr->jmohep[2*(index-1)]; } //!< Get index of 1st mother
    int    last_parent(const int index )  const { return m_hepevtptr->jmohep[2*(index-1)+1]; } //!< Get index of last mother
    int    first_child(const int index )  const { return m_hepevtptr->jdahep[2*(index-1)]; } //!< Get index of 1st daughter
    int    last_child(const int index )  const  { return m_hepevtptr->jdahep[2*(index-1)+1]; } //!< Get index of last daughter
    double px(const int index ) const           { return m_hepevtptr->phep[5*(index-1)];   } //!< Get X momentum
    double py(const int index )  const          { return m_hepevtptr->phep[5*(index-1)+1];   } //!< Get Y momentum
    double pz(const int index ) const           { return m_hepevtptr->phep[5*(index-1)+2];   } //!< Get Z momentum
    double e(const int index )  const           { return m_hepevtptr->phep[5*(index-1)+3];   } //!< Get Energy
    double m(const int index )  const           { return m_hepevtptr->phep[5*(index-1)+4];   } //!< Get generated mass
    double x(const int index )  const           { return m_hepevtptr->vhep[4*(index-1)];   } //!< Get X Production vertex
    double y(const int index )  const           { return m_hepevtptr->vhep[4*(index-1)+1];   } //!< Get Y Production vertex
    double z(const int index )   const          { return m_hepevtptr->vhep[4*(index-1)+2];   } //!< Get Z Production vertex
    double t(const int index )   const          { return m_hepevtptr->vhep[4*(index-1)+3];   } //!< Get production time
    int    number_parents(const int index ) const;                                   //!< Get number of parents
    int    number_children(const int index ) const;                                  //!< Get number of children from the range of daughters
    int    number_children_exact(const int index ) const;                            //!< Get number of children by counting
    void set_event_number( const int evtno )       { *(m_hepevtptr->nevhep) = evtno;         } //!< Set event number
    void set_number_entries( const int noentries ) { *(m_hepevtptr->nhep) = noentries;       } //!< Set number of entries
    void set_status( const int index, const int status ) { m_hepevtptr->isthep[index-1] = status; } //!< Set status code
    void set_id(const int index, const int id )         { m_hepevtptr->idhep[index-1] = id;      } //!< Set PDG particle id
    void set_parents( const int index, const int firstparent, const int lastparent );              //!< Set parents
    void set_children( const int index, const int firstchild, const int lastchild );               //!< Set children
    void set_momentum( const int index, const double px, const double py, const double pz, const double e );   //!< Set 4-momentum
    void set_mass( const int index, double mass );                                     //!< Set mass
    void set_position( const int index, const double x, const double y, const double z, const double t );      //!< Set position in time-space
};


} // namespace HepMC3
#endif
