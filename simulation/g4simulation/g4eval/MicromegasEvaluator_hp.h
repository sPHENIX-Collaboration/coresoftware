#ifndef G4EVAL_MICROMEGASEVALUATOR_HP_H
#define G4EVAL_MICROMEGASEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>
#include <gsl/gsl_rng.h>
#include <micromegas/MicromegasDefs.h>
#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <map>
#include <memory>

class PHG4CylinderGeomContainer;
class PHG4HitContainer;
class TrkrHitSetContainer;
class TrkrHit;

class MicromegasEvaluator_hp : public SubsysReco
{
  public:

  //! constructor
  MicromegasEvaluator_hp( const std::string& = "MicromegasEvaluator_hp" );

  //! global initialization
  virtual int Init(PHCompositeNode*);

  //! run initialization
  virtual int InitRun(PHCompositeNode*);

  //! event processing
  virtual int process_event(PHCompositeNode*);

  //! end of processing
  virtual int End(PHCompositeNode*);

  // tile information
  class TileStruct
  {
    public:

    using List = std::vector<TileStruct>;

    //! default constructor
    TileStruct() = default;

    //! constructor
    TileStruct( int layer, int tile ):
      _layer( layer ),
      _tile( tile )
    {}

    //! layer
    uint _layer = 0;

    //! tile
    uint _tile = 0;

    //! total deposited energy in tile
    float _edep_total = 0;

    //! total ionization energy in tile
    float _eion_total = 0;

    //! total number of electron (after amplification)
    float _electrons_total = 0;

    //! total number of fired strips
    int _strips_total = 0;
  };

  // g4hit information to be stored in the tree
  class G4HitStruct
  {
    public:

    using List = std::vector<G4HitStruct>;

    //! layer
    uint _layer = 0;

    //! tile
    uint _tile = 0;

    // position
    float _x = 0;
    float _y = 0;
    float _z = 0;
    float _t = 0;

    //! deposited energy
    float _edep = 0;

    //! ionization energy
    float _eion = 0;

    //! number of primary electrons
    uint _nprimary = 0;

    //! total number of electrons
    uint _nelectron = 0;

  };

  //! hit information to be stored in tree
  class HitStruct
  {
    public:

    using List = std::vector<HitStruct>;

    //! layer
    uint _layer = 0;

    //! number of hits belonging to the cluster
    uint _tile = 0;

    //! strip
    uint _strip = 0;

    //! energy
    float _energy = 0;

    //! ADC counts
    uint _adc = 0;

  };

  //! track container
  class Container: public PHObject
  {

    public:

    //! constructor
    explicit Container() = default;

    //! copy constructor
    explicit Container(const Container &) = delete;

    //! assignment operator
    Container& operator = ( const Container& ) = delete;

    //! reset
    virtual void Reset();

    //!@name accessors
    //@{

    const TileStruct::List& tiles() const
    { return _tiles; }

    const G4HitStruct::List& g4hits() const
    { return _g4hits; }

    const HitStruct::List& hits() const
    { return _hits; }

    //@}

    //!@name modifiers
    //@{

    //! find tile matching given layer and tile id
    /* create and store new one if not found */
    TileStruct& findTile( const int layer, const int tile );

    void addG4Hit( const G4HitStruct& hit )
    { _g4hits.push_back( hit ); }

    void addHit( const HitStruct& hit )
    { _hits.push_back( hit ); }

    void clearTiles()
    { _tiles.clear(); }

    void clearHits()
    { _hits.clear(); }

    void clearG4Hits()
    { _g4hits.clear(); }

    //@}

    private:

    //! tiles
    TileStruct::List _tiles;

    //! g4hits array
    G4HitStruct::List _g4hits;

    //! hits array
    HitStruct::List _hits;

    ClassDef(Container,1)

  };

  enum Flags
  {
    EvalG4Hits = 1<<0,
    EvalHits = 1<<1
  };

  //! set flags. Should be a bitwise or of Flags enum
  void set_flags( int flags )
  { m_flags = flags; }

  private:

  //! load nodes
  int load_nodes( PHCompositeNode* );

  //! evaluate g4hits
  void evaluate_g4hits();

  //! evaluate hits
  void evaluate_hits();

  //! cluster array
  Container* m_container = nullptr;

  // flags
  int m_flags = EvalG4Hits | EvalHits;

  //! gemometry
  PHG4CylinderGeomContainer* m_geonode = nullptr;

  //! g4hit container
  PHG4HitContainer* m_g4hits_micromegas = nullptr;

  //! hitset container
  TrkrHitSetContainer* m_hitsetcontainer = nullptr;

  //! rng de-allocator
  class Deleter
  {
    public:
    //! deletion operator
    void operator() (gsl_rng* rng) const { gsl_rng_free(rng); }
  };

  //! random generator that conform with sPHENIX standard
  /*! using a unique_ptr with custom Deleter ensures that the structure is properly freed when parent object is destroyed */
  std::unique_ptr<gsl_rng, Deleter> m_rng;

};

#endif  // G4EVAL_MicromegasEvaluator_HP_H
