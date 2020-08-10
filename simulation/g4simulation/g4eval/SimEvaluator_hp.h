#ifndef G4EVAL_SIMEVALUATOR_HP_H
#define G4EVAL_SIMEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>

#include <map>
#include <set>
#include <vector>

class EventHeader;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;

class SimEvaluator_hp : public SubsysReco
{
  public:

  /// constructor
  SimEvaluator_hp( const std::string& = "SimEvaluator_hp" );

  /// global initialization
  virtual int Init(PHCompositeNode*);

  /// run initialization
  virtual int InitRun(PHCompositeNode*);

  /// event processing
  virtual int process_event(PHCompositeNode*);

  /// end of processing
  virtual int End(PHCompositeNode*);

  // event information
  class EventStruct
  {

    public:

    using List = std::vector<EventStruct>;

    //! impact parameter, when relevant
    float _bimp = 0;

    //! reaction plane angle, when relevant
    float _rplane = 0;

    //!@name used keep track of the number of pileup events
    //@{
    int _nevt = 0;
    int _nevt_active = 0;
    int _nevt_bg = 0;
    //@}
  };

  // vertex information
  class VertexStruct
  {
    public:
    using List = std::vector<VertexStruct>;
    int _embed = 0;
    float _x = 0;
    float _y = 0;
    float _z = 0;
    float _t = 0;
    bool _is_main_vertex = false;
  };

  // particle information
  class ParticleStruct
  {
    public:
    using List = std::vector<ParticleStruct>;

    int _charge = 0;
    int _pid = 0;
    int _embed = 0;
    bool _is_primary = false;
    int64_t _mask = 0;

    float _px = 0;
    float _py = 0;
    float _pz = 0;
    float _pt = 0;
    float _p = 0;
    float _eta = 0;

  };

  // particle information
  class G4HitStruct
  {
    public:
    using List = std::vector<G4HitStruct>;
    int _embed = 0;
    int _detid = 0;
    float _x = 0;
    float _y = 0;
    float _z = 0;
    float _t = 0;
  };

  class Container: public PHObject
  {

    public:

    /// constructor
    explicit Container() = default;

    /// copy constructor
    explicit Container(const Container &) = delete;

    /// assignment operator
    Container& operator = ( const Container& ) = delete;

    /// reset
    virtual void Reset();

    ///@name accessors
    //@{

    const EventStruct::List& events() const
    { return _events; }

    const VertexStruct::List& vertexList() const
    { return _vertex_list; }

    const ParticleStruct::List& particleList() const
    { return _particle_list; }

    const G4HitStruct::List& g4hits() const
    { return _g4hits; }

    //@}

    ///@name modifiers
    //@{

    void addEvent( const EventStruct& event )
    { _events.push_back( event ); }

    void addVertex( const VertexStruct& vertex )
    { _vertex_list.push_back( vertex ); }

    void addParticle( const ParticleStruct& particle )
    { _particle_list.push_back( particle ); }

    void addG4Hit( const G4HitStruct& hit )
    { _g4hits.push_back( hit ); }

    void clearEventList()
    { _events.clear(); }

    void clearVertexList()
    { _vertex_list.clear(); }

    void clearParticleList()
    { _particle_list.clear(); }

    void clearG4Hits()
    { _g4hits.clear(); }

    //@}

    private:

    /// event struct
    EventStruct::List _events;

    //* vertex list
    VertexStruct::List _vertex_list;

    //* particles
    ParticleStruct::List _particle_list;

    //* hits
    G4HitStruct::List _g4hits;

    ClassDef(Container,1)

  };

  enum Flags
  {
    EvalEvent = 1<<0,
    EvalVertices = 1<<1,
    EvalParticles = 1<<2,
    EvalHits = 1<<3,
    PrintVertices = 1<<4
  };

  /// set flags. Should be a bitwise or of Flags enum
  void set_flags( int flags )
  { m_flags = flags; }

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// print TPC geometry
  void print_tpc( PHCompositeNode* );

  /// fill MC track map
  void fill_g4particle_map();

  /// genevent
  void check_genevent();

  /// fill event struct
  void fill_event();

  /// fill vertices
  void fill_vertices();

  /// fill particles
  void fill_particles();

  /// fill hits
  void fill_hits();

  /// print vertices
  void print_vertices();

  // get embedded id for given g4hit
  int get_embed(PHG4Hit*);

  // get embedded id for given g4track
  int get_embed(PHG4Particle*) const;

  //* data container
  Container* m_container = nullptr;

  // flags
  // int m_flags = EvalEvent | EvalVertices | EvalParticles | EvalHits;
  int m_flags = EvalEvent | EvalVertices | EvalParticles;

  PHG4HitContainer* m_g4hits_mvtx = nullptr;
  PHG4HitContainer* m_g4hits_intt = nullptr;
  PHG4HitContainer* m_g4hits_tpc = nullptr;
  PHG4HitContainer* m_g4hits_micromegas = nullptr;

  //* hep event
  PHHepMCGenEventMap* m_geneventmap = nullptr;

  //* truth information
  PHG4TruthInfoContainer* m_g4truthinfo = nullptr;

  //! event header
  EventHeader* m_eventheader = nullptr;

  // map trk id to embed id
  using EmbedMap = std::map<int,int>;
  EmbedMap m_g4embed_map;
  
  // map trk_id to layer mask
  using G4ParticleMap = std::map<int,int64_t>;
  G4ParticleMap m_g4particle_map;

};

#endif  // G4EVAL_SIMEVALUATOR_HP_H
