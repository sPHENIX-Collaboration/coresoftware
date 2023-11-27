#ifndef JETBASE_JetContainerv1__h
#define JETBASE_JetContainerv1__h
#include "JetContainer.h"

#include "Jet.h"

class JetContainerv1 : public JetContainer
{
public:
    // Manage class creation, copy, reset, destruction
    JetContainerv1();
    ~JetContainerv1() override;
    void identify(std::ostream& os=std::cout) const override;
    explicit JetContainerv1(const JetContainer& jets);
    JetContainerv1& operator=(const JetContainer& jets);
    void Reset() override;
    TClonesArray* clone_data() const override 
    { return (TClonesArray*) m_clones->Clone(); };

    // status of jet contents
    bool   empty() const override { return m_njets==0; };
    size_t size()  const override { return m_njets; };

    // adding/access jets
    Jet* add_jet()                   override; // Add a new jet to the TClonesArray and return the pointer
    Jet* get_jet(unsigned int index) override; // Get get at location. 
    Jet* get_UncheckedAt(unsigned int index) override; // Get get at location.
                                                
    // convenience shortcuts of get_{jet,UncheckedAt}
    inline Jet* operator()(int index) override { return get_jet(index); }; // synonym for get_jet()
    inline Jet* operator[](int index) override { return get_UncheckedAt(index); }; // get jet, don't check for length
                                                                                     //
    // ----------------------------------------------------------------------------------------
    //  Interface for adding properties to jets, and getting the index to those properties
    // ----------------------------------------------------------------------------------------
    // The optional properties in each Jet are stored in a vector of floats,
    // e.g. { Rg_value, mg_value, area, .... }
    // The JetContainer generates the jet and keeps a map of which properties in which location.
    // Default values in the Jet Properties are set to NAN
    // ----------------------------------------------------------------------------------------
    std::map<Jet::PROPERTY, Jet::PROPERTY> property_indices() const override { return m_pindex; };
    bool has_property(Jet::PROPERTY prop) const override { return m_pindex.find(prop) != m_pindex.end(); };
    size_t size_properties() const override { return m_pindex.size(); };
    size_t add_property(Jet::PROPERTY) override; // add property, if not already there, to jet prop vectors
    size_t add_property(std::set<Jet::PROPERTY>)    override; // same as above, convenience for other code using ::set
    virtual Jet::PROPERTY property_index(Jet::PROPERTY) override; // get the propery index
    void print_property_types(std::ostream& os = std::cout) const override; // print the order of properties in jet

    
    // ---------------------------------------------------------------------------------------
    //  Add ability to loop over jets
    // ---------------------------------------------------------------------------------------
    // Use:
    // ```
    //     for (jet : *jet_containter) { 
    //     // jet will increment as Jet* through all jets
    //     ... }
    // ```
    Jet::IterJetTCA begin() override;
    Jet::IterJetTCA end()   override;
    // ---------------------------------------------------------------------------------------

    // -legacy-set-parameters-----------------------------------------------------------------
    void set_algo(Jet::ALGO algo) override { m_algo = algo; };
    Jet::ALGO get_algo() const override { return m_algo; };

    void  set_par (float par) override { set_jetpar_R(par); }
    float get_par () const override { return get_jetpar_R(); }

    void  set_jetpar_R (float par) override { m_jetpar_R = par; }
    float get_jetpar_R () const override { return m_jetpar_R; }
    
    // set access to source identifiers ------------------------------------------

    bool empty_src() const override { return m_src.empty(); }
    void insert_src(Jet::SRC src) override { m_src.insert(src); }

    ConstSrcIter begin_src() const override { return m_src.begin(); }
    ConstSrcIter find_src(Jet::SRC src) const override { return m_src.find(src); }
    ConstSrcIter end_src() const override { return m_src.end(); }

    SrcIter begin_src() override { return m_src.begin(); }
    SrcIter find_src(Jet::SRC src) override { return m_src.find(src); }
    SrcIter end_src() override { return m_src.end(); }

    void print_jets(std::ostream&) override; // print the order of properties in jet

    void set_rho_median(float _) override { m_RhoMedian = _; };
    float get_rho_median() const override { return m_RhoMedian; };

private:
    std::string str_Jet_PROPERTY(Jet::PROPERTY) const;

    TClonesArray* m_clones {nullptr}; // TClonesArray of Jet objects
    size_t m_njets {0}; // size of jet_array

    // properties contained in vectors of all jets
    std::map<Jet::PROPERTY, Jet::PROPERTY /*really index for vectors in jets*/> m_pindex {}; // indices of properties in each jet property vector
    size_t m_psize {0}; // size of p_index and p_vec

    void resize_jet_pvecs();

    // status
    Jet::ALGO m_algo {Jet::NONE};
    float m_jetpar_R {0.4}; 

    std::set<Jet::SRC> m_src;      //< set of sources (clusters, towers, etc)
    
    float m_RhoMedian {NAN};

    ClassDefOverride(JetContainerv1, 1);
};

#endif
