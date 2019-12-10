
#ifndef HijJet4_h
#define HijJet4_h

extern "C" { void* hijjet4_address_(void); }

/**
@class HijJet4
@brief  Class definition for HijJet4, which is used
     to modify the Hijing HIJJET4 common.
*/

class HijJet4 {
public:
    HijJet4();
    ~HijJet4();
    
    int&	ndr	(void);
    int&    	iadr	(int i, int j);
    int&    	kfdr	(int i);
    float&  	pdr	(int i, int k);
    
    void	init	(void);

    // return common array lengths
    inline int	lenI() const {return _lenI;}
    inline int	lenJ() const {return _lenJ;}
    inline int	lenK() const {return _lenK;}

private: 

    // Lengths of array in HiMain2 common
    static const int _lenI	= 900;
    static const int _lenJ	= 2;
    static const int _lenK	= 5;

    struct HIJJET4;
    friend struct HIJJET4;

    struct HIJJET4
    {
	int    	ndr;
	int	iadr	[_lenJ][_lenI];
	int    	kfdr	[_lenI];
	float  	pdr	[_lenK][_lenI];
    };

    int  _dummy;
    float  _realdummy;

    static HIJJET4* _hijjet4;
};

// set pointer to zero at start
HijJet4::HIJJET4* HijJet4::_hijjet4 =0;

inline void
HijJet4::init(void)
{ if (!_hijjet4) _hijjet4 = static_cast<HIJJET4*>(hijjet4_address_()); }

// Constructor
inline
HijJet4::HijJet4()
    : _dummy		(-999),
      _realdummy	(-999.)
{}

// Destructor
inline
HijJet4::~HijJet4()
{}

inline int&
HijJet4::ndr	(void)
{
    init(); // check COMMON is initialized
    return _hijjet4->ndr;
}

inline int&
HijJet4::iadr	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _dummy;

    return _hijjet4->iadr[j-1][i-1];
}

inline int&
HijJet4::kfdr    (int i)
{
    init(); // check COMMON is initialized
    if(i < 1 || i > lenI()) return _dummy;
    return _hijjet4->kfdr[i-1];
}

inline float&
HijJet4::pdr	(int i, int k)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	k < 1 || k > lenK() ) return _realdummy;

    return _hijjet4->pdr[k-1][i-1];
}

#endif
