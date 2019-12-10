
#ifndef HijJet2_h
#define HijJet2_h

extern "C" { void* hijjet2_address_(void); }
/**
@class HijJet2

@brief Class definition for HijJet2, which is used
      to modify the Hijing HIJJET2 common.
*/
class HijJet2 {
public:
    HijJet2();
    ~HijJet2();
    
    int&	nsg	(void);
    int&    	njsg	(int i);
    int&    	iasg	(int i, int k);
    int&    	k1sg	(int i, int j);
    int&    	k2sg	(int i, int j);
    float&  	pxsg	(int i, int j);
    float&  	pysg	(int i, int j);
    float&  	pzsg	(int i, int j);
    float&  	pesg	(int i, int j);
    float&  	pmsg	(int i, int j);
    
    void	init	(void);

    // return common array lengths
    inline int	lenI() const {return _lenI;}
    inline int	lenJ() const {return _lenJ;}
    inline int	lenK() const {return _lenK;}

private: 

    // Lengths of array in HiMain2 common
    static const int _lenI	= 900;
    static const int _lenJ	= 100;
    static const int _lenK	= 3;

    struct HIJJET2;
    friend struct HIJJET2;

    struct HIJJET2
    {
	int    	nsg;
	int    	njsg	[_lenI];
	int	iasg	[_lenK][_lenI];
	int	k1sg	[_lenJ][_lenI];
	int	k2sg	[_lenJ][_lenI];
	float  	pxsg	[_lenJ][_lenI];
	float  	pysg	[_lenJ][_lenI];
	float  	pzsg	[_lenJ][_lenI];
	float  	pesg	[_lenJ][_lenI];
	float  	pmsg	[_lenJ][_lenI];
    };

    int  _dummy;
    float  _realdummy;

    static HIJJET2* _hijjet2;
};

// set pointer to zero at start
HijJet2::HIJJET2* HijJet2::_hijjet2 =0;

inline void
HijJet2::init(void)
{ if (!_hijjet2) _hijjet2 = static_cast<HIJJET2*>(hijjet2_address_()); }

// Constructor
inline
HijJet2::HijJet2()
    : _dummy		(-999),
      _realdummy	(-999.)
{}

// Destructor
inline
HijJet2::~HijJet2()
{}

inline int&
HijJet2::nsg	(void)
{
    init(); // check COMMON is initialized
    return _hijjet2->nsg;
}

inline int&
HijJet2::njsg    (int i)
{
    init(); // check COMMON is initialized
    if(i < 1 || i > lenI()) return _dummy;
    return _hijjet2->njsg[i-1];
}

inline int&
HijJet2::iasg	(int i, int k)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	k < 1 || k > lenK() ) return _dummy;

    return _hijjet2->iasg[k-1][i-1];
}

inline int&
HijJet2::k1sg	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _dummy;

    return _hijjet2->k1sg[j-1][i-1];
}

inline int&
HijJet2::k2sg	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _dummy;

    return _hijjet2->k2sg[j-1][i-1];
}

inline float&
HijJet2::pxsg	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet2->pxsg[j-1][i-1];
}

inline float&
HijJet2::pysg	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet2->pysg[j-1][i-1];
}

inline float&
HijJet2::pzsg	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet2->pzsg[j-1][i-1];
}

inline float&
HijJet2::pesg	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet2->pesg[j-1][i-1];
}

inline float&
HijJet2::pmsg	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet2->pmsg[j-1][i-1];
}

#endif
