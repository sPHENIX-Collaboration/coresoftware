
#ifndef HijJet1_h
#define HijJet1_h

extern "C" { void* hijjet1_address_(void); }
/**
@brief   Class definition for HijJet1, which is used
      to modify the Hijing HIJJET1 common.
*/
class HijJet1 {
public:
    HijJet1();
    ~HijJet1();
    
    int&    	npj	(int i);
    int&    	kfpj	(int i, int j);
    float&  	pjpx	(int i, int j);
    float&  	pjpy	(int i, int j);
    float&  	pjpz	(int i, int j);
    float&  	pjpe	(int i, int j);
    float&  	pjpm	(int i, int j);
    int&    	ntj	(int i);
    int&    	kftj	(int i, int j);
    float&  	pjtx	(int i, int j);
    float&  	pjty	(int i, int j);
    float&  	pjtz	(int i, int j);
    float&  	pjte	(int i, int j);
    float&  	pjtm	(int i, int j);
    
    void	init	(void);

    // return common array lengths
    inline int	lenI() const {return _lenI;}
    inline int	lenJ() const {return _lenJ;}

private: 

    // Lengths of array in HiMain2 common
    static const int _lenI	= 300;
    static const int _lenJ	= 500;

    struct HIJJET1;
    friend struct HIJJET1;

    struct HIJJET1
    {
	int    	npj	[_lenI];
	int    	kfpj	[_lenJ][_lenI];
	float  	pjpx	[_lenJ][_lenI];
	float  	pjpy	[_lenJ][_lenI];
	float  	pjpz	[_lenJ][_lenI];
	float  	pjpe	[_lenJ][_lenI];
	float  	pjpm	[_lenJ][_lenI];
	int    	ntj	[_lenI];
	int    	kftj	[_lenJ][_lenI];
	float  	pjtx	[_lenJ][_lenI];
	float  	pjty	[_lenJ][_lenI];
	float  	pjtz	[_lenJ][_lenI];
	float  	pjte	[_lenJ][_lenI];
	float  	pjtm	[_lenJ][_lenI];
    };

    int  _dummy;
    float  _realdummy;

    static HIJJET1* _hijjet1;
};

// set pointer to zero at start
HijJet1::HIJJET1* HijJet1::_hijjet1 =0;

inline void
HijJet1::init(void)
{ if (!_hijjet1) _hijjet1 = static_cast<HIJJET1*>(hijjet1_address_()); }

// Constructor
inline
HijJet1::HijJet1()
    : _dummy		(-999),
      _realdummy	(-999.)
{}

// Destructor
inline
HijJet1::~HijJet1()
{}

inline int&
HijJet1::npj	(int i)
{
    init(); // check COMMON is initialized
    if(i < 1 || i > lenI()) return _dummy;
    return _hijjet1->npj[i-1];
}

inline int&
HijJet1::kfpj	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _dummy;

    return _hijjet1->kfpj[j-1][i-1];
}

inline float&
HijJet1::pjpx	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet1->pjpx[j-1][i-1];
}

inline float&
HijJet1::pjpy	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet1->pjpy[j-1][i-1];
}

inline float&
HijJet1::pjpz	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet1->pjpz[j-1][i-1];
}

inline float&
HijJet1::pjpe	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet1->pjpe[j-1][i-1];
}

inline float&
HijJet1::pjpm	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet1->pjpm[j-1][i-1];
}

inline int&
HijJet1::ntj	(int i)
{
    init(); // check COMMON is initialized
    if(i < 1 || i > lenI()) return _dummy;
    return _hijjet1->ntj[i-1];
}

inline int&
HijJet1::kftj	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _dummy;

    return _hijjet1->kftj[j-1][i-1];
}

inline float&
HijJet1::pjtx	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet1->pjtx[j-1][i-1];
}

inline float&
HijJet1::pjty	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet1->pjty[j-1][i-1];
}

inline float&
HijJet1::pjtz	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet1->pjtz[j-1][i-1];
}

inline float&
HijJet1::pjte	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet1->pjte[j-1][i-1];
}

inline float&
HijJet1::pjtm	(int i, int j)
{
    init(); // check COMMON is initialized
    if( i < 1 || i > lenI() ||
	j < 1 || j > lenJ() ) return _realdummy;

    return _hijjet1->pjtm[j-1][i-1];
}

#endif
