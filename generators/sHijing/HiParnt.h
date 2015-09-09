
#ifndef HiParnt_h
#define HiParnt_h

extern "C" { void* hiparnt_address_(void); }

/**
@class HiParnt
@brief      Class definition for HiParnt, which is used
      to modify the Hijing HIPARNT common.
*/
class HiParnt {
public:
    HiParnt();
    ~HiParnt();
    
    float& 	hipr1	(int n);
    int&	ihpr2	(int n);
    float& 	hint1	(int n);
    int&	ihnt2	(int n);
    
    void 	init	(void);

    // return common array lengths
    inline int lenHipr1() const {return _lenHipr1;}
    inline int lenIhpr2() const {return _lenIhpr2;}
    inline int lenHint1() const {return _lenHint1;}
    inline int lenIhnt2() const {return _lenIhnt2;}

private: 

    // Lengths of the COMMONS
    static const int _lenHipr1 = 100;
    static const int _lenIhpr2 = 50;
    static const int _lenHint1 = 100;
    static const int _lenIhnt2 = 50;

    struct HIPARNT;
    friend struct HIPARNT;

    struct HIPARNT
    {
	float 	hipr1[_lenHipr1];
	int	ihpr2[_lenIhpr2];
	float 	hint1[_lenHint1];
	int	ihnt2[_lenIhnt2];
    };

    int _dummy;
    float _realdummy;
    static HIPARNT* _hiparnt;
};

// set pointer to zero at start
HiParnt::HIPARNT* HiParnt::_hiparnt =0;

inline void
HiParnt::init(void)
{ if (!_hiparnt) _hiparnt = static_cast<HIPARNT*>(hiparnt_address_()); }

inline 
HiParnt::HiParnt() 
    : _dummy		(-999),
      _realdummy	(-999.)
{}

inline 
HiParnt::~HiParnt() 
{}

inline float&
HiParnt::hipr1	(int n)
{
    init(); // check COMMON is initialized
    if(n < 1 || n > lenHipr1()) return _realdummy;
    return _hiparnt->hipr1[n-1];
}

inline int&
HiParnt::ihpr2	(int n)
{
    init(); // check COMMON is initialized
    if(n < 1 || n > lenIhpr2()) return _dummy;
    return _hiparnt->ihpr2[n-1];
}

inline float&
HiParnt::hint1	(int n)
{
    init(); // check COMMON is initialized
    if(n < 1 || n > lenHint1()) return _realdummy;
    return _hiparnt->hint1[n-1];
}

// access ihnt2 in common
inline int&
HiParnt::ihnt2	(int n)
{
    init(); // check COMMON is initialized
    if(n < 1 || n > lenIhnt2()) return _dummy;
    return _hiparnt->ihnt2[n-1];
}

#endif
