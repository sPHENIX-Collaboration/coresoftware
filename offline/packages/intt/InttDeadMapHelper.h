#ifndef INTT_DEAD_MAP_HELPER_H
#define INTT_DEAD_MAP_HELPER_H


class InttDeadMapHelper
{
public:
	typedef int		InttDeadMap_HitKey_t;
	typedef long double	InttDeadMap_Double_t;
	typedef long long	InttDeadMap_Long_t;

	enum InttDeadMap_Status_e
	{
		InttDeadMap_Status_GOOD,
		InttDeadMap_Status_LOWER,
		InttDeadMap_Status_UPPER
	};

	struct InttDeadMap_Point_s
	{
		InttDeadMap_HitKey_t	hitkey;
		InttDeadMap_Long_t	counts;
		InttDeadMap_Status_e	status;
	};

	InttDeadMapHelper();

	InttDeadMapHelper(InttDeadMapHelper const&);
	InttDeadMapHelper& operator=(InttDeadMapHelper const&);

	InttDeadMapHelper(InttDeadMap_Long_t const&, InttDeadMap_Double_t const&, InttDeadMap_Double_t const&);

	~InttDeadMapHelper();

	void gen_points(InttDeadMap_Double_t const&);

//commented while debugging
//protected:
	struct InttDeadMap_Point_s* get_point(InttDeadMap_HitKey_t const&, InttDeadMap_Long_t const&, InttDeadMap_Long_t const&);

	int count_type_i();
	int check_type_i();

	void quicksort_hitkey(InttDeadMap_Long_t const&, InttDeadMap_Long_t const&);
	void quicksort_counts(InttDeadMap_Long_t const&, InttDeadMap_Long_t const&);


	InttDeadMap_Long_t	N;
	InttDeadMap_Point_s*	points;

	InttDeadMap_Double_t	n_z;
	InttDeadMap_Double_t	n_t;
	InttDeadMap_Double_t	p;
	InttDeadMap_Double_t	q;

	InttDeadMap_Double_t	mu;
	InttDeadMap_Double_t	sg;

	InttDeadMap_Long_t	n_g;
	InttDeadMap_Long_t	n_l;
	InttDeadMap_Long_t	n_u;

	InttDeadMap_Long_t	i_lower;
	InttDeadMap_Long_t	i_upper;
	InttDeadMap_Long_t	D;
};

#endif
