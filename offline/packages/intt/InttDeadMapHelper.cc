#include "InttDeadMapHelper.h"

#include <cmath>    // for ceil, sqrtl, floor, log, sqrt, cos, erfl, sin
#include <cstdio>   // for printf
#include <cstdlib>  // for rand, EXIT_FAILURE, EXIT_SUCCESS, RAND_MAX
#include <iostream>  // for operator<<, basic_ostream, basic_ostream<>::__os...

InttDeadMapHelper::InttDeadMapHelper()
{
	N = 0;
	points = nullptr;

	n_z = 0.0;
	n_t = 0.0;
	p = 0.0;
	q = 0.0;

	mu = 0;
	sg = 0;

	n_g = 0;
	n_l = 0;
	n_u = 0;

	i_lower = 0;
	i_upper = 0;
	D = 0;
}

InttDeadMapHelper::InttDeadMapHelper(InttDeadMapHelper const& _o)
{
	N = _o.N;
	points = new struct InttDeadMap_Point_s[N];
	for(InttDeadMap_Long_t i = 0; i < N; ++i)
	{
		points[i] = _o.points[i];
	}

	n_z = _o.n_z;
	n_t = _o.n_t;

	p = _o.p;
	q = _o.q;

	mu = _o.mu;
	sg = _o.sg;

	n_g = _o.n_g;
	n_l = _o.n_l;
	n_u = _o.n_u;

	i_lower = _o.i_lower;
	i_upper = _o.i_upper;
	D = _o.D;
}

InttDeadMapHelper& InttDeadMapHelper::operator=(InttDeadMapHelper const& _o)
{
	if(this == &_o)return *this;

	N = _o.N;
	points = new struct InttDeadMap_Point_s[N];
	for(InttDeadMap_Long_t i = 0; i < N; ++i)
	{
		points[i] = _o.points[i];
	}

	n_z = _o.n_z;
	n_t = _o.n_t;

	p = _o.p;
	q = _o.q;

	mu = _o.mu;
	sg = _o.sg;

	n_g = _o.n_g;
	n_l = _o.n_l;
	n_u = _o.n_u;

	i_lower = _o.i_lower;
	i_upper = _o.i_upper;
	D = _o.D;

	return *this;
}

InttDeadMapHelper::InttDeadMapHelper(InttDeadMap_Long_t const& _N, InttDeadMap_Double_t const& _n_z, InttDeadMap_Double_t const& _n_t)
{
	N = _N;
	points = new struct InttDeadMap_Point_s[N];

	n_z = _n_z;
	n_t = _n_t;
	p = 1.0 - erfl(n_z / sqrtl(2.0));
	q = sqrtl(p * (1.0 - p));

	mu = 0;
	sg = 0;

	n_g = 0;
	n_l = 0;
	n_u = 0;

	i_lower = 0;
	i_upper = N;
	D = 0;
}

InttDeadMapHelper::~InttDeadMapHelper()
{
	if(points)delete[] points;
}

void InttDeadMapHelper::quicksort_hitkey(InttDeadMap_Long_t const& _l, InttDeadMap_Long_t const& _u)
{
	if(_l < 0)return;
	if(_u < 1)return;
	if(_u <= _l)return;
	if(!points)return;

	InttDeadMap_Long_t i_pivot = _l - 1;

	{
		InttDeadMap_Long_t i;

		struct InttDeadMap_Point_s pivot = points[_u - 1];
		struct InttDeadMap_Point_s temp;

		for(i = _l; i < _u - 1; ++i)
		{
			if(points[i].hitkey <= pivot.hitkey)
			{
				++i_pivot;
				temp = points[i];
				points[i] = points[i_pivot];
				points[i_pivot] = temp;
			}
		}

		++i_pivot;
		temp = points[i_pivot];
		points[i_pivot] = points[_u - 1];
		points[_u - 1] = temp;
	}

	quicksort_hitkey(_l, i_pivot);
	quicksort_hitkey(i_pivot + 1, _u);
}

void InttDeadMapHelper::quicksort_counts(InttDeadMap_Long_t const& _l, InttDeadMap_Long_t const& _u)
{
	if(_l < 0)return;
	if(_u < 1)return;
	if(_u <= _l)return;
	if(!points)return;

	InttDeadMap_Long_t i_pivot = _l - 1;

	{
		InttDeadMap_Long_t i;

		struct InttDeadMap_Point_s pivot = points[_u- 1];
		struct InttDeadMap_Point_s temp;

		for(i = _l; i < _u - 1; ++i)
		{
			if(points[i].counts <= pivot.counts)
			{
				++i_pivot;
				temp = points[i];
				points[i] = points[i_pivot];
				points[i_pivot] = temp;
			}
		}

		++i_pivot;
		temp = points[i_pivot];
		points[i_pivot] = points[_u - 1];
		points[_u - 1] = temp;
	}

	quicksort_counts(_l, i_pivot);
	quicksort_counts(i_pivot + 1, _u);
}


struct InttDeadMapHelper::InttDeadMap_Point_s* InttDeadMapHelper::get_point(InttDeadMap_HitKey_t const& hitkey, InttDeadMap_Long_t const& _l, InttDeadMap_Long_t const& _u)
{
	if(!points)					return nullptr;
	if(points[(_u + _l) / 2].hitkey == hitkey)	return &(points[(i_upper + i_lower) / 2]);

	if(_u - _l <= 1)				return nullptr;

	if(hitkey < points[(_u + _l) / 2].hitkey)	return get_point(hitkey, _l, (_u + _l) / 2);
	if(hitkey > points[(_u + _l) / 2].hitkey)	return get_point(hitkey, (_u + _l) / 2, _u);

	return nullptr;
}

int InttDeadMapHelper::count_type_i()
{
	if(!points)return EXIT_FAILURE;

	//Alternatively, instead of a two-pass and then classify,
	//I can fit to a Gaussian and use a goodness-of-fit parameter
	InttDeadMap_Long_t i = 0;

	InttDeadMap_Double_t lower;
	InttDeadMap_Double_t upper;

	mu = 0;
	sg = 0;

	n_g = 0;
	n_l = 0;
	n_u = 0;

	for(i = i_lower; i < i_upper; ++i)
	{
		mu += points[i].counts;
	}
	mu /= (i_upper - i_lower);

	for(i = i_lower; i < i_upper; ++i)
	{
		sg += (points[i].counts - mu) * (points[i].counts - mu);
	}
	sg /= (i_upper - i_lower);
	sg = sqrtl(sg);

	lower = mu - n_z * sg;
	upper = mu + n_z * sg;

	for(i = i_lower; i < i_upper; ++i)
	{
		if(points[i].counts < lower)
		{
			++n_l;
			points[i].status = InttDeadMap_Status_LOWER;
			continue;
		}

		if(upper < points[i].counts)
		{
			++n_u;
			points[i].status = InttDeadMap_Status_UPPER;
			continue;
		}

		++n_g;
		points[i].status = InttDeadMap_Status_GOOD;
	}

	return EXIT_SUCCESS;
}

int InttDeadMapHelper::check_type_i()
{
	if(!points)return EXIT_FAILURE;

	printf("check_type_i\n");

	if(i_upper - i_lower <= 1)return EXIT_FAILURE;

	count_type_i();

	{
		int flag = 0;

		InttDeadMap_Long_t d;

		InttDeadMap_Double_t lower = p * (i_upper - i_lower) - n_t * q * sqrtl(i_upper - i_lower);
		InttDeadMap_Double_t upper = p * (i_upper - i_lower) + n_t * q * sqrtl(i_upper - i_lower);

		d = 0;
		if(n_l + n_u < floor(lower))
		{
			d = floor(lower) - n_l - n_u;
		}
		if(ceil(upper) < n_l + n_u)
		{
			d = n_l + n_u - ceil(upper);
		}
		if(n_l - n_u > ceil((upper - lower) / 2.0))
		{
			d = n_l - n_u - ceil((upper - lower) / 2.0);
		}
		if(n_u - n_l > ceil((upper - lower) / 2.0))
		{
			d = n_u - n_l - ceil((upper - lower) / 2.0);
		}

		std::cout << "\t" << floor(lower) << " < " << n_l << " + " << n_u << " < " << ceil(upper) << std::endl;
		std::cout << "\t" << "mu: " << mu << "\tsg: " << sg << std::endl;
		std::cout << "\t" << "lower: " << i_lower << "\tupper: " << i_upper << std::endl;
		std::cout << "\td_l: " << (mu - points[i_lower].counts) / sg;
		std::cout << "\td_u: " << (points[i_upper - 1].counts - mu) / sg << std::endl;
		if(d == 0)return EXIT_SUCCESS;

		lower = mu - n_z * sg;
		upper = mu + n_z * sg;

		flag = 0;
		flag = flag || points[i_lower].status != InttDeadMap_Status_GOOD;
		flag = flag || points[i_upper - 1].status != InttDeadMap_Status_GOOD;

		D = 0;
		while(d)
		{
			if(i_upper - i_lower <= 1)break;

			if(points[i_upper - 1].counts - mu < mu - points[i_lower].counts)
			{
				if(flag && points[i_lower].status == InttDeadMap_Status_GOOD)break;
				points[i_lower].status = InttDeadMap_Status_LOWER;
				++i_lower;
			}
			else
			{
				if(flag && points[i_upper - 1].status == InttDeadMap_Status_GOOD)break;
				points[i_upper - 1].status = InttDeadMap_Status_UPPER;
				--i_upper;
			}

			--d;
			++D;
		}
	}

	return check_type_i();
}

void InttDeadMapHelper::gen_points(InttDeadMap_Double_t const& frac)
{
	if(!points)return;

	InttDeadMap_Long_t n_bad = ceil(frac * N / 2.0);

	mu = 4121.1324151;
	sg = 124.14124;

	InttDeadMap_Double_t z[2];
	InttDeadMap_Double_t u[2];

	InttDeadMap_Long_t i;
	InttDeadMap_HitKey_t h;
	InttDeadMap_Double_t c;

	for(i = 0; i < N; ++i)
	{
		if(i % 2 == 0)
		{
			u[0] = (InttDeadMap_Double_t)rand() / (InttDeadMap_Double_t)RAND_MAX;
			u[1] = (InttDeadMap_Double_t)rand() / (InttDeadMap_Double_t)RAND_MAX;

			z[0] = sqrt(-2.0 * log(u[0])) * cos(2.0 * 3.14159265358979 * u[1]);
			z[1] = sqrt(-2.0 * log(u[0])) * sin(2.0 * 3.14159265358979 * u[1]);
		}

		c = mu + sg * z[i % 2];

		points[i].hitkey = i;
		points[i].counts = c;
		points[i].status = InttDeadMap_Status_GOOD;
	}

	c = 999999.999;
	for(i = 0; i < n_bad; ++i)
	{
		h = rand() % N;
		points[h].counts = c;
	}

	c = 0.0;
	for(i = 0; i < n_bad; ++i)
	{
		h = rand() % N;
		points[h].counts = c;
	}
}



