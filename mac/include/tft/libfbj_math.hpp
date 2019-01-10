// body file: libfbj_math.cpp

#ifndef LIBFBJ_MATH
#define LIBFBJ_MATH

#include <boost/math/distributions.hpp>
#include <set>
#include <map>
#include <vector>
#include <cfloat>			// for DBL_MAX
#include <cmath>			// for std::fabs std::isnan() std::isinf() [always use std::, http://stackoverflow.com/questions/19022561/porting-isnan-to-c11]
#include <algorithm>		// for std::max [template func, safer (no type conversion) and faster (potential inlining) than fmax]
#include "libfbj_base.hpp"	// for each_element each_interval exit_error
#include "libfbj_rng.hpp"	// for rng
#include "libfbj_fet.hpp"	// for fet

// round() in cmath not right. round(1)=0.5. So I write my own
inline double my_round(double number)
{
	return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

// ------- stat distr: qxxx = quantile function = inv of cdf; pdf is density function; cdf is distribution function (<=); _q is >= -------

// students t, nu = d.f = sample size - 1
double cdf_t_q(double t,double nu);
double cdf_t_p(double t,double nu);
double cdf_t_2sided_pv(double t, double nu);
double qstudents_t(double p, double nu);
// chi-square distribution, nu is double, it makes a difference
double qchisq_1df(double p);
double qchisq(double p,double nu);
double cdf_chisq_1df_q(double x);
double cdf_chisq_1df_p(double x);
double cdf_chisq_q(double x,double nu);
double cdf_chisq_p(double x,double nu);
// F distribution
double cdf_f_q(double x, double nu1, double nu2);
double cdf_f_p(double x, double nu1, double nu2);
// standard normal distribution
double pdf_norms(double z);
double cdf_norms(double z);		// same as cdf_norms_p()
double cdf_norms_q(double z);
double cdf_norms_p(double z);
double cdf_norms_2sided_pv(double z);
double qnorms(double p);
// normal distribution
double cdf_norm_q(double mu, double sd, double z);
double cdf_norm_p(double mu, double sd, double z);
double cdf_norm_2sided_pv(double mu, double sd, double z);
double qnorm(double mu, double sd, double p);
// binomial distribution
double pdf_binomial     (double k, double p, double n); // prob[get   k from B(n,p)], k should be an int and >=0
double cdf_binomial_le  (double k, double p, double n); // prob[get 0-k from B(n,p)], k should be an int and >=0
double cdf_binomial_ge  (double k, double p, double n); // prob[get k-n from B(n,p)], k should be an int and >=0
double cdf_binomial_tail(double k, double p, double n); // the smaller tail, 2 * cdf_binomial_tail() = 2-sided p-value

// for multinomial proportion calculation
double prior_strength(int N, int k);

// --------------------- truncated normal distribution ---------------------

double pdf_tnorm(double mu, double sd, double a, double b, double x);
double cdf_tnorm(double mu, double sd, double a, double b, double x);
double pdf_tnorm_UpperTail(double mu, double sd, double a, double x);
double cdf_tnorm_UpperTail(double mu, double sd, double a, double x);
double pdf_tnorm_LowerTail(double mu, double sd, double b, double x);
double cdf_tnorm_LowerTail(double mu, double sd, double b, double x);

// --------------------- Kolmogorov distribution ---------------------

// The D returned by Kolmogorov_Smirnov() is correct, as tested with R's ks.test & STATA's ksmirnov.
// However, the p-value (1-Kolmogorov_K / Kolmogorov_PVal) doesn't match, and all programs give a different P value too.
// They are R STATA GeorgeMarsaglia and an online calculator (http://www.physics.csbsju.edu/stats/KS-test.n.plot_form.html)
// I like 1-Kolmogorov_K because it's very close to STATA's exact P value, and R's ks.test

// Wang et al. (2007) AJHG 81:1278. Correction: pp1279 ES(S)=max{..} => max|..|
template <typename CONTAINER> // vector<double>
double Kolmogorov_Smirnov_like(CONTAINER& stat1, CONTAINER& stat0, bool TwoSided=true)
{
	typedef typename CONTAINER::value_type ElementT;
	std::multimap< ElementT,int,std::greater<ElementT> > combined;
	for (each_element(stat1,it)) combined.insert(std::pair<ElementT,int>(*it,1));
	for (each_element(stat0,it)) combined.insert(std::pair<ElementT,int>(*it,0));
	int N=combined.size();
	double miss = 1/((double)N-stat1.size());
	double Nr=0; for (each_element(stat1,it)) Nr+=*it;
	double max_s=0, score=0;
	for (each_element(combined,it))
	{
		if (it->second)	score += it->first/Nr;
		else			score -= miss;
		double this_score;
		if (TwoSided) this_score=std::fabs(score); else this_score=score;
		if (this_score>max_s) max_s=this_score;
	}
	return max_s;
}

// Wang et al. (2007) AJHG 81:1278. Correction: pp1279 ES(S)=max{..} => max|..|
template <typename CONTAINER> // vector<double>
double Kolmogorov_Smirnov(CONTAINER& stat1, CONTAINER& stat0, bool TwoSided=true)
{
	typedef typename CONTAINER::value_type ElementT;
	std::multimap< ElementT,int,std::greater<ElementT> > combined;
	for (each_element(stat1,it)) combined.insert(std::pair<ElementT,int>(*it,1));
	for (each_element(stat0,it)) combined.insert(std::pair<ElementT,int>(*it,0));
	int N=combined.size();
	double Nr = stat1.size();
	double hits = 1/ Nr;
	double miss = 1/(N-Nr);
	double max_s=0, score=0;
	for (each_element(combined,it))
	{
		if (it->second)	score += hits;
		else			score -= miss;
		double this_score;
		if (TwoSided) this_score=std::fabs(score); else this_score=score;
		if (this_score>max_s) max_s=this_score;
	}
	return max_s;
}

// Wang et al. (2007) AJHG 81:1278. Correction: pp1279 ES(S)=max{..} => max|..|
// combined is a container of statistics from both groups
// idx is a container of indices about whether the element belong to a test group (0/1,false/true)
template <typename CONTAINER_STA, typename CONTAINER_IDX> // multiset<double>, vector<int>
double Kolmogorov_Smirnov_UsingIndex(CONTAINER_STA& stats, CONTAINER_IDX& index, bool TwoSided=true)
{
	int N=stats.size();
	double Nr = 0; for (each_element(index,it)) if (*it) ++Nr;
	double hits = 1/ Nr;
	double miss = 1/(N-Nr);
	double max_s=0, score=0;
	typename CONTAINER_STA::iterator its(stats.begin());
	typename CONTAINER_IDX::iterator iti(index.begin());
	for (; its!=stats.end()&&iti!=index.end(); ++its,++iti)
	{
		if (*iti)	score += hits;
		else		score -= miss;
		double this_score;
		if (TwoSided) this_score=std::fabs(score); else this_score=score;
		if (this_score>max_s) max_s=this_score;
	}
	return max_s;
}

// George Marsaglia, Wai Wan Tsang, Jingbo Wang. Evaluating Kolmogorov's Distribution. J Stat Software (2003) 8:18.
double Kolmogorov_K(int n, double d);
double Kolmogorov_K(double n1, double n2, double d);

// George Marsaglia, Wai Wan Tsang, Jingbo Wang. Evaluating Kolmogorov's Distribution. J Stat Software (2003) 8:18.
double Kolmogorov_PVal(double n, double d);

// For two-sample test, n = n.x*n.y/(n.x + n.y). Ref: R::ks.test
double Kolmogorov_PVal(double n1, double n2, double d);

// George Marsaglia, Wai Wan Tsang, Jingbo Wang. Evaluating Kolmogorov's Distribution. J Stat Software (2003) 8:18.
double Kolmogorov_PVal_NegLog(double n, double d);

// For two-sample test, n = n.x*n.y/(n.x + n.y). Ref: R::ks.test
double Kolmogorov_PVal_NegLog(double n1, double n2, double d);

// ------------------------- basic functions -----------------------------

template< class T >
const T& min_allow_nan(const T& a, const T& b)
{
	if (std::isnan(b)) return a;
	if (std::isnan(a)) return b;
	return std::min(a,b);
}

template< class T >
const T& max_allow_nan(const T& a, const T& b)
{
	if (std::isnan(b)) return a;
	if (std::isnan(a)) return b;
	return std::max(a,b);
}

template< class T >
const T& min_among(const T& a, const T& b, const T& c)
{
	return min_allow_nan(a,min_allow_nan(b,c));
}

template< class T >
const T& max_among(const T& a, const T& b, const T& c)
{
	return max_allow_nan(a,max_allow_nan(b,c));
}

template< class T >
const T& min_among(const T& a, const T& b, const T& c, const T& d)
{
	return min_allow_nan(a,min_allow_nan(b,min_allow_nan(c,d)));
}

template< class T >
const T& max_among(const T& a, const T& b, const T& c, const T& d)
{
	return max_allow_nan(a,max_allow_nan(b,max_allow_nan(c,d)));
}

template< class T >
const T& min_among(const T& a, const T& b, const T& c, const T& d, const T& e)
{
	return min_allow_nan(a,min_allow_nan(b,min_allow_nan(c,min_allow_nan(d,e))));
}

template< class T >
const T& max_among(const T& a, const T& b, const T& c, const T& d, const T& e)
{
	return max_allow_nan(a,max_allow_nan(b,max_allow_nan(c,max_allow_nan(d,e))));
}

// http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the larger value
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<T>::epsilon()
							* std::max(std::abs(x), std::abs(y))
							* ulp;
}

template <typename T>
inline bool IsInBetween(T& a, T& x, T& y)
{
	if (x<=y)	return (x<=a && a<=y);
	else		return (y<=a && a<=x);
}

// the following 4 functions from
// http://www.mathsisfun.com/combinatorics/combinations-permutations.html
// http://www.mathsisfun.com/combinatorics/combinations-permutations-calculator.html
// A permutation is an ordered combination. Unless specified, repetition is not allowed.
// n is the number of items to choose from, r is number of items one can choose.
// n must >=r , but these function does not check for the range of r for faster computation.
// it's very important to use long double inside the functions, but the return could be double.
long double num_combinations(int n, int r); // num_combinations = n! / r! (n-r)! = n(n-1)..(r+1) / (n-r)(n-r-1)..1
long double num_permutations(int n, int r); // num_permutations = n! / (n-r)! = n(n-1)..(n-r+1)
long double num_combinations_rep(int n, int r); // num_combinations with repitition = (n+r-1)! / r! (n-1)!
long double num_permutations_rep(int n, int r); // num_permutations with repitition = n^r
double log10_num_combinations(int n, int r); // num_combinations = n! / r! (n-r)! = n(n-1)..(r+1) / (n-r)(n-r-1)..1
double log10_num_permutations(int n, int r); // num_permutations = n! / (n-r)! = n(n-1)..(n-r+1)
double log10_num_combinations_rep(int n, int r); // num_combinations with repitition = (n+r-1)! / r! (n-1)!
double log10_num_permutations_rep(int n, int r); // num_permutations with repitition = n^r
//  cout << num_combinations(5, 3) << '\n';		// 10
//  cout << num_permutations(5, 3) << '\n';		// 60
//  cout << num_combinations_rep(5, 3) << '\n';	// 35
//  cout << num_permutations_rep(5, 3) << '\n';	// 125
//  cout << num_combinations(5, 2) << '\n';		// 10
//  cout << num_permutations(5, 2) << '\n';		// 20
//  cout << num_combinations_rep(5, 2) << '\n';	// 15
//  cout << num_permutations_rep(5, 2) << '\n';	// 25

// the following 2 taken from
// Understanding the relationship between risks and odds ratios.
// Shrier I, Steele R. Clin J Sport Med. 2006 Mar;16(2):107-10.
// P0 is prevalence of disease in the control group
// (I think the author means unexposed group, because RR = P1 / P0)
double convert_OR_to_RR(double OR, double P0);
double convert_RR_to_OR(double RR, double P0);

// ------------------------------------------- sample size and power -------------------------------------------

// The following 2 functions are based on the article
// A simple method of sample size calculation for linear and logistic regression.
// Hsieh FY, Bloch DA & Larsen MD. Statist. Med. 17, 1623Ð1634 (1998) formula 1 and 2
// none of them were tested yet

double sample_size_logistic_regression_cont(double P1, double beta_star, double alpha=0.05, double power=0.8);

// P1 and P2 are the event rates at X=0 and 1, respectively. B is the proportion of the sample with X=1.
double sample_size_logistic_regression_bin(double P1, double P2, double B, double alpha=0.05, double power=0.8);

// Kelsey et al. Methods in Observational Epidemiology (1986) Table 10-15
// online http://www.sph.emory.edu/~cdckms/sample%20size%202%20grps%20case%20control.html
//     r+1    pb(1-pb)(Zbeta+Zalpha/2)^2
// n = --- x ---------------------------
//      r            (p1-p2)^2
// r is ratio of controls to cases for CC / ratio of unexposed to exposed samlpes for cohort
// p1 is % of cs exposed for CC / % of exposed samples develop the dz for cohort
// p0 is % of ct exposed for CC / % of unexposed samples develop the dz for cohort
// OR = odds ratio
// return n, number of cases for CC / number of exposed individuals studied for cohort
double sample_size_Kelsey_bin(double p0, double OR, double r=1, double alpha=0.05, double power=0.8);

// sd = standard deviation of the studied variable in the population, di = d* (diff in % or means)
double sample_size_Kelsey_cont(double sd, double di, double r=1, double alpha=0.05, double power=0.8);

// the following 2 functions: give n, return power
double power_Kelsey_bin(double p0, double OR, int n, double r=1, double alpha=0.05);
double power_Kelsey_cont(double sd, double di, int n, double r=1, double alpha=0.05);

// taken from the javascript from http://www.stat.ubc.ca/~rollin/stats/ssize/b2.html
// Reference: The cal are the customary ones based on the normal approximation to the binomial distribution.
// Ref: Categorical Data - Estimation of Sample Size and Power for Comparing Two Binomial
// Proportions in Bernard Rosner's Fundamentals of Biostatistics.
int sample_size_proportion_2samples(double p1, double p2, int side=2, double alpha=0.05, double power=0.8);

struct genetic_power_analysis_return_type {
	double power_dom;
	double power_rec;
	double power_add;
};

// power of caes-control SNP association study on binary trait
// Menashe I, Rosenberg PS, Chen BE. PGA: power cal for CC genet assoc analyses. BMC Genet (2008) 13;9:36.
// results not right, need Debug.
genetic_power_analysis_return_type power_CC_SNP_assoc(double pd, double p1, double Dp, double r2,
													double pr, double rr1, double rr2,
													int n1, int n0,
													double alpha=0.05);

// ------------------------------------------- Cochran_Armitage_test_for_trend -------------------------------------------

// http://en.wikipedia.org/wiki/Cochran-Armitage_test_for_trend
template <typename T1,typename T2,typename T3>
void Cochran_Armitage_test_for_trend(T1 N[][2], T2 t[], int k, T3& x2, T3& p);

// http://en.wikipedia.org/wiki/Cochran-Armitage_test_for_trend
template <typename T1,typename T2,typename T3>
void Cochran_Armitage_test_for_trend(T1 N[], T2 t[], int k, T3& x2, T3& p);

// ------------------------------------------- chi-square test -------------------------------------------

// likelihodd ratio test for a contingency table, from sas
template <typename T1,typename T2, typename T3>
void chi_square_test_lr(T1* obs, T2* expct,int num, int df, T3& x2, T3& p);

template <typename T1,typename T2, typename T3>
void chi_square_test_pearson(T1* obs, T2* expct,int num, int df, T3& x2, T3& p);

// http://v8doc.sas.com/sashtml/stat/chap28/sect19.htm , make more sense than wikipedia
template <typename T1,typename T2, typename T3>
void chi_square_test_yate(T1* obs, T2* expct,int num, int df, T3& x2, T3& p);

// ------------------------------------------- contingency table -------------------------------------------
// all statistic function use T* so that both C array and C++ vector can be used for input. This is also the
// reason why I don't put the data inside the class, because otherwise I need to specify how to store the data.
// However, because of this user needs to make sure prepare_from() is called before using any stat functions!
// I don't put ContingencyTableStatType inside struct or class for shorter names, but it pollutes the namespace.
// All functions have tested for VAL & ASE, but pvl is note tested because SAS doesn't show it. Be careful.

enum ContingencyTableStatType
{
	ChiSq_Ps, ChiSq_Yt, ChiSq_LR, UncertCR, UncertRC, UncertSm, LambdaCR, LambdaRC, LambdaSm,	// nominal
	KendllTb, StuartTc, GK_gamma, SomerDcr, SomerDrc, PearsonR, SpearmnR, JTtrendS,	MHChiSqT,	// ordinal
	BowkersQ, KappaSmp,																			// RxR
	CochranQ,																					// other
	NumStats																					// # of stats
};

template <typename T>
class ContingencyTable
{
public:
	std::vector<double> val, ase, pvl, dof; // value, Asymptotic Standard Error, p-value, degree of freedom
	void InitAllStats() { val.assign(NumStats,-2); ase.assign(NumStats,0); pvl.assign(NumStats,-1); dof.assign(NumStats,0); }
	void InitStat(ContingencyTableStatType type) { val[type]=-2; ase[type]=0; pvl[type]=-1; dof[type]=0; }
	
	unsigned nr, nc, nd;			// #row, #col, #data
	double n;						// sum of all cells
	std::vector<double> sumr, sumc; // sum of row, col
	ContingencyTable():nr(0),nc(0),nd(0),n(0) {}
	void prepare_from(T* data, unsigned r, unsigned c)
	{
		nr=r; nc=c; nd=r*c;
		sumr.assign(nr,0);
		sumc.assign(nc,0);
		n=0;
		for (unsigned ij=0,i=0; i<nr; ++i)
			for (unsigned  j=0; j<nc; ++j, ++ij)
			{
				n      +=data[ij];
				sumr[i]+=data[ij];
				sumc[j]+=data[ij];
			}
		InitAllStats();
	}
};

// has segmentation fault: 11
template <typename T>
void RxC_both_nominal(T* N, ContingencyTable<T>& ctab);

// http://v8doc.sas.com/sashtml/stat/chap28/sect18.htm
template <typename T>
void RxC_both_ordinal(T* N, ContingencyTable<T>& ctab);

template <typename T>
void RxR_agreement(T* N, ContingencyTable<T>& ctab);

// http://en.wikipedia.org/wiki/Cochran's_Q_test
// tested with http://support.sas.com/documentation/cdl/en/procstat/63104/HTML/default/viewer.htm#procstat_freq_sect034.htm
template <typename T>
void RxC_BinDV_CochranQ(T* N, ContingencyTable<T>& ctab);

// Chi-squared and Fisher–Irwin tests of two-by-two tables with small sample recommendations. Ian Campbell. Statist. Med. 2007; 26:3661–3675
template <typename T1,typename T2>
inline void test_2x2(T1 table[2][2],T2& p)
{
	double expct[2][2];
	double total_0x = table[0][0]+table[0][1]; if (!total_0x) {p=1; return;}
	double total_1x = table[1][0]+table[1][1]; if (!total_1x) {p=1; return;}
	double total_x0 = table[0][0]+table[1][0]; if (!total_x0) {p=1; return;}
	double total_x1 = table[0][1]+table[1][1]; if (!total_x1) {p=1; return;}
	double total_xx = total_0x+total_1x;
	expct[0][0]=total_0x*total_x0/total_xx;
	expct[0][1]=total_0x*total_x1/total_xx;
	expct[1][0]=total_1x*total_x0/total_xx;
	expct[1][1]=total_1x*total_x1/total_xx;
	
	if (expct[0][0]<1 || expct[0][1]<1 || expct[1][0]<1 || expct[1][1]<1)
	{
		p=Fishers_exact_test_2x2(table,false);
	}
	else
	{	// x2= (ad-bc)^2(a+b+c+d-1) / (a+b)(c+d)(b+d)(a+c)
		double denominator = total_0x * total_1x * total_x0 * total_x1; // cannot be 0
		double x2 = table[0][0] * table[1][1]  -  table[0][1] * table[1][0];
		x2 = x2 * x2 * (total_xx-1) / denominator;
		p=cdf_chisq_1df_q(x2);
	}
}

// Chi-squared and Fisher–Irwin tests of two-by-two tables with small sample recommendations. Ian Campbell. Statist. Med. 2007; 26:3661–3675
template <typename T1,typename T2>
inline void test_2x2_array(T1 table[4],T2& p)
{
	double expct[4];
	double total_0x = table[0]+table[1]; if (!total_0x) {p=1; return;}
	double total_1x = table[2]+table[3]; if (!total_1x) {p=1; return;}
	double total_x0 = table[0]+table[2]; if (!total_x0) {p=1; return;}
	double total_x1 = table[1]+table[3]; if (!total_x1) {p=1; return;}
	double total_xx = total_0x+total_1x;
	expct[0]=total_0x*total_x0/total_xx;
	expct[1]=total_0x*total_x1/total_xx;
	expct[2]=total_1x*total_x0/total_xx;
	expct[3]=total_1x*total_x1/total_xx;
	
	if (expct[0]<1 || expct[1]<1 || expct[2]<1 || expct[3]<1)
	{
		p=Fishers_exact_test_2x2(table,false);
	}
	else
	{	// x2= (ad-bc)^2(a+b+c+d-1) / (a+b)(c+d)(b+d)(a+c)
		double denominator = total_0x * total_1x * total_x0 * total_x1; // cannot be 0
		double x2 = table[0] * table[3]  -  table[1] * table[2];
		x2 = x2 * x2 * (total_xx-1) / denominator;
		p=cdf_chisq_1df_q(x2);
	}
}

template <typename T1,typename T2>
inline void test_2x2_array_yate(T1 table[4],T2& p)
{
	double expct[4];
	double total_0x = table[0]+table[1]; if (!total_0x) {p=1; return;}
	double total_1x = table[2]+table[3]; if (!total_1x) {p=1; return;}
	double total_x0 = table[0]+table[2]; if (!total_x0) {p=1; return;}
	double total_x1 = table[1]+table[3]; if (!total_x1) {p=1; return;}
	double total_xx = total_0x+total_1x;
	expct[0]=total_0x*total_x0/total_xx;
	expct[1]=total_0x*total_x1/total_xx;
	expct[2]=total_1x*total_x0/total_xx;
	expct[3]=total_1x*total_x1/total_xx;
	double x2;
	chi_square_test_yate(table, expct, 4, 1, x2, p);
}

// ------------------------------------------- Wilcoxon_signed_rank_test -------------------------------------------

template <typename T1,typename T2>
void Wilcoxon_signed_rank_test(const std::multiset<T1>& z, T2& sumW, T2& p);

// ------------------------------------------- general functions for containers w/o sorting -------------------------------------------
// These functions are robust to NaN, empty, all NaN situations.

template <typename container>
inline int N(const container& s)
{
	int result=0;
	for (each_element(s,it))
		if (!std::isnan(*it)) ++result;
	return result;
}

// problem: the sum could go out of bound
template <typename CONTAINER_T>
inline double sum(const CONTAINER_T& s)
{
	double result=0;
	for (each_element(s,it))
		if (!std::isnan(*it)) result += *it;
	return result;
}

// problem: the sum could go out of bound
template <typename container, typename N_t>
inline void sum_and_N(const container& s, double& a, N_t& n)
{
	a=0; n=0;
	for (each_element(s,it)) if (!std::isnan(*it)) { a+=*it; ++n; }
}

// problem: the sum could go out of bound
template <typename container, typename N_t>
inline void fast_mean_and_N(const container& s, double& m, N_t& n)
{
	m=0; n=0;
	for (each_element(s,it)) if (!std::isnan(*it)) { m+=*it; ++n; }
	if (!n) m  = std::numeric_limits<double>::signaling_NaN();
	else	m /= n;
}

template <typename container, typename N_t>
inline void safe_mean_and_N(const container& s, double& m, N_t& n)
{
	double offset = 0;	for (each_element(s,it)) if (!std::isnan(*it)) { offset =*it; break; }
	m=0; n=0;			for (each_element(s,it)) if (!std::isnan(*it)) { m+=*it-offset; ++n; }
	if (!n) { m=std::numeric_limits<double>::signaling_NaN(); return; }
	m /= n;
	m += offset;
}

// the caveat of the above function is always return double. This function can return other type, but cannot check NaN
template <typename container, typename m_t, typename N_t>
inline void safe_mean_noNaN(const container& s, m_t& m, N_t& n)
{
	m_t offset = 0;	for (each_element(s,it)) { offset =*it; break; }
	m=0; n=0;		for (each_element(s,it)) { m+=*it-offset; ++n; }
	if (!n) return;
	m /= n;
	m += offset;
}

template <typename container>
inline double tss(const container& s, double m) //  total sum of squares (TSS) about m
{
	if (std::isnan(m)) return std::numeric_limits<double>::signaling_NaN();
	double result=0;
	int size=0;
	for (each_element(s,it)) if (!std::isnan(*it)) { result += (*it-m)*(*it-m); ++size; }
	if (size) return result;
	else return std::numeric_limits<double>::signaling_NaN();
}

template <typename container>
inline double tss(const container& s) //  total sum of squares (TSS) about mean
{
	int size=0;
	double m=0;
	safe_mean_and_N(s,m,size);
	if (!size) return std::numeric_limits<double>::signaling_NaN();
	double result=0;
	for (each_element(s,it)) if (!std::isnan(*it)) { result += (*it-m)*(*it-m); }
	return result;
}

template <typename container>
inline double tsn(const container& s, double m, double n) //  total sum of ^n about m
{
	if (std::isnan(m)) return std::numeric_limits<double>::signaling_NaN();
	double result=0;
	int size=0;
	for (each_element(s,it)) if (!std::isnan(*it)) { result += pow(*it-m,n); ++size; }
	if (size) return result;
	else return std::numeric_limits<double>::signaling_NaN();
}

template <typename container>
inline double tsn(const container& s, double n) //  total sum of ^n about mean
{
	int size=0;
	double m=0;
	safe_mean_and_N(s,m,size);
	if (!size) return std::numeric_limits<double>::signaling_NaN();
	double result=0;
	for (each_element(s,it)) if (!std::isnan(*it)) { result += pow(*it-m,n); }
	return result;
}

template <typename container>
inline double population_variance(const container& s) // sigma squared
{
	int size=0;
	double m=0;
	safe_mean_and_N(s,m,size);
	if (size==0) return std::numeric_limits<double>::signaling_NaN();
	return tss(s,m)/size;
}

template <typename container>
inline double sample_variance(const container& s) // s squared
// estimate of population variance based on sampled data, should be the default for var
{
	int size=0;
	double m=0;
	safe_mean_and_N(s,m,size);
	if (size<=1) return std::numeric_limits<double>::signaling_NaN();
	return tss(s,m)/(size-1);
}

template <typename container>
inline void standardize(container& c)
{
	int size = 0;
	double m = 0;
	safe_mean_and_N(c,m,size);
	if (size==0) return;
	double s = sqrt(tss(c,m)/size);
	for (each_element(c,it)) if (!std::isnan(*it)) *it = (*it-m)/s;
}

// ------------------------------------------- general functions for containers w/o sorting, require [], good for Eigen -------------------------------------------

template <typename container>
inline int array_N(const container& s)
{
	int result=0;
	for (int i=s.size()-1; i>=0; --i)
		if (!std::isnan(s[i])) ++result;
	return result;
}

// problem: the sum could go out of bound
template <typename CONTAINER_T>
inline double array_sum(const CONTAINER_T& s)
{
	double result=0;
	for (int i=s.size()-1; i>=0; --i)
		if (!std::isnan(s[i])) result += s[i];
	return result;
}

// problem: the sum could go out of bound
template <typename container, typename N_t>
inline void array_sum_and_N(const container& s, double& a, N_t& n)
{
	a=0; n=0;
	for (int i=s.size()-1; i>=0; --i)
		if (!std::isnan(s[i])) { a+=s[i]; ++n; }
}

// problem: the sum could go out of bound
template <typename container, typename N_t>
inline void array_fast_mean_and_N(const container& s, double& m, N_t& n)
{
	m=0; n=0;
	for (int i=s.size()-1; i>=0; --i)
		if (!std::isnan(s[i])) { m+=s[i]; ++n; }
	if (!n) m  = std::numeric_limits<double>::signaling_NaN();
	else	m /= n;
}

template <typename container, typename N_t>
inline void array_safe_mean_and_N(const container& s, double& m, N_t& n)
{
	double offset = 0;	for (int i=s.size()-1; i>=0; --i) if (!std::isnan(s[i])) { offset =s[i]; break; }
	m=0; n=0;			for (int i=s.size()-1; i>=0; --i) if (!std::isnan(s[i])) { m+=s[i]-offset; ++n; }
	if (!n) { m=std::numeric_limits<double>::signaling_NaN(); return; }
	m /= n;
	m += offset;
}

template <typename container>
inline double array_tss(const container& s, double m) //  total sum of squares (TSS) of data about m (mean)
{
	if (std::isnan(m)) return std::numeric_limits<double>::signaling_NaN();
	double result=0;
	int size=0;
	for (int i=s.size()-1; i>=0; --i) if (!std::isnan(s[i])) { result += (s[i]-m)*(s[i]-m); ++size; }
	if (!size) return std::numeric_limits<double>::signaling_NaN();
	return result;
}

template <typename container>
inline double array_tss(const container& s) //  total sum of squares (TSS) of data about m (mean)
{
	int size=0;
	double m=0;
	array_safe_mean_and_N(s,m,size);
	if (!size) return std::numeric_limits<double>::signaling_NaN();
	double result=0;
	for (int i=s.size()-1; i>=0; --i) if (!std::isnan(s[i])) { result += (s[i]-m)*(s[i]-m); }
	return result;
}

template <typename container>
inline double array_tsn(const container& s, double m, double n) //  total sum of squares (TSS) of data about m (mean)
{
	if (std::isnan(m)) return std::numeric_limits<double>::signaling_NaN();
	double result=0;
	int size=0;
	for (int i=s.size()-1; i>=0; --i) if (!std::isnan(s[i])) { result += pow(s[i]-m,n); ++size; }
	if (!size) return std::numeric_limits<double>::signaling_NaN();
	return result;
}

template <typename container>
inline double array_tsn(const container& s, double n) //  total sum of squares (TSS) of data about m (mean)
{
	int size=0;
	double m=0;
	array_safe_mean_and_N(s,m,size);
	if (!size) return std::numeric_limits<double>::signaling_NaN();
	double result=0;
	for (int i=s.size()-1; i>=0; --i) if (!std::isnan(s[i])) { result += pow(s[i]-m,n); }
	return result;
}

template <typename container>
inline double array_population_variance(const container& s) // sigma squared
{
	int size=0;
	double m=0;
	array_safe_mean_and_N(s,m,size);
	if (size==0) return std::numeric_limits<double>::signaling_NaN();
	return array_tss(s,m)/size;
}

template <typename container>
inline double array_sample_variance(const container& s) // s squared
// estimate of population variance based on sampled data, should be the default for var
{
	int size=0;
	double m=0;
	array_safe_mean_and_N(s,m,size);
	if (size<=1) return std::numeric_limits<double>::signaling_NaN();
	return array_tss(s,m)/(size-1);
}

template <typename container>
inline void array_standardize(container& c)
{
	int size = 0;
	double m = 0;
	array_safe_mean_and_N(c,m,size);
	if (size==0) return;
	double s = sqrt(array_tss(c,m)/size);
	for (int i=c.size()-1; i>=0; --i) if (!std::isnan(c[i])) { c[i] = (c[i]-m)/s; }
}

// ------------------------------------------- general functions for containers -------------------------------------------

// https://class.stanford.edu/courses/Engineering/CS144/Introduction_to_Computer_Networking/wiki/HRP258/example-r-classwork-solutions-using-r/calculating-inner-quartile-range-r/
// http://stat.ethz.ch/R-manual/R-patched/library/stats/html/quantile.html
// http://en.wikipedia.org/wiki/Quartile
// The quantile method is equal to SAS method 5, the default method of SAS.
// quantile median MAD all require that the contaier is already sorted.
// The first 3 functions is fast, but they can handle only container that has [].
// The next 3 functions is slower, but can handle set/multiset, hence the name.
// BEWARE: quantile median MAD do not check size(). If s.empty(), it will cause segmentaion fault!

template <typename CONTAINER_T> // SORTED vector/deque/array
inline double quantile(const CONTAINER_T& s, double p) // old name: quantile_SortedVec
{
	assert(p>0);
	assert(p<1);
	
	// get data size. If there're NaN, they must be at the end for sorted vector or set
	int N=s.size();
	for (each_element_rev(s,it)) { if (std::isnan(*it)) --N; else break; }
	if (!N) return std::numeric_limits<double>::signaling_NaN();

	// get the desire location j. If g==0, then mid-way between j and j+1, otherwise j+1
	int j=int(N*p);
	double g=N*p-j;
	
	if (g < std::numeric_limits<double>::epsilon() * N * p * 2)
		return (s[j-1]+s[j])/2.0;
	else
		return s[j];
}

template <typename CONTAINER_T>
inline double median(const CONTAINER_T& s)
{
	return quantile(s,0.5);
}

// previously MAD()=MedianMAD(), but it too frequenctly returns 0, so I changed to MAD()=MeanMAD() below.
template <typename container>				// SORTED vector/deque/array
inline double MedianMAD(const container& s)	// median absolute deviation (previously MAD)
{
	double med = median(s);					// median
	std::deque<double> abs_dev;				// absolute deviation
	for (each_element(s,it)) if (!std::isnan(*it)) abs_dev.push_back(std::fabs(*it-med));
	std::sort(abs_dev.begin(),abs_dev.end());
	return median(abs_dev);
}

template <typename container>			// SORTED vector/deque/array
inline double MAD(const container& s)	// MEAN absolute deviation
{
	double m;
	int N;
	safe_mean_and_N(s,m,N);				// mean
	std::deque<double> abs_dev;			// absolute deviation
	for (each_element(s,it)) if (!std::isnan(*it)) abs_dev.push_back(std::fabs(*it-m));
	safe_mean_and_N(abs_dev,m,N);
	return m;
}

// Container can be set, multiset, SORTED vector/deque/array. But vec/deq/array should use quantile(), hence the name.
// However, multiset is not robust to nan. The location of nan depends on the order of insertion, not always the end or the beginning!
// Therefore, these functions can only be used when there's no nan.
template <typename CONTAINER_T>
inline double quantile_Multiset(const CONTAINER_T& s, double p) // old_name = quantile_SAS_method5
{
	assert(p>0);
	assert(p<1);
	
	// get data size. If there're NaN, they must be at the end for sorted vector or set
	int N=s.size();
	for (each_element_rev(s,it)) { if (std::isnan(*it)) --N; else break; }
	if (!N) return std::numeric_limits<double>::signaling_NaN();
	
	// get the desire location j. If g==0, then mid-way between j and j+1, otherwise j+1
	int j=int(N*p);
	double g=N*p-j;
	double result=std::numeric_limits<double>::signaling_NaN();
	// below see http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon :
	// "the machine epsilon has to be scaled to the magnitude of the larger value
    // and multiplied by the desired precision in ULPs (units in the last place)"
	if ( g < std::numeric_limits<double>::epsilon() * N * p * 2 )
	{
		typename CONTAINER_T::const_iterator it;
		int i;
		for (i=1, it=s.begin(); i<j; ++i, ++it) ;
		result=*it;
		++it;
		result+=*it;
		result/=2.0;
	}
	else
	{
		typename CONTAINER_T::const_iterator it;
		int i;
		for (i=1, it=s.begin(); i<=j; ++i, ++it) ;
		result=*it;
	}
	return result;
}

template <typename CONTAINER_T>
inline double median_Multiset(const CONTAINER_T& s)
{
	return quantile_Multiset(s,0.5);
}

template <typename container>			// SORTED vector/deque/array
inline double MAD_Multiset(const container& s)	// median absolute deviation
{
	double med = median_Multiset(s);	// median
	std::deque<double> abs_dev;			// absolute deviation
	for (each_element(s,it)) if (!std::isnan(*it)) abs_dev.push_back(std::fabs(*it-med));
	std::sort(abs_dev.begin(),abs_dev.end());
	return median(abs_dev);
}

// ------------------------------------------- Values class -------------------------------------------

// Features: 1) lazy computation of statistics only when get() is called
//           2) never do the same statistics again
//           3) automatically do the required hidden step, such as sort() before MAD(), which is easy to forget
//           4) robust to missing values as NaN
//           5) the calculation of mean is safe (not easily go out of bound if there're too many variables)
//           6) fast computation even with large data (use deque instead of vector or multiset)
//           7) modifying the data will clear the internal statistics results automatically
//           8) robust to empty or all NaN
// BEWARE:   1) remove_xxx() do_abp() will remove NaN, but sort() only puts NaN at the end.
// BELOW:    1) The 1st line do not sort data, the 2nd line will sort. The 3rd is for adjusted box plot, will sort too.
//           2) M3=third central moment. SKEW is Sample skewness = N/[(N-1)(N-2)] SUM[((Xi-Mean)/s)^3].

struct STAT {
	enum StatType { MIN,MAX,MEAN,SUM,N,TSS,TSS0,RMS,SPL_VAR,POP_VAR,SPL_SD,POP_SD,PRD,M3,SKEW,ADJ_SKEW,SEM,POP_VMR,CV,RSD,UCV,URSD,NUM_NAN,MEAN_ABS,MED_ABS,MEAN_POS,MED_POS,
		IQR,QCD,MAD,BWMV,MEDIAN,PCT001,PCT005,PCT01,PCT05,PCT1,PCT2,PCT2p5,PCT5,PCT9,PCT10,PCT25,PCT75,PCT90,PCT91,PCT95,PCT97p5,PCT98,PCT99,PCT995,PCT999,PCT9995,PCT9999,SD_LO,SD_UP,SORTED,
		ABP_MC,ABP_LF,ABP_UF,ABP_LO,ABP_UO,Adil_LF,Adil_UF,Kimber_LF,Kimber_UF,
		NUM_STAT
	};
	static const std::vector<std::string> StatName;
};

STAT::StatType	to_StatType(const std::string&   name);
std::string		to_StatName(const STAT::StatType type);

template <typename T>
class Values {
public:
	typedef T												value_type;
	typedef typename std::deque<T>::iterator				iterator;
	typedef typename std::deque<T>::reverse_iterator		reverse_iterator;
	typedef typename std::deque<T>::const_iterator			const_iterator;
	typedef typename std::deque<T>::const_reverse_iterator	const_reverse_iterator;
	Values()						{}
	Values(size_t n, const T& val)	{ data.assign(n,val); }
	
	// get or clear result
	int		num_ge	(const T& threshold);
	int		num_le	(const T& threshold);
	double	get		(STAT::StatType type);
	double	pct		(const double& p) { sort(); return quantile(data, p/100); } // percentile p=(0,100)
	double	qt		(const double& q) { sort(); return quantile(data, q); }		// quantile   q=(0,1)
	void	clear_result()	{ num_ge_cnt.clear(); num_le_cnt.clear(); ranks.clear(); result.clear(); }
	
	// Change of data that clear the internal results. Be careful whether SORTED should be maintained.
	template <typename It>
	void	push_back(It start, It end)		{ clear_result(); for (It i=start; i!=end; ++i) data.push_back(*i); }
	void	clear()							{ clear_result(); data.clear();		 }
	void	push_back(const T& t)			{ clear_result(); data.push_back(t); }
	void	push_front(const T& t)			{ clear_result(); data.push_front(t);}
	void	pop_back()						{ clear_result(); data.pop_back();   }
	void	pop_front()						{ clear_result(); data.pop_front();  }
	void	assign(size_t n, const T& val)	{ clear_result(); data.assign(n,val);}
	void	remove_smallest(int n)			{ if (n>0) { remove_nan(); while (n--) data.pop_front(); clear_result(); result[STAT::SORTED]=1; } }
	void	remove_largest(int n)			{ if (n>0) { remove_nan(); while (n--) data.pop_back();  clear_result(); result[STAT::SORTED]=1; } }
	void	remove_extremes(int l, int r)	{ remove_smallest(l); remove_largest(r); }
	void	remove_outliers()				{ remove_extremes(get(STAT::ABP_LO), get(STAT::ABP_UO)); }
	void	standardize() {
		T m = get(STAT::MEAN);
		T s = get(STAT::POP_SD);
		for (auto &x:data) if (!std::isnan(x)) x=(x-m)/s;
		double sorted_bkup = get(STAT::SORTED);
		clear_result();
		if (sorted_bkup==1) result[STAT::SORTED]=sorted_bkup;
	}
//	typename std::deque<T>::iterator erase(const_iterator position)						{ clear_result(); return data.erase(position);	 }
//	typename std::deque<T>::iterator erase(const_iterator first, const_iterator last)	{ clear_result(); return data.erase(first,last); }
	
	// change of data that won't change the internal results
	void	remove_nan();
	void	sort();
	void	rank();
	
	// access data
	typename std::deque<T>::const_iterator			begin()   const	{ return data.cbegin(); }
	typename std::deque<T>::const_iterator			end()	  const	{ return data.cend();   }
	typename std::deque<T>::const_iterator			cbegin()  const	{ return data.cbegin(); }
	typename std::deque<T>::const_iterator			cend()    const	{ return data.cend();   }
	typename std::deque<T>::const_reverse_iterator	rbegin()  const { return data.crbegin();}
	typename std::deque<T>::const_reverse_iterator	rend()    const { return data.crend();  }
	typename std::deque<T>::const_reverse_iterator	crbegin() const { return data.crbegin();}
	typename std::deque<T>::const_reverse_iterator	crend()   const { return data.crend();  }
	bool	 empty()				const { return data.empty(); }
	size_t	 size()					const { return data.size();  }
	
	// access data, may change results, but doesn't clear_result() for fast access (not safe)
	const T& operator[] (size_t i) const	{ return data[i];		}
	const T& front() const					{ return data.front();	}
	const T& back() const					{ return data.back();	}
	T& operator[] (size_t i)				{ return data[i];		}
	T& front()								{ return data.front();	}
	T& back()								{ return data.back();	}
	const std::deque<T>& get_data()			{ return data;			}
	
private:
	// data
	std::deque<T>						data;
	// internal results
	std::deque<double>					ranks;
	std::map<STAT::StatType,double>		result;
	std::map< T,int,std::greater<T> >	num_ge_cnt;
	std::map< T,int >					num_le_cnt;
	// functions
	void do_min();
	void do_max();
	void do_rms();
	void do_pop_vmr();
	void do_cv();
	void do_mean_abs();
	void do_med_abs();
	void do_mean_pos();
	void do_med_pos();
	void do_n();
	void do_spl_var();
	void do_pop_var();
	void do_prd();
	void do_m3();
	void do_adj_skew();
	void do_skew();
	void do_sem();
	void do_abp();
	void do_kimber();
	void do_bwmv();
	void do_sd1sided();
};

// --------------------- t test / Kruskal-Wallis ---------------------

template <typename T>
double two_sample_t_same_sd(Values<T>& sample1, Values<T>& sample2)
{
	double n1 = sample1.get(STAT::N);
	double n2 = sample2.get(STAT::N);
	double df = n1+n2-2;
	if (df<1) return std::numeric_limits<double>::signaling_NaN();
	double s1 = sample1.get(STAT::SPL_SD);
	double s2 = sample2.get(STAT::SPL_SD);
	double m1 = sample1.get(STAT::MEAN);
	double m2 = sample2.get(STAT::MEAN);
	double Sp = sqrt(((n1-1)*s1*s1 + (n2-1)*s2*s2)/(n1+n2-2));
	double t = (m1-m2) / (Sp*sqrt(1/n1+1/n2));
	double pv = cdf_t_2sided_pv(t,df);
	return pv;
}

// allow any number of groups. If just two groups, it is equivalent to Mann-Whitney U-test.
// assumed that all groups have a distribution with the same shape
// alternative hypothesis: at least one population median of one group is different from the population median of at least one other group
// If each group has <5 samples, the KW test will has little power and the p-value will be always > 0.05.
// Some say "small" means "total sample size is seven or less". -- http://www.graphpad.com/guides/prism/6/statistics/index.htm?how_the_kruskal-wallis_test_works.htm
template <typename T1, typename T2>
double Kruskal_Wallis(std::deque< std::pair<T1,T2> >& data, // input, nonNaN value->group
					  std::vector<int>& ni,					// output, number of observations in group i
					  std::vector<double>& ri,				// output, average rank of observations in group i
					  std::vector<T2>& gi)					// output, ID of group i
{
	std::sort(data.begin(),data.end());
	std::deque<double> ranks;
	ranks.assign(data.size(),0);
	int tie_bgn=0, tie_end=0;
	for (size_t i=0; i<data.size(); ++i)
	{
		ranks[i]=i+1;
		if (i)
		{
			if (data[i].first==data[i-1].first)
			{
				++tie_end;
				if (i+1==data.size())
				{
					double mean = (tie_bgn+tie_end+2)/2.0;
					for (int j=tie_bgn; j<=tie_end; ++j) ranks[j]=mean;
				}
			}
			else
			{
				if (tie_bgn!=tie_end)
				{
					double mean = (tie_bgn+tie_end+2)/2.0;
					for (int j=tie_bgn; j<=tie_end; ++j) ranks[j]=mean;
				}
				tie_bgn=i;
				tie_end=i;
			}
		}
	}
	int N=data.size();
	double rb = (N+1)/2.0;
	ni.clear();
	ri.clear();
	gi.clear();
	std::vector<double> rij_rb;
	std::map<T2,int> group2coord;
	for (size_t i=0; i<data.size(); ++i)
	{
		T2& group = data[i].second;
		int coord = -1;
		typename std::map<T2,int>::iterator it = group2coord.find(group);
		if (it==group2coord.end())	{ coord=group2coord.size(); group2coord[group]=coord; ni.push_back(0); ri.push_back(0); gi.push_back(group); rij_rb.push_back(0); }
		else						{ coord=it->second; }
		ri[coord]+=ranks[i];
		ni[coord]+=1;
		rij_rb[coord] += ((ranks[i]-rb)*(ranks[i]-rb));
	}
	int df=ni.size()-1;
	if (df>0)
	{
		double nom=0;
		double den=0;
		// bool small_sample = false;
		for (size_t g=0; g<ni.size(); ++g)
		{
			// if (ni[g]<5) small_sample=true;
			ri[g]/=ni[g];
			nom += (ri[g]-rb)*(ri[g]-rb)*ni[g];
			den += rij_rb[g];
		}
		if (den && nom)
		{
			double H = (N-1)*nom/den;
			return cdf_chisq_q(H, df);
		}
		else
			return 1;
	}
	else
		return 1;
}

template <typename T1, typename T2>
double Kruskal_Wallis(std::deque< std::pair<T1,T2> >& data) // input, nonNaN value->group
{
	std::vector<int> ni;				// output, number of observations in group i
	std::vector<double> ri;				// output, average rank of observations in group i
	std::vector<T2> gi;					// output, ID of group i
	return Kruskal_Wallis(data,ni,ri,gi);
}

// return p-value. If sample size is small (n<10) or no observation in one group (n1==0||n2==0), return NaN.
template <typename T1, typename T2>
double ranksum(std::deque< std::pair<T1,T2> >& data, char test) // input, nonNaN value->group(1/2)
{
	std::sort(data.begin(),data.end());
	std::deque<double> ranks;
	ranks.assign(data.size(),0);
	int tie_bgn=0, tie_end=0;
	double sum_k=0;
	for (size_t i=0; i<data.size(); ++i)
	{
		ranks[i]=i+1;
		if (i)
		{
			if (data[i].first==data[i-1].first)
			{
				++tie_end;
				if (i+1==data.size())
				{
					double mean = (tie_bgn+tie_end+2)/2.0;
					for (int j=tie_bgn; j<=tie_end; ++j) ranks[j]=mean;
					double ti=tie_end-tie_bgn+1;
					sum_k += (ti*ti*ti-ti);
				}
			}
			else
			{
				if (tie_bgn!=tie_end)
				{
					double mean = (tie_bgn+tie_end+2)/2.0;
					for (int j=tie_bgn; j<=tie_end; ++j) ranks[j]=mean;
					double ti=tie_end-tie_bgn+1;
					sum_k += (ti*ti*ti-ti);
				}
				tie_bgn=i;
				tie_end=i;
			}
		}
	}
	double n =data.size();
	if (n<10) return std::numeric_limits<double>::signaling_NaN();
	double n1=0, n2=0;
	double R1=0, R2=0;
	for (size_t i=0; i<data.size(); ++i)
	{
		if 		(data[i].second==1)	{ ++n1; R1+=ranks[i]; }
		else if (data[i].second==2)	{ ++n2; R2+=ranks[i]; }
		else exit_error("Group ID for ranksum() must be 1 or 2.");
	}
	if (n1==0||n2==0) return std::numeric_limits<double>::signaling_NaN();
	double U1=R1-n1*(n1+1)/2;
	//double U2=R2-n2*(n2+1)/2;
	double mu=n1*n2/2;
	double sd=sqrt(n1*n2/12*((n+1)-sum_k/(n*(n-1))));
	double z=(U1-mu)/sd;
	double pv=0.5;
	if 		(test=='<') pv=cdf_norms_p(z);
	else if	(test=='>') pv=cdf_norms_q(z);
	else if	(test=='=') pv=cdf_norms_2sided_pv(z);
	else exit_error("wrong test for ranksum()");
	// std::cerr<<n1<<' '<<n2<<' '<<R1<<' '<<R2<<' '<<sum_k<<' '<<mu<<' '<<sd<<' '<<z<<' '<<pv<<std::endl;
	return pv;
}

// --------------------- Histograms ---------------------------

// goodness-of-fit test of how a distribution fits the data
// http://courses.wcupa.edu/rbove/Berenson/10th%20ed%20CD-ROM%20topics/section12_5.pdf
// df=obs.size()-NumParameters-1; pv=cdf_chisq_q(gof,df); NumParameters=2 for N(mu,sd).
template <typename DIST_T, typename FREQ_T, typename VALUE_T, typename OUT_T>
double goodness_of_fit(const DIST_T& dist,				// distribution to be tested (boost::math::normal)
					   const std::vector<FREQ_T>& data,	// obs freq for (lb,ub]. (int/int64_t)
					   const VALUE_T& lb,				// lower bound, sorted (vector/deque<float/double>)
					   const VALUE_T& ub,				// upper bound, sorted (vector/deque<float/double>)
					   OUT_T& obs,						// obs freq for merged bins (vector/deque<float/double>)
					   OUT_T& exp,						// exp freq for merged bins (vector/deque<float/double>)
					   bool left_truncated=false,		// values <lb[0] already removed from "data"
					   bool right_truncated=false)		// values >ub[l] already removed from "data" (l=last)
{
	FREQ_T ss = 0; // sample size
	for (each_element(data,it)) ss+=*it;
	obs.clear(); exp.clear();
	for (unsigned i=0;i<data.size();++i)
	{
		double e = ss * (boost::math::cdf(dist, ub[i]) - boost::math::cdf(dist, lb[i]));
		double o = data[i];
		if (e<1) // need to merge the tail with the next/previous bin
		{
			if (exp.empty()) // left tail
			{
				for (++i; i<data.size(); ++i)
				{
					e = ss * (boost::math::cdf(dist, ub[i]) - boost::math::cdf(dist, lb[i]));
					if (e>=1) // next bin with e>1
					{
						double nxt_e = 0;
						double nxt_o = 0;
						if (left_truncated) nxt_e = ss * (boost::math::cdf(dist, ub[i]) - boost::math::cdf(dist, lb[0]));
						else				nxt_e = ss *  boost::math::cdf(dist, ub[i]);
						for (unsigned j=0; j<=i; j++) nxt_o += data[j];
						obs.push_back(nxt_o);
						exp.push_back(nxt_e);
						break;
					}
				}
			}
			else // right tail
			{
				double& prv_o = obs.back();
				double& prv_e = exp.back();
				if (right_truncated)	prv_e = ss * (boost::math::cdf(dist, ub.back())-boost::math::cdf(dist, lb[i-1]));
				else					prv_e = ss * (1-boost::math::cdf(dist, lb[i-1]));
				for (; i<data.size(); ++i) prv_o+=data[i];
				break;
			}
		}
		else
		{
			obs.push_back(o);
			exp.push_back(e);
		}
	}
	double gof=0;
	for (unsigned i=0; i<obs.size(); ++i)
	{
		double dif = obs[i]-exp[i];
		gof += dif*dif/exp[i];
	}
	return gof;
}

/*/ to debug goodness_of_fit()
int main (int argc, char * const argv[])
{
	vector<double> lb,ub;
	vector<int> data;
	for (double i=-2; i<2; i+=0.1)
	{
		lb.push_back(i);
		ub.push_back(i+0.1);
		data.push_back(4);
	}
	boost::math::normal dist;
	vector<double> obs,exp;
	goodness_of_fit(dist,data,lb,ub,obs,exp);
	return 0;
} //*/

// Simple method: mean value of a fixed-size window.
// Problem 1) No smoothing at both ends, it's good only if both ends are a sequence of 0s.
// Problem 2) It's better to use a bigger window size when the number of observations is smaller.
// Problem 3) May produce 0 counts, which may be a problem for downstream analysis.
template <typename HISTOGR_T, typename CURVE_T> // vector/deque<..>, vector/deque<..>
void SmoothHistogram(const HISTOGR_T& data, CURVE_T& dest, int ws=3)
{
	dest.clear();
	int ms = (ws - 1) / 2;	// window middle size
	int nb = data.size();	// number of bins
	
	for (int i=0; i<ms; ++i)
		dest.push_back(data[i]);
	
	for (int i=ms; i<nb-ms; ++i)
	{
		double sum = 0;
		for (int j=i-ms; j<=i+ms; ++j)	sum += data[j];
		dest.push_back(sum/ws);
	}
	
	for (int i=nb-ms; i<nb; ++i)
		dest.push_back(data[i]);
}

// Automatically adjust window size by counting the total number of observations of the window.
template <typename HISTOGR_T, typename CURVE_T> // vector/deque<..>, vector/deque<..>
void SmoothHistogram_AdjWinSize(const HISTOGR_T& data, CURVE_T& dest, int ws=3, double min_sum=100)
{
	dest.clear();
	int ms = (ws - 1) / 2;	// window middle size
	int nb = data.size();	// number of bins
	int ms_max = nb -1;		// max ms
	
	for (int i=0; i<nb; ++i) // for (int i=ms; i<nb-ms; ++i)
	{
		bool first_trial = true;
		for (int msc=ms; msc<=ms_max; ++msc)// msc is ms current
		{
			if (first_trial) // window never go out of boundary, msc<ms near boundary
			{
				int jbgn, jend, msu = msc;	// ms used
				if (i-msu>=0) jbgn=i-msu;   else { jbgn=0;  msu=i;      jend=i+msu+1; }
				if (i+msu<nb) jend=i+msu+1; else { jend=nb; msu=nb-i-1; jbgn=i-msu;   }
				double sum=0;
				for (int j=jbgn; j<jend; ++j) sum+=data[j];
				if (sum>min_sum) { dest.push_back(sum/(msu+msu+1)); break; }
				else if (msu!=msc) first_trial=false;
			}
			else // window can go out of boundary, assuming data[external_bins]=0
			{
				int jbgn, jend;
				if (i-msc>=0) jbgn=i-msc;   else { jbgn=0;  }
				if (i+msc<nb) jend=i+msc+1; else { jend=nb; }
				double sum=0;
				for (int j=jbgn; j<jend; ++j) sum+=data[j];
				if (sum) { dest.push_back(sum/(msc+msc+1)); break; }
			}
			// For both of them, data point is always at the middle of a window.
		}
	}
}

// Weight for data[j]=wt[|j-i|]. wt[0]=sum of all weights, because wt for data[i] is always 1.
// If data[j] is out of boundary, suppose its value is 0.

template <typename HISTOGR_T>
double WeightedSum(const HISTOGR_T& data, int i, std::vector<double>& weights, int ms)
{
	double sum = data[i]; // weight=1
	for (int j=1;j<=ms;++j) if (i+j<(int)data.size())	sum += data[i+j] * weights[j];
	for (int j=1;j<=ms;++j) if (i-j>=0)					sum += data[i-j] * weights[j];
	return sum;
}

template <typename HISTOGR_T>
double WeightedAverage(const HISTOGR_T& data, int i, std::vector<double>& weights, int ms)
{
	double sum = data[i];	// weight=1
	double twt = 1;			// total weight
	for (int j=1;j<=ms;++j) if (i+j<(int)data.size()) { sum += data[i+j] * weights[j]; twt += weights[j]; }
	for (int j=1;j<=ms;++j) if (i-j>=0)				  { sum += data[i-j] * weights[j]; twt += weights[j]; }
	return sum / twt;
}

struct SmoothHistogramLoggerType {
	std::vector<int> cs;					// current middle-window-size
	std::vector<int> wp;					// weight pointer
	std::vector< std::vector<double> > wt;	// weight vector
	int ng;							// Number of Grades
	int ms;							// Middle Size
};

// Automatically adjust window size by counting the total number of observations of the window.
// DecreaseConfidence means the reliability of the data points are decreasing. Correspondingly,
// 1) The first i<ms bins are directly copied.
// 2) It has weight vectors with increasing weights for far-away bins. Good for the beginning.
// 3) Window size can only go up, never go down. Good for fluctuations.
template <typename HISTOGR_T, typename CURVE_T> // vector/deque<..>, vector/deque<..>
SmoothHistogramLoggerType SmoothHistogram_DecreaseConfidence(const HISTOGR_T& data, CURVE_T& dest, double min_sum=100, int ws=3, int num_grades=20)
{
	dest.clear();
	int ms = (ws - 1) / 2;	// window middle size
	int nb = data.size();	// number of bins
	int ms_max = nb-1;		// max ms
	
	std::vector< std::vector<double> >  wt(num_grades);	// weights
	for (int i=0; i<num_grades; ++i)	// i is parameter
	{
		wt[i].assign(ms+1,1);			// wt[i][0] is sum of weights, initialized with weight for j=0 that is always 1
		for (int j=1; j<=ms; ++j)
		{
			double w = 1/pow(j+1,0.2*i);// 10/(pow(j+1,i)+9);
			wt[i][j] = w;
			wt[i][0]+= w;
			wt[i][0]+= w;
		}
	}
	
	SmoothHistogramLoggerType log;
	log.cs.assign(data.size(),0);
	log.wp.assign(data.size(),num_grades-1);
	log.wt = wt;
	log.ng = num_grades;
	log.ms = ms;
	
	int WtPtr = num_grades-1;
	int msc = ms;	// current ms
	for (int i=0; i<nb; ++i)
	{
		if (i==0) { dest.push_back(data[i]); continue; }
		if (i<ms) { dest.push_back( WeightedAverage(data,i,wt[WtPtr],i) ); continue; }
		while (msc<=ms_max)
		{
			double sum=0, avg=0;
			if (WtPtr)
			{
				sum = WeightedSum(data,i,wt[WtPtr],ms);
				avg = WeightedAverage(data,i,wt[WtPtr],ms);
			}
			else
			{
				int jbgn = std::max(i-msc,0);
				int jend = std::min(i+msc+1,nb);
				for (int j=jbgn; j<jend; ++j) sum+=data[j];
				avg = sum / (msc+msc+1);
			}
			if (sum>min_sum || msc==ms_max) { dest.push_back(avg); log.cs[i]=msc; log.wp[i]=WtPtr; if (WtPtr) --WtPtr; break; }
			if (WtPtr) --WtPtr; else ++msc;
		}
	}
	return log;
}

template <typename HISTOGR_T, typename CURVE_T> // vector/deque<..>, vector/deque<..>
void SmoothHistogram_DecreaseConfidence(const HISTOGR_T& data, CURVE_T& dest, SmoothHistogramLoggerType& log)
{
	int nb = data.size();
	dest.assign(nb,0);
	for (int i=0; i<nb; ++i)
	{
		int msc = log.cs[i];
		int wpt = log.wp[i];
		if (msc)
		{
			if (wpt)
				dest[i] = WeightedAverage(data,i,log.wt[wpt],log.ms);
			else
			{
				int jbgn = std::max(i-msc,0);
				int jend = std::min(i+msc+1,nb);
				double sum=0; for (int j=jbgn; j<jend; ++j) sum+=data[j];
				dest[i] = sum/(msc+msc+1);
			}
		}
		else
			dest[i] = data[i];
	}
}

template <typename SRC_T, typename DST_T> // both vector/deque<..>
void MakeFraction(const SRC_T& src, DST_T& dst) // Normalize by sum to 1
{
	typedef typename SRC_T::value_type ElementT;
	size_t size;
	ElementT sum;
	array_sum_and_N(src,sum,size);
	if (sum!=1 && sum!=0)
	{
		ElementT factor = 1.0/sum;
		dst.assign(size,0);
		for (unsigned i=0; i<size; ++i) dst[i] = src[i] * factor;
	}
	else
	{
		dst=src;
	}
}

template <typename DST_T> // vector/deque<..>
void MakeFraction(DST_T& src) // Normalize by sum to 1
{
	DST_T dst;
	MakeFraction(src,dst);
	src=dst;
}

template <typename SRC_T, typename DST_T> // both vector/deque<..>
void MakeFractionAgrestiCoull(const SRC_T& src, DST_T& dst)
{
	typedef typename SRC_T::value_type ElementT;
	size_t size;
	ElementT sum;
	array_sum_and_N(src,sum,size);
	if (sum!=1 && sum!=0)
	{
		double factor = 1/(sum+3.841458881296);
		std::vector<double> mid (size,0);
		for (unsigned i=0; i<size; ++i) mid[i] = (src[i]+1.920729440648) * factor;
		MakeFraction(mid,dst);
	}
	else
	{
		dst=src;
	}
}

/*/ to run: <program_name> > smooth_histogram.out
int main (int argc, char * const argv[])
{
	vector<int> src;
	vector<float> dst;
	for (Rows_in_File(in,argv[1],1))
	src.push_back(boost::lexical_cast<int>(in[0]));
	smooth_histogram_nz(src,dst,21);
	print_container(dst,program.outf,'\n',true);
	return 1;
} //*/

// Tested by R: library(robustbase); x<-scan("file"); mc(x)
// But it's much slower than R; most of the time it's making and sorting mcij. Don't know R's trick.
template <typename CONTAINER_T, typename VALUE_T> // CONTAINER_T = vector/deque of VALUE_T
double AdjustedBoxPlot(CONTAINER_T& data, VALUE_T& lf, VALUE_T& uf)
{
	// health check up
	int jend=data.size();
	if (jend<4) exit_error("sample size is too small for AdjustedBoxPlot.");
	assert(std::is_sorted(data.begin(),data.end()));
	
	// make Q1 Q3 IQR med medlb medub medk
	VALUE_T Q1  = quantile(data,0.25);	// prv quantile_Multiset
	VALUE_T Q3  = quantile(data,0.75);	// prv quantile_Multiset
	VALUE_T IQR = Q3-Q1;				// Inter Quartile Range
	VALUE_T med = median(data);			// median, type match to container because I need data[i]==med later
	int		medlb = std::lower_bound(data.begin(),data.end(),med) - data.begin();
	int		medub = std::upper_bound(data.begin(),data.end(),med) - data.begin();
	int		medk  = medub-medlb;		// number of data[i]=med
		
	// make mcij, values to calculate MC
	std::deque<float> mcij;		// double is slow
	for (int i=0;i<medub;++i)
		for (int j=medlb;j<jend;++j)
		{
			if (i==j) continue;
			if (data[i]==med && data[j]==med)
			{
				if		(i+j+1< medk)	mcij.push_back(-1);
				else if (i+j+1==medk)	mcij.push_back(0);
				else					mcij.push_back(1);
			}
			else
			{
				mcij.push_back( ((data[j]-med)-(med-data[i]))/((float)data[j]-data[i]) );
			}
		}

	// calculate MC from mcij
	sort(mcij.begin(),mcij.end());
	float MC = median(mcij);		// MedCouple
	if (MC>=0) { lf = Q1 - 1.5 * exp(-3.5*MC) * IQR ; uf = Q3 + 1.5 * exp(4.0*MC) * IQR; }
	else	   { lf = Q1 - 1.5 * exp(-4.0*MC) * IQR ; uf = Q3 + 1.5 * exp(3.5*MC) * IQR; }
	return MC;
}

// --------------------- p_values or rankings ---------------------

// make ranks (1-based) for the right-values of a map, allow ties
template <typename T1, typename T2>
void rank_up_mapped_values(const std::map<T1,T2>& vMap, std::map<T1,double>& rMap); // v=value r=rank

// make ranks (1-based) for the right-values of a map, allow ties
template <typename T1, typename T2, typename T3>
void num_ge_mapped_values(const std::map<T1,T2>& vMap, std::map<T1,T3>& rMap); // v=value r=rank

// make ranks (1-based) for the right-values of a map, allow ties, smallest value ranks top
template <typename Container>
void rank_up(const Container& vMap, std::vector<double>& rMap); // v=value r=rank

// make ranks (1-based) for the right-values of a map, allow ties, largest value ranks top
template <typename Container>
void rank_dn(const Container& vMap, std::vector<double>& rMap); // v=value r=rank

template <typename T1>
void FDR_from_p(const std::map<T1,double>& vMap, std::map<T1,double>& rMap); // v=value r=rank

// likelihood of exponential distribution at a rate lambda
double lik_exp_dist(double x, double lambda=1, double epsilon=0.001);

// from a set of p-values, return a p-value
template <typename T1>
double CPMA_from_p(const std::map<T1,double>& pval);

// http://en.wikipedia.org/wiki/Fisher's_method
// but also see http://www.burns-stat.com/pages/Working/perfmeasrandport.pdf page 14 for "why you should not use Fisher's method"
// use Stouffer's method instead.
// tested with http://www.ncbi.nlm.nih.gov/CBBresearch/Yu/downloads/CoinedPValues.html
template <typename T1>
double Fisher_from_p(const std::map<T1,double>& pval);

// R code: pnorm(sum(qnorm(x)) / sqrt(length(x))) , where x is the vector of individual p-values.
// this gives the same results: 1 - pnorm(sum(qnorm(1-x)) / sqrt(length(x)))
// Stouffer, S. A., Suchman, E. A , DeVinney, L.C., Star, S.A., Williams, R.M. Jr (1949). 
// Adjustment During Army Life. Princeton, NJ, Princeton University Press. (this is original reference)
template <typename T1>
double Stouffer_from_one_sided_p(const std::map<T1,double>& pval);

// see 2011 Optimally weighted Z-test is a powerful method for combining probabilities in meta analysis.
template <typename T1>
double weighted_Stouffer_from_one_sided_p(const std::map<T1,double>& pval, std::map<T1,double>& wght);

// see 2011 Optimally weighted Z-test is a powerful method for combining probabilities in meta analysis.
// this function has not been tested by comparing with known programs
// different direction of association average out, make sure this is what you want, otherwise use the one_sided_p()
template <typename T1>
double Stouffer_from_two_sided_p(const std::map<T1,double>& stat, const std::map<T1,double>& pval);

// this function has not been tested by comparing with known programs
template <typename T1>
double Stouffer_from_z(const std::map<T1,double>& zMap);

// Draw qq-plot. If type="pdf", gnuplot must be 4.6 or later. If version is too early, set type="postscript".
void qq_plot(std::multiset<double> p_values, const std::string& out_prefix, const int range=10, const std::string& type="pdf", const std::string& program="gnuplot", const bool save_data=false);

// calculate the posterior probability given a log10 Bayes factor
double posterior_given_log10BF(double log10_BF);				// with an uninformative prior (pi=0.5)
double posterior_given_log10BF(double log10_BF, double prior);	// with a provided prior

// --------------------- outlier detection ---------------------

struct outlier_info {
	bool success;
	int num_l_ol;
	int num_r_ol;
	double l_thr;
	double r_thr;
	double mean;
	double median;
	double mad;
	double q25pc;
	double q75pc;
	double IQR;
	double stddev;
	double ih_sigval;
	double gr_t2;
	double gr_sigval;
	double dq_pl;
	double dq_pr;
	outlier_info():success(false),num_l_ol(0),num_r_ol(0),l_thr(-DBL_MAX),r_thr(DBL_MAX),\
	mean(0),median(0),mad(0),q25pc(0),q75pc(0),IQR(0),stddev(0),ih_sigval(0),gr_t2(0),gr_sigval(0),dq_pl(0),dq_pr(0){}
	void init() { success=false; num_l_ol=0; num_r_ol=0; l_thr=-DBL_MAX; r_thr=DBL_MAX; }
};

template <typename T>
double cal_outlier_byDQ(const outlier_info& ol_inf, const T& val); // return p-val

template <typename T>
void num_outlier_byDQ(const std::multiset<T>& values, outlier_info& ol_inf, double alpha=0.05); //Dixon's Q test
//http://www.statistics4u.info/fundstat_eng/cc_outlier_tests_dixon.html#

template <typename T>
double cal_outlier_byBP(const outlier_info& ol_inf, const T& val);

template <typename T>
void num_outlier_byBP(const std::multiset<T>& values, outlier_info& ol_inf, double k=1.5); // box plot

template <typename T>
double cal_outlier_byIH(const outlier_info& ol_inf, const T& val); // return Mi

template <typename T>
void num_outlier_byIH(const std::multiset<T>& values, outlier_info& ol_inf, double alpha=0.000465258); // median absolute deviation
// http://itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
//+-------Iglewicz and Hoaglin's modified Z-sore to detect outliers-----+
//| modified Z-score Mi = 0.6745(Yi-Y~)/MAD                             |
//| where, median absolute deviation (MAD) = median(|Yi-Y~|)            |
//| Y~ = median of the data                                             |
//| data points whose Mi > 3.5 should be labeled as potential outliers, |
//| which has an alpha of 0.000465258 in a two-sided test.              |
//+---------------------------------------------------------------------+

template <typename T>
double cal_outlier_byGr(const outlier_info& ol_inf, const T& val); // Grubbs, return G

template <typename T>
void num_outlier_byGr(const std::multiset<T>& values, outlier_info& ol_inf, double alpha=0.000465258); // took this number from above

// ------------------------------------------- ROC analysis -------------------------------------------

// Reference: AUC Optimization vs. Error Rate Minimization.
// Corinna Cortes and Mehryar Mohri. Advances in Neural Information Processing Systems 16 (NIPS 2003)
// Both container require [] operators, and have the same size(). ScoresContainer should have ::value_type
// Labels=0/1 but not checking. Labels & scores can contain NaN. Higher score is positive.
template <typename LablesContainer,      typename ScoresContainer>
double AUC ( const LablesContainer& labels, const ScoresContainer& scores, bool to_add_uncertainties=false)
{
	if (labels.size()!=scores.size()) exit_error("AUC() requires labels & scores have the same size.");
	typedef typename ScoresContainer::value_type score_t;
	typedef size_t cnt_t;						// not double, count >= instead of mid-point for ties
	using namespace std;
	
	// make cnt,m,n
	cnt_t nr = labels.size();					// number of records including missing values
	cnt_t m=0, n=0;								// number of cases/controls w/o missing value
	map<score_t,pair<cnt_t,cnt_t> > cnt;		// number of observations among cs,ct for each score
	for (size_t i=0; i<nr; ++i)
	{
		if (std::isnan(scores[i]))	continue;	// missing value
		if (std::isnan(labels[i]))	continue;	// missing value
		pair<cnt_t,cnt_t>& r = cnt[scores[i]];
		if (labels[i])	{ ++r.first;  ++m; }	// label=1 means case
		else			{ ++r.second; ++n; }	// label=0 means control
	}

	// add uncertainty => AUC is never 0/1
	if (to_add_uncertainties)
	{
		cnt_t median = (m+n)/2;					// find the score around median
		cnt_t total =0;
		for (each_element(cnt,it))
		{
			total += it->second.first;
			total += it->second.second;
			if (total>=median)
			{
				++it->second.first;
				++it->second.second;
				++m; ++n;
				break;
			}
		}
	}
	
	// calculate AUC
	if (!m || !n) return std::numeric_limits<double>::signaling_NaN();
	double a=0;
	a += cnt.begin()->second.first * cnt.begin()->second.second / 2.0;
	for (each_interval(cnt,i1,i2))
	{
		a += i2->second.first * i1->second.second;
		a += i2->second.first * i2->second.second / 2.0;
		i2->second.second += i1->second.second;
	}
	return (double)a/(m*n);
}

// Used double intead of "ScoresContainer::value_type" because some data type doesn't have it, e.g. Eigen::VectorXf
// Removed adding uncertainties (hence the name _lite)
template <typename LablesContainer, typename ScoresContainer>
double AUC_lite ( const LablesContainer& labels, const ScoresContainer& scores)
{
	using namespace std;
	size_t nr = labels.size();					// number of records including missing values
	size_t m=0, n=0;							// number of cases/controls w/o missing value
	map<double,pair<size_t,size_t> > cnt;		// number of observations among cs,ct for each score
	for (size_t i=0; i<nr; ++i)
	{
		if (std::isnan(scores[i]))	continue;	// missing value
		if (std::isnan(labels[i]))	continue;	// missing value
		pair<size_t,size_t>& r = cnt[scores[i]];
		if (labels[i])	{ ++r.first;  ++m; }	// label=1 means case
		else			{ ++r.second; ++n; }	// label=0 means control
	}
	
	if (!m || !n) return std::numeric_limits<double>::signaling_NaN();
	double a=0;
	a += cnt.begin()->second.first * cnt.begin()->second.second / 2.0;
	for (each_interval(cnt,i1,i2))
	{
		a += i2->second.first * i1->second.second;
		a += i2->second.first * i2->second.second / 2.0;
		i2->second.second += i1->second.second;
	}
	return (double)a/(m*n);
}

// Used double intead of "ScoresContainer::value_type" because some data type doesn't have it, e.g. Eigen::VectorXf
// AUC validated with R::pROC and R:ROCR. Algorithm different from AUC() but result is the same.
extern bool ROC_VERBOSE;
template <typename LablesContainer, typename ScoresContainer>
void ROC ( const LablesContainer& labels, const ScoresContainer& scores, std::set< std::pair<double,double> >& curve, double& optmcc, double& optcut)
{
	double optval = std::numeric_limits<double>::max();
	optmcc =-std::numeric_limits<double>::max();
	optcut = std::numeric_limits<double>::signaling_NaN();
	curve.clear();
	using namespace std;
	size_t nr = labels.size();					// number of records including missing values
	size_t m=0, n=0;							// number of cases/controls w/o missing value
	map<double,pair<size_t,size_t> > cnt;		// number of observations among cs,ct for each score
	for (size_t i=0; i<nr; ++i)
	{
		if (std::isnan(scores[i]))	continue;	// missing value
		if (std::isnan(labels[i]))	continue;	// missing value
		pair<size_t,size_t>& r = cnt[scores[i]];
		if (labels[i])	{ ++r.first;  ++m; }	// label=1 means case
		else			{ ++r.second; ++n; }	// label=0 means control
	}
	
	if (!m || !n) return ;
	for (each_interval_rev(cnt,i1,i2))
	{
		i2->second.first  += i1->second.first;
		i2->second.second += i1->second.second;
	}
	curve.insert(pair<double,double>(0,0));
	curve.insert(pair<double,double>(1,1));
	for (each_element_rev(cnt,it))
	{
		double FP = it->second.second;
		double TP = it->second.first;
		double TN = n - FP;
		double FN = m - TP;
		double one_min_spe = FP / (double)n;
		double sensitivity = TP / (double)m;
		double mcc = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
		// method 1: optimized by sens and spec at the same time
		double val = std::abs(sensitivity+one_min_spe-1); //*/
		/*/ method 2: optimized by MCC
		double val = 1-mcc; //*/
		/*/ method 3: distance from (0,1)
		double val = sqrt(one_min_spe*one_min_spe+(1-sensitivity)*(1-sensitivity)); //*/
		if (val<optval) { optval=val; optcut=it->first; optmcc=mcc; }
		curve.insert(pair<double,double>(one_min_spe,sensitivity));
		if (ROC_VERBOSE)
		{
			std::cerr<<"cut="<<it->first<<" sens="<<sensitivity<<" spec="<<1-one_min_spe<<" mcc="<<mcc<<" |sens-spec|="<<val<<std::endl;
		}
	}
}

inline double AUC_of(const std::set< std::pair<double,double> >& curve, double false_positive_cutoff=1)
{
	if (false_positive_cutoff<=0) exit_error("false_positive_cutoff for AUC_of() must be >0 && <=1");
	if (false_positive_cutoff >1) exit_error("false_positive_cutoff for AUC_of() must be >0 && <=1");
	if (curve.empty()) return std::numeric_limits<double>::signaling_NaN();
	double auc=0;
	for (each_interval(curve,i1,i2))
	{
		const double& x1 = i1->first;
		const double& y1 = i1->second;
		double  x2 = i2->first;
		double  y2 = i2->second;
		if (false_positive_cutoff<=x2)
		{
			y2 = y1 + (false_positive_cutoff-x1)*(y2-y1)/(x2-x1);
			x2 = false_positive_cutoff;
			auc += std::abs(x2-x1) * (y2+y1)/2;
			break;
		}
		else
			auc += std::abs(x2-x1) * (y2+y1)/2;
	}
	return auc;
}

// no output curve, but can calculate partial AUC!
template <typename LablesContainer, typename ScoresContainer>
double partial_AUC ( const LablesContainer& labels, const ScoresContainer& scores, double false_positive_cutoff, double& optmcc, double& optcut)
{
	std::set< std::pair<double,double> > curve;
	ROC (labels, scores, curve, optmcc, optcut);
	return AUC_of(curve,false_positive_cutoff);
}

/*/
sc	l1	l2	l3
.1	0	0	0
.2	0	0	1
.3	0	0	0
.3	0	1	0
.4	1	0	0
.5	0	0	0

 l1 AUC=0.8
 l2 AUC=0.4
 l3 AUC=0.2
 
//*/

// ------------------------------------------- ensemble learning -------------------------------------------

// http://en.wikipedia.org/wiki/Ensemble_learning
// Lables_C: labels for each record, any type that fits fptr
// Scores_C: scores for each record, any type that fits fptr
// Weight_C: weights on each model, vector|deque|..<double> that provides []
#include <type_traits> // for std::is_same (c++11)
template <typename Lables_C, typename Scores_C, typename Weight_C>
void train_bayesian_model_combination(const Lables_C& labels, const Scores_C& scores,									// input data
									  double (*fptr)(const Lables_C&,const Scores_C&,const Weight_C&, std::vector<double>&),	// calculate accuracy|AUC
									  Weight_C& weights,// output weights, inititialy should has the correct size but contents will change
									  int TrainIt=200)	// optional argument. smaller is faster, bigger is more precise, 200 is reasonable
{
	if (!std::is_same<typename Weight_C::value_type,double>::value) exit_error("BMC requires that weights are of double type.");
	if (weights.empty()) exit_error("BMC requires that weights are not empty.");
	size_t nr = labels.size();							// number of records
	size_t nm = weights.size();							// number of models
	weights.assign(nm,0);								// initial weights are 0
	double sum_weight = 0;								// initial weights are 0
	double z = -std::numeric_limits<double>::infinity();// z = -infinity
	Weight_C v(nm);										// temporary weights container
	std::vector<double> results;
	for (int i=0; i<TrainIt; ++i)
	{
		for (auto &w : v) w=-log(rng.uniform_gt0lt1());	// draw from a uniform Dirichlet distribution
		MakeFraction(v);								// Normalize to sum to 1
		double x = (*fptr)(labels,scores,v,results);	// accuracy (0,1) of the ensemble weighted by v
		double log_like = nr*(x*log(x)+(1-x)*log(1-x)); // Use x to estimate log_likelihood_i
		if (log_like > z)								// z is used to maintain numerical stability
		{
			for (size_t m=0; m<nm; ++m)
				weights[m] = weights[m] * exp(z-log_like);
			z = log_like;
		}
		double w = exp(log_like - z);
		for (size_t m=0; m<nm; ++m)
			weights[m] = weights[m] * sum_weight / (sum_weight + w) + w * v[m];
		sum_weight = sum_weight + w;
	}
	MakeFraction(weights);								// Normalize to sum to 1
}

// ------------------------------------------- ensemble learning -------------------------------------------
#include <initializer_list>
class Interpolate {
public:
	std::map<double,double>	reference_data;
	Interpolate() {}
	Interpolate(const std::string& filename) { setup(filename); }
	Interpolate(std::initializer_list< std::pair<double,double> >);
	void setup(std::istream& file);			// 2 columns: x y
	void setup(const std::string& filename);// 2 columns: x y
	double solve(double x,int mode=1);		// return nan if x=nan / ref is empty; mode=0 NA, =1 same_as_last, =2 linear
};

#endif
