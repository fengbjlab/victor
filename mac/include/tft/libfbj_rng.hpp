// body file: libfbj_rng.cpp

#ifndef LIBFBJ_RNG
#define LIBFBJ_RNG

#include <string>
#include <map>
#include <vector>
#include <boost/random.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/math/distributions.hpp>

typedef boost::mt19937 RNGType;

// Usage: just the same as boost:mt19937. Instantiate a distribution -- Instantiate a generator -- use the generator.
// 	  boost::uniform_int<> one_to_six( 1, 6 );
//    boost::variate_generator< RNGType&, boost::uniform_int<> > dice(rng, one_to_six);
//    int n=dice();
class rng_by_boost : public RNGType {
private:
	static int lc;
	unsigned readseed;
	std::string seed_filename;
	boost::posix_time::ptime time_begin, time_end; //begin and end time of program running
	boost::uniform_int<> int_one_to_max;
	boost::variate_generator< RNGType&, boost::uniform_int<> > dice_int;
	boost::uniform_01<RNGType&> zeroone;
	// forbid the following two operations because: although having the same seed, output is diff if othr is already used
	rng_by_boost(const rng_by_boost& othr);					
	rng_by_boost& operator=(const rng_by_boost& othr);
public:
	~rng_by_boost();						// if seed is read from file, write a new seed to the same file, plus the run duration
	rng_by_boost();							// seed = a hash value of microsec, sec since 01/01/1970, process ID
	rng_by_boost(unsigned s);				// seed = s
	rng_by_boost(const std::string& f);		// seed = read from file (skipping whitespaces, the first string should be an unsigned)
	void set_seed();						// seed = a hash value of microsec, sec since 01/01/1970, process ID
	void set_seed(unsigned s);				// seed = s
	void set_seed(const std::string& f);	// seed = read from file (skipping whitespaces, the first string should be an unsigned)
	unsigned long uniform_int_ge0leMax();	// return [0,INT_MAX]
	unsigned long uniform_int_ge0ltN(unsigned long n);	// return [0,n)
	double uniform_gt0lt1();				
	double uniform_ge0lt1();				
	int flip(const double& p);				// return 0/1, p:[0-1]
	void memran(unsigned char* b, size_t l);// randomize memory buffer b with length l
};

extern rng_by_boost rng;

// this can be used in: std::random_shuffle(MyVec.begin(),MyVec.end(),MyRandom);
// however, this is not thread-safe
inline int MyRandom(int i) { return rng.uniform_int_ge0ltN(i); }
std::string random_string(int length);

// Usage:
// rng_normal_distr nd;
// nd.setup_par(mean,sigma);
// nd.gen_num();
class rng_normal_distr {
private:
	boost::normal_distribution<>* dist;
	boost::variate_generator<RNGType&, boost::normal_distribution<> >* vgen;
public:
	~rng_normal_distr();
	rng_normal_distr(); // default = no allocation
	rng_normal_distr(const rng_normal_distr& othr);
	rng_normal_distr& operator=(const rng_normal_distr& othr);
	rng_normal_distr(double m,double s);
	void setup_par(double m,double s);
	double gen_num();
};

// the same as gsl_cdf_ugaussian_P(gsl_ran_ugaussian(r)+offset);
inline double sim_pvalue(double offset)
{
	static rng_normal_distr nd(0,1); // cannot use const, don't know whether it's thread-safe
	static const boost::math::normal dist;
	return boost::math::cdf(dist,nd.gen_num()+offset);
}

// adapted from http://snipplr.com/view/5907/stdrandomshuffle-and-boostrandomhpp/ (original code is wrong [N should be N-1]!)
// to use uniform_int_ge0ltN to shuffle a vector:
//     uniform_int_ge0ltN rnd(myvector.size());
//     std::random_shuffle ( myvector.begin(), myvector.end(), rnd); slightly slower than using built-in random generator:
//     std::random_shuffle ( myvector.begin(), myvector.end() ); // 01.107912 vs 00.921865 seconds
//     however, this is not thread-safe
// to use uniform_int_ge0ltN to create a number:
//     uniform_int_ge0ltN rnd(N); || uniform_int_ge0ltN rnd; rnd.setup_par(N);
//     rnd.gen_num();
class uniform_int_ge0ltN
{
private:
	boost::uniform_int<int>* dist;
	boost::variate_generator< RNGType&, boost::uniform_int<int> >* vgen;
public:
	~uniform_int_ge0ltN();
	uniform_int_ge0ltN();
	uniform_int_ge0ltN(const uniform_int_ge0ltN& othr);
	uniform_int_ge0ltN& operator=(const uniform_int_ge0ltN& othr);
	uniform_int_ge0ltN( const int N );
	void setup_par(int N);
	int gen_num();
	std::ptrdiff_t operator()( std::ptrdiff_t arg );
};

// usage:
// wheeling_select ws(ARGUMENTS); || wheeling_select ws; ws.setup_par(ARGUMENTS); // ARGUMENTS can be removed afterward
// ws.gen_num();
// ARGUMENTS should >=0 but this class provides no checking of validity
class wheeling_select { 
private:
	std::vector<double> cumulative;
	boost::uniform_real<>* dist;
	boost::variate_generator<RNGType&, boost::uniform_real<> >* vgen;
public:
	~wheeling_select();
	wheeling_select();
	wheeling_select(const wheeling_select& othr);
	wheeling_select(const double *f,const int& n);
	wheeling_select(const std::vector<double>& v);
	wheeling_select& operator=(const wheeling_select& othr);
	void setup_par(const double *f,const int& n);
	void setup_par(const std::vector<double>& v);
	int gen_num(); // return [0,n-1]
};

#endif
