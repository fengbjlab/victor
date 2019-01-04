// comp flag: -O1
//     notes: -O2 / -O3 will lead to segmentation fault from mlink
// Compile manually in Mac: g++ -O1 -D_MAKE_PROG -isystem /opt/local/include -L/opt/local/static_lib_libc++ -lptoc linkage51_mlink.cpp -o linkage51_mlink

#include <ptoc/ptoc.h>
#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <climits>
#include <cmath>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "linkage51.hpp"
namespace linkage
{
	MlinkResultType defaultResultType = LOD;
}

namespace {
	
	struct _GEOL {} GEOL;
	inline std::istream& operator>> (std::istream& file, _GEOL& g)
	{
		for (char c=file.get(); c!=EOF; c=file.get())
		{
			if (c=='\n') {	if (file.peek()=='\r') file.get(); break; }
			if (c=='\r') {	if (file.peek()=='\n') file.get(); break; }
		}
		return file;
	}
	inline bool is_EndOfLine(char c) {
		return (c=='\n' || c=='\r' || c==EOF);
	}
	inline bool is_blank_row(std::istream & file) {
		return is_EndOfLine(file.peek());
	}
	inline void exit_error(const std::string& s) { std::cerr<<s<<std::endl; exit(1);	}
	
	std::string fmt(const int& num, int width)
	{
		return str( boost::format("%"+boost::lexical_cast<std::string>(width)+"d") % num );
	}
	
	std::string fmt(const double& num, int width, int prec)
	{
		return str( boost::format("%"+boost::lexical_cast<std::string>(width)+"."+boost::lexical_cast<std::string>(prec)+"f") % num );
	}
	
	const real version = 5.1;		/*PRESENT VERSION OF LINKAGE  10 May 1994*/
	
	/* SOME USER DEFINED CONSTANTS */
	/*THE PROGRAM WILL TELL YOU IF THE FOLLOWING TWO CONSTANTS CAN BE REDUCED*/
	/*IF THE PROGRAM TERMINATES WITH AN ERROR IN RECOMBINATION INCREASE MAXNEED*/
	const int maxneed = 19532;		/*MAXIMUM NUMBER OF RECOMBINATION PROBABILITIES*/
	/*THE FOLLOWING SHOULD BE LARGER THAN MININT*/
	const int maxcensor = 10000;	/*MAXIMUM FOR CENSORING ARRAY*/
	const int maxlocus = 10;		/*MAXIMUM NUMBER OF LOCI */ // fbj 7 => 10
	const int maxrec = maxlocus;	/*MAXIMUM POSSIBLE NUMBER OF RECOMB/MEI*/
	const int maxall = 100;			/*MAXIMUM NUMBER OF ALLELES AT A SINGLE LOCUS 15*/
	const int maxhap = 1280;		/*MAXIMUM NUMBER OF HAPLOTYPES = n1 x n2 x ... WHERE ni = actual # OF ALLELES LOCUS I (not maxall)*/ // fbj 256 => 1280 (even 1536 failed)
	const int maxind = 500000;		/*MAXIMUM NUMBER OF INDIVIDUALS previously 5000 */
	const int maxped = 1000;		/*MAXIMUM NUMBER OF PEDIGREES   previously 500 */
	const int maxchild = 50;		/*MAXIMUM NUMBER OF FULLSIBS IN A SIBSHIP previously 20*/
	const int maxloop = 3;			/*MAXIMUM NUMBER OF LOOPS PER PEDIGREE*/
	const real minfreq = 0.0;
	const int affall = 2;			/*DISEASE ALLELE FOR QUANT. TRAITS OR AFFECTION STATUS*/
	
	const int maxseg = 2^(maxlocus-1);/*2 to the power maxlocus-1*/ // fbj 64 => 2^(maxlocus-1)
	const int maxfem = (maxhap/2)*(maxhap+1);	// MAX # OF JOINT GENOTYPES FOR A FEMALE = (maxhap*(maxhap+1))/2 ; fbj 32896 => (maxhap/2)*(maxhap+1)
	const int maxmal = maxfem;					// MAX # OF JOINT GENOTYPES FOR A MALE   = maxfem (AUTOSOMAL) OR maxhap (SEXLINKED)
	
	/* QUANTITATIVE TRAIT */
	const int maxtrait = 2;			/*MAXIMUM NUMBER OF QUANTITATIVE FACTORS AT A SINGLE LOCUS*/
	const real missval = 0.0;		/*MISSING VALUES FOR QUANTITATIVE TRAITS */
	
	/* AFFECTION STATUS */
	const int missaff = 0;			/*MISSING VALUE FOR AFFECTION STATUS*/
	const int affval = 2;			/*CODE FOR AFFECTED INDIVIDUAL*/
	const int maxliab = 120;			/*MAXIMUM NUMBER OF LIABILITY CLASSES*/ // fbj 10 => 120
	
	/* BINARY (FACTOR UNION) SYSTEM */
	const int maxfact = maxall+1;	/*MAXIMUM NUMBER OF BINARY CODES AT A SINGLE LOCUS*/
	/*Must be at least as large as MAXALL*/
	/* OTHERS */
	const real scale = 2.0;			/*SCALE FACTOR*/
	const real scalemult = 3.0;		/*SCALE WEIGHT FOR EACH LOCUS INCLUDED*/
	const boolean fitmodel = true;	/*TRUE IF ESTIMATING PARAMETERS OTHER THAN REC*/
	const boolean dooutput = false;	/*write outfile.dat*/
	const boolean usespeed = true;	/*read speedfile.dat*/
	const boolean score = true;		/*GIVE LOD SCORES*/
	const boolean byfamily = true;	/*GIVE LIKELIHOODS BY FAMILY*/
	const real zerolike = -1.0E20;	/*FOR INCONSISTENT DATA OR RECOMBINATION */
	const int minint = -32767;		/*MINIMUM ALLOWED int*/
	const real log_10=log(10);
	
	/*GRADIENT APPROXIMATIONS*/
	const boolean approximate = false;	/*Do not change*/
	const real epsilon = 1.0E-6;
	
	typedef struct gennurec* gennuptr;
	struct gennurec {
		array2d<1,maxhap,1,maxhap,int> genenumber;
	};
	typedef struct approxrec* approxpnt;
	struct approxrec {
		array2d<1,1,1,1,boolean> approxarray;            /*changed from 1..maxped,1..maxfem*/
	};
	typedef struct censorrec* censorpnt;
	struct censorrec {
		array1d<minint,maxcensor,boolean> censor;
	};
	typedef array1d<1,maxfem,real> genotype;
	typedef array2d<1,maxtrait,1,maxtrait,real> covmatrix;
	typedef array1d<1,maxlocus,unsigned char> hapvector;
	typedef array2d<0,2,1,2,int> mutarray;
	typedef array1d<1,maxtrait,real> thesemeans;
	typedef array2d<0,maxall,1,maxall,thesemeans> means;
	enum pathway {auto_,mauto,sex,msex, last_pathway};
	enum direction {peelup,peeldown, last_direction};
	
	typedef set binset;
	typedef array1d<1,maxall,binset> phenarray;
	
	typedef struct locusvalues* locuspoint;
	typedef struct phenotype* phenpoint;
	enum locustype {affection,quantitative,binary,null, last_locustype};
	struct locusvalues {
		int nallele,format;
		array1d<1,maxall,real> freq;
		locuspoint privlocus;
		locustype which;
		union {
			struct {
				array4d<0,maxall,1,maxall,0,2,1,maxliab,real> pen;
				int nclass;} s_affection;
			struct {int ntrait;
				means pm;
				covmatrix vmat;
				real det,contrait,conmat;} s_quantitative;
			phenarray allele;
			
		};
	};
	struct phenotype {
		locustype which;
		union {
			binset phenf;
			struct {
				array1d<1,maxtrait,real> x;
				boolean missing; } s_quantitative;
			struct {int aff,liability;} s_affection;
			
		};
	};
	
	typedef struct thisperson* ind;
	typedef struct thisarray* genpoint;
	typedef array1d<1,maxlocus,phenpoint> indphen;
	struct thisarray {
		genotype genarray;
	};
	
	typedef array2d<1,maxall,1,maxall,boolean> possvect;
	typedef array1d<1,maxlocus,possvect> possarray;
	
	typedef struct information* infoptr;
	struct information {
		possarray possible;
	};
	
	struct thisperson {
		int id,ped,inloop;
		ind pa,ma,foff,nextpa,nextma;
		genpoint gen;	// fbj: haplotype data, don't know what the real number means. to get: p->gen->genarray[1-nhap*(nhap+1)/2]
		indphen phen;	// fbj: original genotypes, no phase, not filled. to test: p->phen[locus#]->phenf.has(allele#)
		phenpoint privphen;
		boolean unknown,multi,done,up,male,firstpass;
		infoptr store_;
	};
	
	typedef array1d<1,maxseg,int> subhap;
	typedef array1d<1,maxlocus,real> thetarray;
	typedef struct thetavalues* thetapoint;
	typedef array1d<1,maxneed,real> happrob;
	struct thetavalues {
		thetarray theta;
		happrob segprob;
	};
	
	pathway thispath;
	array1d<1,maxped,boolean> informative;
	array1d<1,maxfem,boolean> rare,risk1,risk2;
	array1d<1,maxhap,boolean> riskmale;
	approxpnt approxstruct;
	censorpnt censorstruct;
	int thisc;
	array1d<1,maxchild,boolean> malechild;
	array1d<1,maxchild,genpoint> thischild;
	int nchild;
	array1d<1,maxfem,unsigned short> seghap1,seghap2;
	array1d<1,maxfem,unsigned short> segstart;
	array1d<1,maxfem,unsigned short> probstart,probend;
	array1d<1,maxlocus,boolean> nohom;
	array1d<1,maxlocus,locuspoint> thislocus;
	array1d<1,maxlocus,int> increment,order;
	/*PEOPLE*/
	array1d<0,maxind,ind> person;
	array1d<1,maxped,ind> proband;
	array3d<1,maxped,1,maxloop,1,2,ind> looppers;
	/*MUTATION */
	array1d<1,maxhap,unsigned short> muthap;
	gennuptr gennustruct;
	/*RECOMBINATION*/
	thetapoint maletheta,femaletheta;
	/*FREQUENCIES */
	genpoint hapfreq;
	/*RISK*/
	unsigned char riskall;
	/*OTHERS*/
	int risksys,mutsys,nlocus,which,lastpriv,i; // fbj rm thissyste, not used
	unsigned short nuhap;
	unsigned short fgeno,mgeno;
	int nuped,totperson;
	real segscale,mutmale,mutfemale,like,alike,tlike,finish,inc1,distratio;
	boolean interfer,disequi,sexlink,risk,sexdif,readfemale,
	mapping,dolod,firstapprox,firsttime,lasttime; // fbj rm inconsistent, not used
	
	int thissystem; // fbj rm j, not used
	real scorevalue,holdtheta;
	std::vector<real> log10like_theta5;
	array1d<1,maxlocus,boolean> zeromale,zerofemale;
	
	struct extra_line_data {
		real t,i,f; // theta, incr, finish
	};
	std::queue<extra_line_data> extra_lines;
	
	real mapfunction(real theta1,real theta2)
	/*User defined function giving recombination between
	 flanking markers as a function of recombination
	 between adjacent markers*/
	{
		real mapfunction_result;
		mapfunction_result=(theta1+theta2)/(1+4*theta1*theta2);
		return mapfunction_result;
	}
	
	
	real getdist(real& theta)
	{
		real getdist_result;
		if (theta<0.5)
			getdist_result=-log(1.0-2.0*theta)/2.0;
		else getdist_result=10.0;
		return getdist_result;
	}
	
	
	real invdist(real& dist)
	{
		real invdist_result;
		if (dist!=10.0)
			invdist_result=(1.0-exp(-2*dist))/2.0;
		else invdist_result=0.5;
		return invdist_result;
	}
	
	
	void invert(covmatrix& m, int n, real& det)
	{
		covmatrix v;
		real val;
		int i,j,k;
		
		
		det=1.0;
		for( i=1; i <= n; i ++)
		{
			val=m[i][i];
			for( k=1; k <= i-1; k ++) val=val-v[k][i]*v[k][i];
			det=det*val;
			v[i][i]=sqrt(val);
			for( j=i+1; j <= n; j ++)
			{
				val=m[i][j];
				for( k=1; k <= i-1; k ++) val=val-v[k][i]*v[k][j];
				v[i][j]=val/v[i][i];
				v[j][i]=0;
			}
		}
		for( i=1; i <= n; i ++)
		{
			m[i][i]=1/v[i][i];
			for( j=i+1; j <= n; j ++)
			{
				val=0;
				for( k=i; k <= j-1; k ++) val=val-m[k][i]*v[k][j];
				m[j][i]=val/v[j][j];
				m[i][j]=0;
			}
		}
		for( i=1; i <= n; i ++)
		{
			for( j=1; j <= i; j ++)
			{
				val=0;
				for( k=j; k <= n; k ++) val=val+m[k][i]*m[k][j];
				v[i][j]=val;
				v[j][i]=val;
			}
		}
		m=v;
	}
	
	void recombination();
	static void recombine(thetarray& theta, happrob& segprob, int& here, real& p4, real& p2, real& p3, real& p1);
	
	
	static void scramble(array1d<1,maxlocus,boolean>& thishet, int& there, int& here, hapvector& hap1, array1d<1,maxseg,hapvector>& thishap1, thetarray& theta, array1d<1,maxseg,hapvector>& thishap2, happrob& segprob)
	{
		int whichhap,start,length,i,j,k,nrec;
		real recval,val;
		
		start=0;
		do { start=start+1; } while (!(thishet[start]));
		length=there-here;
		for( i=2; i <= length; i ++)
		{
			hap1=thishap1[i];
			for( j=1; j <= length; j ++)
			{
				nrec=0;
				val=0.5;
				whichhap=1;
				recval=theta[start];
				for( k=start+1; k <= nlocus; k ++)
				{
					if (! thishet[k])
						recval=recval*(1.0-theta[k])+theta[k]*(1.0-recval);
					else
					{
						if (whichhap==1)
						{
							if (thishap1[j][k]==hap1[k])
								val=val*(1-recval);
							else
							{
								nrec=nrec+1;
								val=val*recval;
								whichhap=2;
							}
						}
						else
						{
							if (thishap2[j][k]==hap1[k])
								val=val*(1-recval);
							else
							{
								nrec=nrec+1;
								val=val*recval;
								whichhap=1;
							}
						}
						recval=theta[k];
					}
				}
				there=there+1;
				if (nrec>maxrec)	segprob[there]=0.0;
				else				segprob[there]=val;
			}
		}
	}
	
	
	static void setrec(real val, int& nhap, array1d<1,maxseg,hapvector>& thishap1, hapvector& hap1, array1d<1,maxseg,hapvector>& thishap2, hapvector& hap2, int& there, happrob& segprob)
	{
		nhap=nhap+1;
		thishap1[nhap]=hap1;
		thishap2[nhap]=hap2;
		there=there+1;
		segprob[there]=val;
	}
	
	
	
	static void dointer(array1d<1,maxlocus,boolean>& thishet, int& nhap, array1d<1,maxseg,hapvector>& thishap1, hapvector& hap1, array1d<1,maxseg,hapvector>& thishap2, hapvector& hap2, int& there, happrob& segprob, thetarray& theta, real& p4, real& p2, real& p3, real& p1)
	{
		int i;
		array1d<1,3,boolean> temphet;
		
		for( i=1; i <= nlocus; i ++)
			temphet[i]=thishet[i];
		
		if (temphet[1] && temphet[2] && ! temphet[3])
			
		{
			setrec(0.5-0.5*theta[1], nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*theta[1], nhap, thishap1, hap1, thishap2, hap2, there, segprob);
			setrec(0.5*theta[1], nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5-0.5*theta[1], nhap, thishap1, hap1, thishap2, hap2, there, segprob);
		}
		else
		{
			if (temphet[3] && temphet[2] && ! temphet[1])
				
			{
				setrec(0.5-0.5*theta[nlocus-1], nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*theta[nlocus-1], nhap, thishap1, hap1, thishap2, hap2, there, segprob);
				setrec(0.5*theta[nlocus-1], nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5-0.5*theta[nlocus-1], nhap, thishap1, hap1, thishap2, hap2, there, segprob);
			}
			else
			{
				if (temphet[3] && temphet[1] && ! temphet[2])
				{
					setrec(0.5-0.5*theta[nlocus], nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*theta[nlocus], nhap, thishap1, hap1, thishap2, hap2, there, segprob);
					setrec(0.5*theta[nlocus], nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5-0.5*theta[nlocus], nhap, thishap1, hap1, thishap2, hap2, there, segprob);
				}
				else
				{	if (temphet[1] && temphet[2] && temphet[3])
				{
					setrec(0.5*p4, nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*p2, nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*p3, nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*p1, nhap, thishap1, hap1, thishap2, hap2, there, segprob);
					setrec(0.5*p2, nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*p4, nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*p1, nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*p3, nhap, thishap1, hap1, thishap2, hap2, there, segprob);
					setrec(0.5*p3, nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*p1, nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*p4, nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*p2, nhap, thishap1, hap1, thishap2, hap2, there, segprob);
					setrec(0.5*p1, nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*p3, nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*p2, nhap, thishap1, hap1, thishap2, hap2, there, segprob); setrec(0.5*p4, nhap, thishap1, hap1, thishap2, hap2, there, segprob);
				}
				else       /*not informative*/
					setrec(0.5, nhap, thishap1, hap1, thishap2, hap2, there, segprob);
				}
			}
		}
	}
	
	
	static void nexthet(int i, real val, boolean inphase, thetarray& theta, hapvector& hap1, hapvector& hap2, array1d<1,maxlocus,boolean>& thishet, int& nhap, array1d<1,maxseg,hapvector>& thishap1, array1d<1,maxseg,hapvector>& thishap2, int& there, happrob& segprob)
	{
		real newval,recval;
		
		recval=theta[i];
		do {
			i=i+1;
			hap1[i]=0;
			hap2[i]=0;
			if ((! thishet[i]) && (! (i==nlocus)))
				recval=recval*(1-theta[i])+(1-recval)*theta[i];
		} while (!((i==nlocus) || thishet[i]));
		if (i!=nlocus)
		{
			if (inphase)	newval=val*(1-recval);
			else			newval=val*recval;
			hap1[i]=1;
			hap2[i]=2;
			nexthet(i,newval,true, theta, hap1, hap2, thishet, nhap, thishap1, thishap2, there, segprob);
			hap2[i]=1;
			hap1[i]=2;
			if (! inphase)	newval=val*(1-recval);
			else			newval=val*recval;
			nexthet(i,newval,false, theta, hap1, hap2, thishet, nhap, thishap1, thishap2, there, segprob);
		}
		else
		{
			if (! thishet[i])  setrec(val, nhap, thishap1, hap1, thishap2, hap2, there, segprob);
			else
			{
				if (inphase)	newval=val*(1-recval);
				else			newval=val*recval;
				hap1[i]=1;
				hap2[i]=2;
				setrec(newval, nhap, thishap1, hap1, thishap2, hap2, there, segprob);
				if (! inphase)	newval=val*(1-recval);
				else			newval=val*recval;
				hap2[i]=1;
				hap1[i]=2;
				setrec(newval, nhap, thishap1, hap1, thishap2, hap2, there, segprob);
			}
		}
	}
	
	
	static void getrecprob(int& nhap, int& there, int& here, array1d<1,maxlocus,boolean>& thishet, hapvector& hap1, hapvector& hap2, array1d<1,maxseg,hapvector>& thishap1, array1d<1,maxseg,hapvector>& thishap2, happrob& segprob, thetarray& theta, real& p4, real& p2, real& p3, real& p1)
	{
		int i;
		
		nhap=0;
		there=here;
		i=0;
		do {
			i=i+1;
			if (thishet[i])
			{
				hap1[i]=1;
				hap2[i]=2;
			}
			else
			{
				hap1[i]=0;
				hap2[i]=0;
			}
		} while (!(thishet[i] || (i==nlocus)));
		if		(i==nlocus)	setrec(0.5, nhap, thishap1, hap1, thishap2, hap2, there, segprob);
		else if (interfer)	dointer(thishet, nhap, thishap1, hap1, thishap2, hap2, there, segprob, theta, p4, p2, p3, p1);
		else				nexthet(i,0.5,true, theta, hap1, hap2, thishet, nhap, thishap1, thishap2, there, segprob);
		if ((nhap>1) && ! interfer)  scramble(thishet, there, here, hap1, thishap1, theta, thishap2, segprob);
		here=there;
	}
	
	
	static void gethet(int& system, array1d<1,maxlocus,boolean>& thishet, int& nhap, int& there, int& here, hapvector& hap1, hapvector& hap2, array1d<1,maxseg,hapvector>& thishap1, array1d<1,maxseg,hapvector>& thishap2, happrob& segprob, thetarray& theta, real& p4, real& p2, real& p3, real& p1)
	{
		int newsystem;
		
		newsystem=system+1;
		thishet[system]=false;
		if (system!=nlocus)	gethet(newsystem, thishet, nhap, there, here, hap1, hap2, thishap1, thishap2, segprob, theta, p4, p2, p3, p1);
		else				getrecprob(nhap, there, here, thishet, hap1, hap2, thishap1, thishap2, segprob, theta, p4, p2, p3, p1);
		thishet[system]=true;
		if (system!=nlocus)	gethet(newsystem, thishet, nhap, there, here, hap1, hap2, thishap1, thishap2, segprob, theta, p4, p2, p3, p1);
		else				getrecprob(nhap, there, here, thishet, hap1, hap2, thishap1, thishap2, segprob, theta, p4, p2, p3, p1);
	}
	
	
	
	static void recombine(thetarray& theta, happrob& segprob, int& here, real& p4, real& p2, real& p3, real& p1)
	{
		int there,nhap,system;
		array1d<1,maxlocus,boolean> thishet;
		hapvector hap1,hap2;
		array1d<1,maxseg,hapvector> thishap1,thishap2;
		
		/*RECOMBINE*/
		here=0;
		system=1;
		gethet(system, thishet, nhap, there, here, hap1, hap2, thishap1, thishap2, segprob, theta, p4, p2, p3, p1);
	}
	
	
	
	static void getfemaletheta()
	{
		real dist;
		int ntheta,i;
		
		
		if (interfer)	ntheta=nlocus;
		else			ntheta=nlocus-1;
		for( i=1; i <= ntheta; i ++)
		{
			dist=getdist(maletheta->theta[i])*distratio;
			femaletheta->theta[i]=invdist(dist);
		}
	}
	
	
	void recombination()
	{
		int i,here;
		real p1,p2,p3,p4;
		thetarray oldtheta;
		
		
		if (interfer)
		{
			thetavalues& with = *maletheta;
			
			if (mapping)  with.theta[nlocus]=mapfunction(with.theta[nlocus-1],with.theta[1]);
			oldtheta=with.theta;
			if ((! mapping) && (! dolod))
			{
				for( i=1; i <= nlocus; i ++)
					oldtheta[i]=1/(1+exp(oldtheta[i]));
				with.theta[1]=oldtheta[1]+oldtheta[nlocus];
				with.theta[nlocus-1]=oldtheta[nlocus-1]+oldtheta[nlocus];
				with.theta[nlocus]=oldtheta[1]+oldtheta[nlocus-1];
				p1=oldtheta[1];
				p2=oldtheta[nlocus-1];
				p3=oldtheta[nlocus];
				p4=1.0-p1-p2-p3;
			}
			else
			{
				p1=(with.theta[1]+with.theta[nlocus]-with.theta[nlocus-1])/2.0;
				p2=(with.theta[nlocus-1]+with.theta[nlocus]-with.theta[1])/2.0;
				p3=(with.theta[nlocus-1]+with.theta[1]-with.theta[nlocus])/2.0;
				p4=1.0-p1-p2-p3;
			}
			recombine(with.theta,with.segprob, here, p4, p2, p3, p1);
		}
		else
			recombine(maletheta->theta,maletheta->segprob, here, p4, p2, p3, p1);
		if (sexdif)
		{
			if (! readfemale)
			{
				if (interfer && ! dolod)
				{
					thetavalues& with = *maletheta;
					
					with.theta[1]=oldtheta[1]+oldtheta[nlocus];
					with.theta[nlocus-1]=oldtheta[nlocus-1]+oldtheta[nlocus];
					with.theta[nlocus]=oldtheta[1]+oldtheta[nlocus-1];
				}
				getfemaletheta();
			}
			if (interfer)
			{
				thetavalues& with = *femaletheta;
				
				if (mapping)  with.theta[nlocus]=mapfunction(with.theta[nlocus-1],with.theta[1]);
				if (readfemale && ! mapping && ! dolod)
				{
					oldtheta=with.theta;
					for( i=1; i <= nlocus; i ++)
						oldtheta[i]=1/(1+exp(oldtheta[i]));
					with.theta[1]=oldtheta[1]+oldtheta[nlocus];
					with.theta[nlocus-1]=oldtheta[nlocus-1]+oldtheta[nlocus];
					with.theta[nlocus]=oldtheta[1]+oldtheta[nlocus-1];
					p1=oldtheta[1];
					p2=oldtheta[nlocus-1];
					p3=oldtheta[nlocus];
					p4=1.0-p1-p2-p3;
				}
				else
				{
					p1=(with.theta[1]+with.theta[nlocus]-with.theta[nlocus-1])/2.0;
					p2=(with.theta[nlocus-1]+with.theta[nlocus]-with.theta[1])/2.0;
					p3=(with.theta[nlocus-1]+with.theta[1]-with.theta[nlocus])/2.0;
					p4=1.0-p1-p2-p3;
				}
				recombine(with.theta,with.segprob, here, p4, p2, p3, p1);
			}
			else
				recombine(femaletheta->theta,femaletheta->segprob, here, p4, p2, p3, p1);
		}
		if (firsttime)
			if (here<maxneed)
				if (dooutput) std::cout << "Maxneed can be reduced to " << here << std::endl;
	}
	
	void getlocations();
	
	
	static boolean checkrare(hapvector& hap1, hapvector& hap2)
	{
		int i;
		boolean check;
		
		boolean checkrare_result;
		check=false;
		for( i=1; i <= nlocus; i ++)
		{
			if (nohom[i])
			{
				locusvalues& with = *thislocus[i];
				if ( (with.freq[hap1[i]]<minfreq) || (with.freq[hap2[i]]<minfreq) )  check=true;
			}
		}
		checkrare_result=check;
		return checkrare_result;
	}
	
	
	static void checkrisk(boolean& riskhet,boolean& riskhom, hapvector& hap1, hapvector& hap2)
	{
		riskhet=false;
		riskhom=false;
		if		( (hap1[risksys]==riskall) && (hap2[risksys]==riskall))
			riskhom=true;
		else if (((hap1[risksys]!=riskall) && (hap2[risksys]==riskall))
				 ||   ((hap2[risksys]!=riskall) && (hap1[risksys]==riskall)))
			riskhet=true;
	}
	
	
	static int gethapn(hapvector& hap)
	{
		int i,n;
		
		int gethapn_result;
		n=1;
		for( i=1; i <= nlocus; i ++)
			n=n+increment[i]*(hap[i]-1);
		gethapn_result=n;
		return gethapn_result;
	}
	
	static void domalerisk(hapvector& hap1);
	
	
	static void setrisk(hapvector& hap1)
	{
		int n;
		
		n=gethapn(hap1);
		if (hap1[risksys]==riskall)  riskmale[n]=true;
		else riskmale[n]=false;
	}
	
	
	
	static void getriskhap(int system, hapvector& hap1)
	{
		int i;
		
		{
			locusvalues& with = *thislocus[system];
			for( i=1; i <= with.nallele; i ++)
			{
				hap1[system]=i;
				if (system!=nlocus)	getriskhap(system+1, hap1);
				else				setrisk(hap1);
			}
		}
	}
	
	
	static void domalerisk(hapvector& hap1)
	{
		getriskhap(1, hap1);
	}
	
	static void domutation(hapvector& hap1);
	
	
	static void setmutation(hapvector& hap1)
	{
		int i,n;
		
		n=gethapn(hap1);
		if (hap1[mutsys]==thislocus[mutsys]->nallele)
			muthap[n]=n;
		else
		{
			i=hap1[mutsys];
			hap1[mutsys]=thislocus[mutsys]->nallele;
			muthap[n]=gethapn(hap1);
			hap1[mutsys]=i;
		}
	}
	
	
	static void getmuthap(int system, hapvector& hap1)
	{
		int i;
		
		
		{
			locusvalues& with = *thislocus[system];
			for( i=1; i <= with.nallele; i ++)
			{
				hap1[system]=i;
				if (system!=nlocus)	getmuthap(system+1, hap1);
				else				setmutation(hap1);
			}
		}
	}
	
	
	static void domutation(hapvector& hap1)
	{
		getmuthap(1, hap1);
	}
	
	
	static void setnumbers(int& ngene, int& here, int& there, int& nseg, hapvector& hap1, hapvector& hap2, boolean& rarepresent, boolean& riskhet, boolean& riskhom, int& thisseg)
	{
		int nhap1,nhap2;
		
		ngene=ngene+1;
		
		segstart[ngene]=here+1;
		probstart[ngene]=there+1;
		probend[ngene]=there+nseg;
		
		there=there+nseg;
		
		nhap1=gethapn(hap1);
		nhap2=gethapn(hap2);
		gennustruct->genenumber[nhap1][nhap2]=ngene;
		gennustruct->genenumber[nhap2][nhap1]=ngene;
		
		if (minfreq!=0.0)
		{
			if (rarepresent)	rare[ngene]=true;
			else				rare[ngene]=false;
		}
		else					rare[ngene]=false;
		if (risk)
		{
			risk1[ngene]=riskhet;
			risk2[ngene]=riskhom;
		}
		
		thisseg=thisseg+1;
		seghap1[thisseg]=nhap1;
		seghap2[thisseg]=nhap2;
	}
	
	
	static void hapscr(int system,int nscramble, array1d<1,maxlocus,boolean>& thishet, int& ngene, int& here, int& there, int& nseg, hapvector& hap1, hapvector& hap2, boolean& rarepresent, boolean& riskhet, boolean& riskhom, int& thisseg)
	{
		int i,j;
		
		if (thishet[system])	nscramble=nscramble+1;
		if (system!=nlocus)		hapscr(system+1,nscramble, thishet, ngene, here, there, nseg, hap1, hap2, rarepresent, riskhet, riskhom, thisseg);
		else					setnumbers(ngene, here, there, nseg, hap1, hap2, rarepresent, riskhet, riskhom, thisseg);
		if (nscramble>1)
		{
			if (hap1[system]!=hap2[system])
			{
				i=hap1[system];
				j=hap2[system];
				hap1[system]=j;
				hap2[system]=i;
				if (system!=nlocus)	hapscr(system+1,nscramble, thishet, ngene, here, there, nseg, hap1, hap2, rarepresent, riskhet, riskhom, thisseg);
				else				setnumbers(ngene, here, there, nseg, hap1, hap2, rarepresent, riskhet, riskhom, thisseg);
				hap1[system]=i;
				hap2[system]=j;
			}
		}
	}
	
	
	
	static void sethap(int system, array1d<1,maxlocus,boolean>& thishet, hapvector& hap1, hapvector& hap2, boolean& rarepresent, boolean& riskhet, boolean& riskhom, int& there, int& start, int& thisseg, int& here, int& ngene, int& nseg)
	{
		int i,j;
		
		{
			locusvalues& with = *thislocus[system];
			if (thishet[system])
				for( i=1; i <= with.nallele-1; i ++)
				{
					hap1[system]=i;
					for( j=i+1; j <= with.nallele; j ++)
					{
						hap2[system]=j;
						if (system!=nlocus)  sethap(system+1, thishet, hap1, hap2, rarepresent, riskhet, riskhom, there, start, thisseg, here, ngene, nseg);
						else
						{
							rarepresent=checkrare(hap1, hap2);
							if (risk) checkrisk(riskhet,riskhom, hap1, hap2);
							there=start;
							thisseg=here;
							hapscr(1,0, thishet, ngene, here, there, nseg, hap1, hap2, rarepresent, riskhet, riskhom, thisseg);
							here=here+nseg;
						}
					}
				}
			else
				for( i=1; i <= with.nallele; i ++)
				{
					hap1[system]=i;
					hap2[system]=i;
					if (system!=nlocus)  sethap(system+1, thishet, hap1, hap2, rarepresent, riskhet, riskhom, there, start, thisseg, here, ngene, nseg);
					else
					{
						rarepresent=checkrare(hap1, hap2);
						if (risk) checkrisk(riskhet,riskhom, hap1, hap2);
						thisseg=here;
						there=start;
						hapscr(1,0, thishet, ngene, here, there, nseg, hap1, hap2, rarepresent, riskhet, riskhom, thisseg);
						here=here+nseg;
					}
				}
		}
	}
	
	
	
	static void starthap(int& nseg, int& nhet, array1d<1,maxlocus,boolean>& thishet, hapvector& hap1, hapvector& hap2, boolean& rarepresent, boolean& riskhet, boolean& riskhom, int& there, int& start, int& thisseg, int& here, int& ngene)
	{
		int i;
		
		nseg=1;
		for( i=2; i <= nhet; i ++) nseg=nseg*2;
		sethap(1, thishet, hap1, hap2, rarepresent, riskhet, riskhom, there, start, thisseg, here, ngene, nseg);
		start=there;
	}
	
	
	
	static void gethet1(int system, array1d<1,maxlocus,boolean>& thishet, int& nseg, int& nhet, hapvector& hap1, hapvector& hap2, boolean& rarepresent, boolean& riskhet, boolean& riskhom, int& there, int& start, int& thisseg, int& here, int& ngene)
	{
		thishet[system]=false;
		if (system!=nlocus)  gethet1(system+1, thishet, nseg, nhet, hap1, hap2, rarepresent, riskhet, riskhom, there, start, thisseg, here, ngene);
		else starthap(nseg, nhet, thishet, hap1, hap2, rarepresent, riskhet, riskhom, there, start, thisseg, here, ngene);
		thishet[system]=true;
		nhet=nhet+1;
		if (system!=nlocus)  gethet1(system+1, thishet, nseg, nhet, hap1, hap2, rarepresent, riskhet, riskhom, there, start, thisseg, here, ngene);
		else starthap(nseg, nhet, thishet, hap1, hap2, rarepresent, riskhet, riskhom, there, start, thisseg, here, ngene);
		nhet=nhet-1;
	}
	
	void getlocations()
	{
		int ngene,nseg,here,there,start,nhet,
		thisseg;
		boolean rarepresent,riskhom,riskhet;
		hapvector hap1,hap2;
		array1d<1,maxlocus,boolean> thishet;
		
		nhet=0;
		here=0;
		there=0;
		ngene=0;
		start=0;
		gethet1(1, thishet, nseg, nhet, hap1, hap2, rarepresent, riskhet, riskhom, there, start, thisseg, here, ngene);
		if (mutsys!=0)  domutation(hap1);
		if (sexlink && risk)  domalerisk(hap1);
	}
	
	
	static void respond()
	/*Include file to LINKAGE programs*/
	{
		std::cerr << "  Press Enter key to continue or Ctrl-C to abort" << std::endl;
		input >> NL;
	}
	
	
	static void inputerror(int nerror,int par1,int par2)
	{
		// this line added on 2015-07-09
		throw 200+nerror;
		
		std::cerr << "Fatal error detected in procedure inputdata" << std::endl;
		switch (nerror) {
			case 0: std::cerr << "Number of loci " << fmt(par1,2) << " exceeds the constant maxlocus" << std::endl; break;
			case 1: std::cerr << "Number of loci read " << fmt(par1,2) << ". Less than minimum of 1" << std::endl; break;
			case 2: std::cerr << "Error detected reading loci order. Locus number " << fmt(par2,2) << " in position " << fmt(par1,2) << " exceeds number of loci" << std::endl; break;
			case 3: std::cerr << "Error detected reading loci order. Illegal locus number " << fmt(par2,2) << " in position " << fmt(par1,2) << std::endl; break;
			case 4: std::cerr << "Error detected reading loci order. Locus number repeated in positions " << fmt(par1,2) << " and " << fmt(par2,2) << std::endl; break;
			case 5: std::cerr << "Error detected reading locus description. Illegal locus type " << fmt(par2,2) << " for locus " << fmt(par1,2) << std::endl; break;
			case 6: std::cerr << "Error detected reading locus description for system " << fmt(par1,2) << ". Number of alleles  " << fmt(par1,2) << " exceeds maxall" << std::endl; break;
			case 7: std::cerr << "Error detected reading locus description for system " << fmt(par1,2) << ". Illegal number of alleles  " << fmt(par2,2) << std::endl; break;
			case 8: std::cerr << "Error detected reading locus description for system " << fmt(par1,2) << ". Number of factors  " << fmt(par2,2) << " exceeds maxfact" << std::endl; break;
			case 9: std::cerr << "Error detected reading locus description for system " << fmt(par1,2) << ". Illegal number of factors  " << fmt(par2,2) << std::endl; break;
			case 10: std::cerr << "Error detected reading locus description for system " << fmt(par1,2) << ". Alleles not codominant" << std::endl; break;
			case 11: std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Illegal code for sex " << fmt(par2,2) << std::endl; break;
			case 12: std::cerr << "Error detected reading pedigree record at pedigree " << fmt(par1,2) << ". > maxped" << std::endl; break;
			case 13: std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". > maxind" << std::endl; break;
			case 14: std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Illegal binary factor code " << fmt(par2,2) << std::endl; break;
			case 15: std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". No allelic pair for genotype" << std::endl; break;
			case 16: std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Allele number " << fmt(par2,2) << " exceeds maxall" << std::endl; break;
			case 17: std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Illegal allele number " << fmt(par2,2) << std::endl; break;
			case 18: std::cerr << "Number of systems after factorization (" << fmt(par1,3) << ") exceeds maxsystem" << std::endl; break;
			case 19: std::cerr << "Number of systems after factorization (" << fmt(par1,3) << ") less than minimum of 1" << std::endl; break;
			case 20: std::cerr << "Number of recombination types (" << fmt(par1,3) << ") exceeds maxrectype" << std::endl; break;
			case 21: std::cerr << "Number of recombination types (" << fmt(par1,3) << ") less than minimum of 1" << std::endl; break;
			case 22: std::cerr << "End of file detected in tempdat by procedure readthg before all data found" << std::endl; break;
			case 23: std::cerr << "Error detected reading iterated locus in datafile. Value (" << fmt(par1,3) << ") greater than nlocus" << std::endl; break;
			case 24: std::cerr << "Error detected reading iterated locus in datafile. Illegal value (" << fmt(par1,3) << ')' << std::endl; break;
			case 25: std::cerr << "Number of iterated parameters greater then maxn" << std::endl; break;
			case 26: std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Liability class (" << fmt(par2,2) << ") exceeds nclass" << std::endl; break;
			case 27: std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Illegal liability class (" << fmt(par2,2) << ')' << std::endl; break;
			case 28: std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Liability classes (" << fmt(par2,3) << ") exceed maxliab" << std::endl; break;
			case 29: std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Illegal number of liability classes (" << fmt(par2,3) << ')' << std::endl; break;
			case 30: std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Penetrance out of range" << std::endl; break;
			case 31: std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Number of traits (" << fmt(par2,3) << ") exceeds maxtrait" << std::endl; break;
			case 32: std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". No. traits (" << fmt(par2,3) << ')' << std::endl; break;
			case 33: std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Variance must be positive" << std::endl; break;
			case 34: std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Illegal variance multiplier" << std::endl; break;
			case 35: std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Risk allele " << fmt(par2,3) << ") > nallele" << std::endl; break;
			case 36: std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Illegal risk allele (" << fmt(par2,3) << ')' << std::endl; break;
			case 37: std::cerr << "Error detected reading datafile. Risk locus " << fmt(par2,3) << ") exceeds nlocus" << std::endl; break;
			case 38: std::cerr << "Error detected reading datafile. Illegal value for risk locus " << fmt(par2,3) << ')' << std::endl; break;
			case 39: std::cerr << "Error detected reading datafile. Mutation locus " << fmt(par2,3) << ") exceeds nlocus" << std::endl; break;
			case 40: std::cerr << "Error detected reading datafile. Illegal value for mutation locus " << fmt(par2,3) << ')' << std::endl; break;
			case 41: std::cerr << "Error detected reading datafile. Linkage disequilbirium is not allowed with this program" << std::endl; break;
			case 42: std::cerr << "Locus " << fmt(par1,5) << " in lod score list exceeds nlocus " << fmt(par2,5) << std::endl; break;
			case 43: std::cerr << "Illegal locus number " << fmt(par1,5) << " in lod score list" << std::endl; break;
			case 44: std::cerr << "Actual number of loops in pedigree data exceeds constant MAXLOOP = " << fmt(par1,1) << std::endl; break; /*changed*/
		}
		respond(); /*changed*/
	}
	
	
	static void inputwarning(int nwarning,int par1,int par2)
	{
		// this line added on 2015-07-09
		throw 100+nwarning;

		std::cerr << "Warning number from procedure inputdata" << std::endl;
		switch (nwarning) {
			case 0: std::cerr << "Illegal sex difference parameter " << fmt(par1,2) << " Parameter should be 0, 1, or 2" << std::endl; break;
			case 1: std::cerr << "Illegal interference parameter " << fmt(par1,2) << " Lack of interference assumed" << std::endl; break;
			case 2: std::cerr << "Illegal sex difference parameter " << fmt(par1,2) << " Parameter must be 0 with sex-linked data" << std::endl; break;
			case 3: std::cerr << "Non-standard affection status" << fmt(par2,4) << " interpreted as normal in pedigree record" << fmt(par1,5) << std::endl; break;
		}
		respond(); /*changed*/
	}
	
	
	static void readspeed(std::istream& speedfile)
	{
		int i,a,b,sys;
		char ch;
		
		while (! is_blank_row(speedfile))
		{
			ch = speedfile.get(); // previously "speedfile >> ch;" but works only with ptoc::text, not with istream
			if ((ch=='i') || (ch=='I'))
			{
				speedfile >> ch >> i;
				person[i]->unknown=true;
				person[i]->store_ = new information;
				{
					information& with = *person[i]->store_;
					for( sys=1; sys <= nlocus; sys ++)
						for( a=1; a <= maxall; a ++)
							for( b=1; b <= maxall; b ++)
								with.possible[sys][a][b]=false;
				}
			}
			else
			{
				information& with = *person[i]->store_;
				speedfile >> sys >> a >> b;
				if (sys<=nlocus)  with.possible[order[sys]][a][b]=true;
			}
			speedfile >> GEOL;
		}
	}
	
	static void readped(std::istream& ipedfile, array1d<1,maxlocus,int>& nfactor);
	static void getphenotype(std::istream& ipedfile, ind& p, array1d<1,maxlocus,int>& nfactor);
	
	
	static void readbin(std::istream& ipedfile, phenpoint& phen, locuspoint thislocus, array1d<1,maxlocus,int>& nfactor, int& system, ind& p)
	{
		int i,j;
		{
			phenotype& with = *phen;
			with.which=binary;
			with.phenf=set::of(eos);
			for( i=1; i <= nfactor[system]; i ++)
			{
				ipedfile >> j;
				if ((j!=0) && (j!=1))
					inputerror(14,p->id,j);
				if (j==1)
					with.phenf=with.phenf+set::of(i, eos);
			}
		}
	}
	
	
	static void readnumber(std::istream& ipedfile, phenpoint& phen, locuspoint thislocus, ind& p)
	{
		int i,j;
		
		{
			phenotype& with = *phen;
			with.which=binary;
			with.phenf=set::of(eos);
			for( i=1; i <= 2; i ++)
			{
				ipedfile >> j;
				if (j>maxall)
					inputerror(16,p->id,j);
				if (j<0)
					inputerror(17,p->id,j);
				if (j!=0)
					with.phenf=with.phenf+set::of(j, eos);
			}
		}
	}
	
	
	
	static void readaff(std::istream& ipedfile, phenpoint& phen, locuspoint thislocus, ind& p)
	{
		int thisval;
		{
			phenotype& with = *phen;
			with.which=affection;
			ipedfile >> thisval;
			if (thisval==missaff)
				with.s_affection.aff=0;
			else
			{	if (thisval==affval)
				with.s_affection.aff=2;
			else
			{
				if (thisval!=1)  inputwarning(3,p->id,thisval);
				with.s_affection.aff=1;
			}
			}
			if (thislocus->s_affection.nclass==1)	with.s_affection.liability=1;
			else									ipedfile >> with.s_affection.liability;
			if (with.s_affection.liability>thislocus->s_affection.nclass)  inputerror(26,p->id,with.s_affection.liability);
			if (with.s_affection.liability<=0)  inputerror(27,p->id,with.s_affection.liability);
		}
	}
	
	
	static void readquan(std::istream& ipedfile, phenpoint& phen, locuspoint thislocus, ind& p)
	{
		int i;
		real xval;
		
		{
			phenotype& with = *phen;
			if ((! sexlink) || (! p->male))
			{
				with.which=quantitative;
				{
					locusvalues& with1 = *thislocus;
					for( i=1; i <= with1.s_quantitative.ntrait; i ++)
						ipedfile >> with.s_quantitative.x[i];
					with.s_quantitative.missing=true;
					for( i=1; i <= with1.s_quantitative.ntrait; i ++)
						if (with.s_quantitative.x[i]!=missval)
							with.s_quantitative.missing=false;
				}
			}
			else
			{
				with.which=affection;
				ipedfile >> xval;
				{
					locusvalues& with1 = *thislocus;
					if (xval==missval)
						with.s_affection.aff=missaff;
					else
					{
						if (xval==affall)	with.s_affection.aff=affall;
						else				with.s_affection.aff=-11;
					}
					with.s_affection.liability=1;
					for( i=2; i <= with1.s_quantitative.ntrait; i ++) ipedfile >> xval;
				}
			}
		}
	}
	
	
	static void getphenotype(std::istream& ipedfile, ind& p, array1d<1,maxlocus,int>& nfactor)
	{
		int thisread,system;
		{
			thisperson& with = *p;
			for( thisread=1; thisread <= nlocus; thisread ++)
			{
				system=order[thisread];
				with.phen[system]=nil;
				if (thislocus[system]->which!=null)
					with.phen[system] = new phenotype;
				switch (thislocus[system]->which) {
					case quantitative:	readquan(ipedfile, with.phen[system],thislocus[system], p); break;
					case affection:		readaff(ipedfile, with.phen[system],thislocus[system], p); break;
					case binary:
						if (thislocus[system]->format==3) readnumber(ipedfile, with.phen[system],thislocus[system], p);
						else readbin(ipedfile, with.phen[system],thislocus[system], nfactor, system, p);
						break;
					case null: exit_error("thislocus[system]->which = null (fbj)");
					case last_locustype: exit_error("thislocus[system]->which = last_locustype (fbj)");
				}
				if (lastpriv==system)
				{
					with.privphen = new phenotype;
					switch (thislocus[system]->privlocus->which) {
						case quantitative: readquan(ipedfile, with.privphen,thislocus[system]->privlocus, p); break;
						case affection: readaff(ipedfile, with.privphen,thislocus[system]->privlocus, p); break;
						case binary: exit_error("thislocus[system]->which = binary (fbj)");
						case null: exit_error("thislocus[system]->which = null (fbj)");
						case last_locustype: exit_error("thislocus[system]->which = last_locustype (fbj)");
					}
				}
			}
		}
	}
	
	
	static void getinformative(array1d<1,maxped,int>& startped, array1d<1,maxped,int>& endped)
	{
		int i,j,k,l,m,count,nchild;
		ind child;
		
		if (fitmodel || risk)  for( i=1; i <= nuped; i ++) informative[i]=true;
		else
			for( i=1; i <= nuped; i ++)
			{
				informative[i]=false;
				for( j=startped[i]; j <= endped[i]; j ++)
					if (person[j]->foff!=nil)
					{
						nchild=0;
						child=person[j]->foff;
						do {
							nchild=nchild+1;
							if (person[j]->male)	child=child->nextpa;
							else					child=child->nextma;
						} while (!(child==nil));
						count=0;
						if ((nchild>1) || ((nchild==1) && (person[j]->pa!=nil)))
							for( k=1; k <= nlocus; k ++)
							{
								phenotype& with = *person[j]->phen[k];
								{
									locusvalues& with1 = *thislocus[k];
									if (with1.which!=binary)  count=count+1;
									else
									{
										if (with.phenf==set::of(eos))  count=count+1;
										else
										{
											l=0;
											for( m=1; m <= with1.nallele; m ++)
												if (with.phenf.has(m))  l=l+1;
											if (l>1)  count=count+1;
										}
									}
								}
							}
						if (count>1)  informative[i]=true;
					}
			}
	}
	
	
	static void getind(std::istream& ipedfile, int& id, int& sequence)
	
	{               /*getind*/
		ipedfile >> id;
		if (id!=0)
		{
			id=id+sequence;
			if (id>maxind)
				inputerror(13,id,id);
			if (person[id]==nil)  person[id] = new thisperson;
		}
	}                 /*getind*/
	
	
	static void multimarriage(ind& p)
	
	{
		ind q,child;
		
		/*multimarriage*/
		if (p->foff!=nil)
		{
			thisperson& with = *p;
			
			if (with.male)	q=with.foff->ma;
			else			q=with.foff->pa;
			child=with.foff;
			p->multi=false;
			do {
				if (with.male)
				{
					with.multi=q==child->ma;
					child=child->nextpa;
				}
				else {
					with.multi=q==child->pa;
					child=child->nextma;
				}
			} while (!((child==nil) || (with.multi)));
		}
		else p->multi=false;
	}                 /*multimarriage*/
	
	
	static void readped(std::istream& ipedfile, array1d<1,maxlocus,int>& nfactor)
	{
		int i,newid,sex,profield,newped,sequence,nuperson,thisone,thisped;
		array1d<1,maxped,int> startped,endped;
		ind holdloop;
		
		
		/*readped*/
		for( i=0; i <= maxind; i ++) person[i]=nil;
		sequence=0;
		nuperson=0;
		nuped=1;
		for( i=1; i <= maxloop; i ++)
		{
			looppers[nuped][i][1]=nil;
			looppers[nuped][i][2]=nil;
		}
		proband[nuped]=nil;
		ipedfile >> newped;
		thisped=newped;
		startped[1]=1;
		while (! is_blank_row(ipedfile))
		{
			nuperson=nuperson+1;
			getind(ipedfile, thisone, sequence);
			if (proband[nuped]==nil)  proband[nuped]=person[thisone];
			{
				thisperson& with = *person[thisone];
				
				with.ped=thisped;
				with.id=thisone;
				getind(ipedfile, newid, sequence);
				with.pa=person[newid];
				getind(ipedfile, newid, sequence);
				with.ma=person[newid];
				getind(ipedfile, newid, sequence);
				with.foff=person[newid];
				getind(ipedfile, newid, sequence);
				with.nextpa=person[newid];
				getind(ipedfile, newid, sequence);
				with.nextma=person[newid];
				ipedfile >> sex;
				if ((sex!=1) && (sex!=2))
					inputerror(11,with.id,sex);
				if (sex==1)	with.male=true;
				else		with.male=false;
				with.unknown=false;
				with.inloop=0;
			}
			ipedfile >> profield;
			if ((profield - 1) > maxloop)  inputerror(44,maxloop,0);     /*changed*/
			if (profield==1)
				proband[nuped]=person[thisone];
			else if ((profield>1) && (profield-1<=maxloop))
			{	if (looppers[nuped][profield-1][2]==nil)
				looppers[nuped][profield-1][2]=person[thisone];
			else
				looppers[nuped][profield-1][1]=person[thisone];
			}
			getphenotype(ipedfile,person[thisone], nfactor);
			ipedfile >> GEOL;
			if (! is_blank_row(ipedfile))  ipedfile >> newped;
			if (newped!=thisped)
			{
				sequence=nuperson+sequence;
				endped[nuped]=sequence;
				nuperson=0;
				nuped=nuped+1;
				if (nuped>maxped)
					inputerror(12,newped,nuped);
				startped[nuped]=sequence+1;
				for( i=1; i <= maxloop; i ++)
				{
					looppers[nuped][i][1]=nil;
					looppers[nuped][i][2]=nil;
				}
				proband[nuped]=nil;
				thisped=newped;
			}
		}
		totperson=sequence+nuperson;
		endped[nuped]=totperson;
		for( newped=1; newped <= nuped; newped ++)
		{
			if ((looppers[newped][1][2]!=nil) && (looppers[newped][1][1]==nil))
				looppers[newped][1][1]=proband[newped];
			for( i=1; i <= maxloop; i ++)
				if (looppers[newped][i][1]==nil)
					looppers[newped][i][2]=nil;
				else
				{
					looppers[newped][i][1]->inloop=i;
					looppers[newped][i][2]->inloop=i;
					if ((looppers[newped][i][1]->pa==nil) && (looppers[newped][i][2]->pa!=nil))
					{
						holdloop=looppers[newped][i][1];
						looppers[newped][i][1]=looppers[newped][i][2];
						looppers[newped][i][2]=holdloop;
					}
				}
		}
		for( thisone=1; thisone <= totperson; thisone ++) multimarriage(person[thisone]);
		getinformative(startped, endped);
	}                 /*readped*/
	
	
	static void getbin(std::istream& datafile, locuspoint& locus, int& system, array1d<1,maxlocus,int>& nfactor)
	{
		int i,j,k;
		
		datafile >> nfactor[system] >> GEOL;
		if (nfactor[system]>maxfact)  inputerror(8,system,nfactor[system]);
		if (nfactor[system]<=0)  inputerror(9,system,nfactor[system]);
		{
			locusvalues& with = *locus;
			
			for( i=1; i <= with.nallele; i ++) with.allele[i]=set::of(eos);
			for( i=1; i <= with.nallele; i ++)
				for( j=1; j <= nfactor[system]; j ++)
				{
					datafile >> k;
					if (k==1)
						with.allele[i]=with.allele[i]+set::of(j, eos);
				}
		}
		datafile >> GEOL;
	}
	
	
	
	static void getnumber(locuspoint& locus, int& system)
	{
		int i;
		
		{
			locusvalues& with = *locus;
			for( i=1; i <= with.nallele; i ++)
				with.allele[i]=set::of(i, eos);
		}
	}
	
	// these variables are used by ch_pen_xxx()
	double maxpen=0, minpen=1;
	int affection_locus=0;
	array4d<0,maxall,1,maxall,0,2,1,maxliab,real> bkup_pen;
	
	static void getpen(std::istream& datafile, locuspoint& locus, int& system)
	{
		int i,j,l;
		{
			locusvalues& with = *locus;
			datafile >> with.s_affection.nclass >> GEOL;
			if (with.s_affection.nclass>maxliab)	inputerror(28,system,with.s_affection.nclass);
			if (with.s_affection.nclass<=0)			inputerror(29,system,with.s_affection.nclass);
			for( l=1; l <= with.s_affection.nclass; l ++)
			{
				for( i=1; i <= with.nallele; i ++)
					for( j=i; j <= with.nallele; j ++)
					{
						double p;
						datafile >> p;
						if ((i!=1||j!=1) && p>maxpen) maxpen=p;
						if ((i!=1||j!=1) && p<minpen) minpen=p;
						if (p<0 || p>1.0)  inputerror(30,system,system);
						with.s_affection.pen[i][j][2][l] = p;
						with.s_affection.pen[i][j][1][l] = 1-p;
						with.s_affection.pen[i][j][0][l] = 1.0;
						with.s_affection.pen[j][i][2][l] = p;
						with.s_affection.pen[j][i][1][l] = 1-p;
						with.s_affection.pen[j][i][0][l] = 1.0;
					}
				datafile >> GEOL;
				for( i=1; i <= with.nallele; i ++) with.s_affection.pen[0][i][0][l]=1.0;
				if (sexlink)
				{
					for( i=1; i <= with.nallele; i ++)
					{
						double p;
						datafile >> p;
						if ((i!=1||j!=1) && p>maxpen) maxpen = p;
						if ((i!=1||j!=1) && p<minpen) minpen = p;
						if (p<0 || p>1.0)  inputerror(30,system,system);
						with.s_affection.pen[0][i][2][l] = p;
						with.s_affection.pen[0][i][1][l] = 1-p;
					}
					datafile >> GEOL;
				}
			}
			affection_locus = system;
			bkup_pen = with.s_affection.pen;
		}
	}
	
	void ch_pen_dom(double scale)
	{
		int i,j,l;
		{
			locusvalues& with = *(thislocus[affection_locus]);
			for( l=1; l <= with.s_affection.nclass; l ++)
			{
				for( i=1; i <= with.nallele; i ++)
					for( j=i; j <= with.nallele; j ++)
					{
						if (i==1 && j==1) continue;
						double p =  bkup_pen[i][j][2][l] * scale;
						with.s_affection.pen[i][j][2][l] = p;
						with.s_affection.pen[i][j][1][l] = 1-p;
						//						with.s_affection.pen[i][j][0][l] = 1.0; // constant, no need to assign values again
						with.s_affection.pen[j][i][2][l] = p;
						with.s_affection.pen[j][i][1][l] = 1-p;
						//						with.s_affection.pen[j][i][0][l] = 1.0; // constant, no need to assign values again
					}
				//				for( i=1; i <= with.nallele; i ++) with.s_affection.pen[0][i][0][l]=1.0; // constant, no need to assign values again
				if (sexlink)
				{
					for( i=1; i <= with.nallele; i ++)
					{
						double p =  bkup_pen[0][i][2][l] * scale;
						with.s_affection.pen[0][i][2][l] = p;
						with.s_affection.pen[0][i][1][l] = 1-p;
					}
				}
			}
		}
	}
	
	void ch_pen_dom(double ph, double pen)
	{
		double p;
		int i,j,l;
		{
			locusvalues& with = *(thislocus[affection_locus]);
			for( l=1; l <= with.s_affection.nclass; l ++)
			{
				for( i=1; i <= with.nallele; i ++)
					for( j=i; j <= with.nallele; j ++)
					{
						if (i==1 && j==1)	p = ph;
						else				p = pen;
						with.s_affection.pen[i][j][2][l] = p;
						with.s_affection.pen[i][j][1][l] = 1-p;
						//						with.s_affection.pen[i][j][0][l] = 1.0; // constant, no need to assign values again
						with.s_affection.pen[j][i][2][l] = p;
						with.s_affection.pen[j][i][1][l] = 1-p;
						//						with.s_affection.pen[j][i][0][l] = 1.0; // constant, no need to assign values again
					}
				//				for( i=1; i <= with.nallele; i ++) with.s_affection.pen[0][i][0][l]=1.0; // constant, no need to assign values again
				if (sexlink)
				{
					for( i=1; i <= with.nallele; i ++)
					{
						if (i==1)	p = ph;
						else		p = pen;
						with.s_affection.pen[0][i][2][l] = p;
						with.s_affection.pen[0][i][1][l] = 1-p;
					}
				}
			}
		}
	}
	
	static void getquan(std::istream& datafile, locuspoint& locus, boolean privelege, int& system, int& nupriv)
	/*Get information of a quantitative locus. Privelege says whether it is
	 a priveleged locus or not*/
	{
		int i,j,k;
		
		{
			locusvalues& with = *locus;
			datafile >> with.s_quantitative.ntrait >> GEOL;
			if (with.s_quantitative.ntrait>maxtrait)  inputerror(31,system,with.s_quantitative.ntrait);
			if (with.s_quantitative.ntrait<=0)  inputerror(32,system,with.s_affection.nclass);
			for( i=1; i <= with.s_quantitative.ntrait; i ++)
				for( j=1; j <= with.nallele; j ++)
					for( k=j; k <= with.nallele; k ++)
					{
						datafile >> with.s_quantitative.pm[j][k][i];
						with.s_quantitative.pm[k][j][i]=with.s_quantitative.pm[j][k][i];
					}
			datafile >> GEOL;
			if ((! privelege) || (nupriv==lastpriv))
				
			{
				for( i=1; i <= with.s_quantitative.ntrait; i ++)
					for( j=i; j <= with.s_quantitative.ntrait; j ++)
					{
						datafile >> with.s_quantitative.vmat[i][j];
						if ((i==j) && (with.s_quantitative.vmat[i][j]<=0.0))  inputerror(33,system,system);
						with.s_quantitative.vmat[j][i]=with.s_quantitative.vmat[i][j];
					}
				datafile >> GEOL;
				invert(with.s_quantitative.vmat,with.s_quantitative.ntrait,with.s_quantitative.det);
				with.s_quantitative.det=1/sqrt(with.s_quantitative.det);
				datafile >> with.s_quantitative.conmat >> GEOL;
				if (with.s_quantitative.conmat<=0)
					inputerror(34,system,system);
				with.s_quantitative.conmat=1/with.s_quantitative.conmat;
				with.s_quantitative.contrait=1.0;
				for( i=1; i <= with.s_quantitative.ntrait; i ++) with.s_quantitative.contrait=with.s_quantitative.contrait*with.s_quantitative.conmat;
				with.s_quantitative.contrait=sqrt(with.s_quantitative.contrait);
			}
		}
	}
	
	
	static void getlocus(std::istream& datafile, int system, int& whichtype, int& nupriv, array1d<1,maxlocus,int>& nfactor)
	{
		int j;
		
		thislocus[system] = new locusvalues;
		{
			locusvalues& with = *thislocus[system];
			with.privlocus=nil;
			datafile >> whichtype >> with.nallele;
			if ((whichtype<0) && (whichtype>4))
				inputerror(5,system,whichtype);
			if (with.nallele>maxall)
				inputerror(6,system,with.nallele);
			if (with.nallele<=0)
				inputerror(7,system,with.nallele);
			switch (whichtype) {
				case 0 : with.which=quantitative;	break;
				case 1 : with.which=affection;		break;
				case 2 : with.which=binary;			break;
				case 3 : with.which=binary;			break;
			}
			with.format=whichtype;
			if (lastpriv==0)
				datafile >> GEOL;
			else
			{
				datafile >> whichtype >> GEOL;
				if ((whichtype==0) || (whichtype==1))
				{
					nupriv=nupriv+1;
					with.privlocus = new locusvalues;
					with.privlocus->nallele=with.nallele;
					{
						locusvalues& with1 = *with.privlocus;
						switch (whichtype) {
							case 0 : with1.which=quantitative;	break;
							case 1 : with1.which=affection;		break;
						}
					}
				}
			}
			if (! disequi)
			{
				for( j=1; j <= with.nallele; j ++) datafile >> with.freq[j];
				datafile >> GEOL;
			}
			switch (with.which) {
				case binary :
					if (with.format==3)  getnumber(thislocus[system],system);
					else getbin(datafile,thislocus[system],system, nfactor); break;
				case affection :		getpen(datafile,thislocus[system], system); break;
				case quantitative :		getquan(datafile,thislocus[system],false, system, nupriv); break;
				case null:				exit_error("thislocus[system]->which = null (fbj)");
				case last_locustype:	exit_error("thislocus[system]->which = last_locustype (fbj)");
			}
			if (with.privlocus!=nil)
				switch (with.privlocus->which) {
					case affection :	getpen(datafile,with.privlocus, system); break;
					case quantitative :	getquan(datafile,with.privlocus,true, system, nupriv); break;
					case null:			exit_error("thislocus[system]->which = null (fbj)");
					case last_locustype:exit_error("thislocus[system]->which = last_locustype (fbj)");
					case binary:		exit_error("thislocus[system]->which = binary (fbj)");
				}
			if ((nupriv==lastpriv) && (lastpriv!=0))
				lastpriv=system;
			if (risk)
			{
				if (system==risksys) datafile >> riskall >> GEOL;
				if (riskall>thislocus[system]->nallele)  inputerror(35,system,riskall);
				//		if (riskall<0)  inputerror(36,system,riskall);
			}
		}
	}
	
	
	static void gettheta(std::istream& datafile, thetapoint& sex)
	{
		thetarray oldtheta;
		int i;
		
		sex = new thetavalues;
		{
			thetavalues& with = *sex;
			
			if ((sex==maletheta) || readfemale)
				
			{
				for( i=1; i <= nlocus-1; i ++) datafile >> with.theta[i];
				if (interfer && (! mapping))
					datafile >> with.theta[nlocus];
				datafile >> GEOL;
			}
			else datafile >> distratio >> GEOL;
			/*     FOR j:=1 TO maxneed DO segprob[j]:=0.0;*/
			if (interfer && ! mapping)
			{
				oldtheta=with.theta;
				with.theta[1]=(oldtheta[1]+oldtheta[nlocus]-oldtheta[nlocus-1])/2.0;
				with.theta[nlocus-1]=(oldtheta[nlocus-1]+oldtheta[nlocus]-oldtheta[1])/2.0;
				with.theta[nlocus]=(oldtheta[1]+oldtheta[nlocus-1]-oldtheta[nlocus])/2.0;
				for( i=1; i <= nlocus; i ++)
				{
					if (with.theta[i]>0.0)
					{
						if (with.theta[i]<0.999)	with.theta[i]=log(1.0/with.theta[i]-1.0);
						else						with.theta[i]=-6.9;/*=ln(1/0.999-1.0)*/
					}
					else with.theta[i]=9.21;/*=ln(1/0.0001-1.0)*/
				}
			}
		}
	}
	
	
	static void readloci(std::istream& datafile, array1d<1,maxlocus,int>& nfactor)
	{
		int i,j,coupling,autosomal,independent,difference,whichtype,
		nupriv;
		
		lastpriv=0;
		datafile >> nlocus >> risksys >> autosomal >> GEOL;
		/*Replace the above line with the next when using epistasis*/
		/*readln(datafile,nlocus,risksys,autosomal,lastpriv);*/
		if (nlocus>maxlocus)  inputerror(0,nlocus,nlocus);
		if (nlocus<=0)  inputerror(1,nlocus,nlocus);
		if (risksys>maxlocus)  inputerror(37,risksys,risksys);
		if (risksys<0)  inputerror(38,risksys,risksys);
		risk=risksys!=0;
		sexlink= autosomal==1;
		datafile >> mutsys >> mutmale >> mutfemale >> coupling >> GEOL;
		if (mutsys>maxlocus)  inputerror(39,mutsys,mutsys);
		if (mutsys<0)  inputerror(40,mutsys,mutsys);
		if (coupling==1)	disequi=true;
		else				disequi=false;
		if (disequi)  hapfreq = new thisarray;
		for( i=1; i <= nlocus; i ++)
		{
			datafile >> j;
			if (j>nlocus)	inputerror(2,i,j);
			if (j<=0)		inputerror(3,i,j);
			order[j]=i;
		}
		for( i=1; i <= nlocus; i ++)
			for( j=1; j <= i-1; j ++)
				if (order[i]==order[j])  inputerror(4,i,j);
		datafile >> GEOL;
		if (mutsys!=0)
			mutsys=order[mutsys];
		if (risksys!=0)
			risksys=order[risksys];
		nupriv=0;
		for( i=1; i <= nlocus; i ++) getlocus(datafile, order[i], whichtype, nupriv, nfactor);
		increment[nlocus]=1;
		for( i=nlocus-1; i >= 1; i --)
			increment[i]=increment[i+1]*thislocus[i+1]->nallele;
		fgeno=1;
		for( j=1; j <= nlocus; j ++) fgeno=thislocus[j]->nallele*fgeno;
		mgeno=fgeno;
		nuhap=fgeno;
		for( i=1; i <= nlocus; i ++)
			nohom[i]=false;
		if (disequi)
		{
			for( i=1; i <= mgeno; i ++) datafile >> hapfreq->genarray[i];
			datafile >> GEOL;
		}
		else
		{
			for( i=1; i <= nlocus; i ++)
			{ locusvalues& with = *thislocus[i];
				if ((with.which==affection) || (with.which==quantitative))
					if (with.freq[affall]<minfreq)  nohom[i]=true;}
		}
		fgeno=(fgeno*(fgeno+1)) / 2;
		if (! sexlink)
			
			mgeno=fgeno;
		datafile >> difference;
		if ((difference<0) || (difference>2))
		{
			inputwarning(0,difference,difference);
			difference=0;
		}
		sexdif= difference!=0;
		readfemale=difference==2;
		datafile >> independent >> GEOL;
		if ((independent<0) || (independent>2))
		{
			inputwarning(1,independent,independent);
			independent=0;
		}
		interfer= independent!=0;
		mapping= independent==2;
		gettheta(datafile,maletheta);
		if (sexdif)
			
			gettheta(datafile,femaletheta);
		else femaletheta=maletheta;
		if (sexlink && sexdif)
		{
			inputwarning(2,difference,difference);
			sexdif=false;
			readfemale=false;
		}
		if (! sexlink)
			if (mutsys==0)  thispath=auto_;
			else thispath=mauto;
			else if (mutsys==0)  thispath=sex;
			else thispath=msex;
		segscale=scale;
		for( i=1; i <= nlocus; i ++) segscale=segscale*scalemult;
		
		datafile >> which >> inc1 >> finish >> GEOL;
		
		while (! is_blank_row(datafile))
		{
			extra_line_data d;
			datafile >> d.t >> d.i >> d.f >> GEOL;
			if (d.i!=0) extra_lines.push(d);
			else break;
		}
		
		if (dooutput)
		{
			std::cout << "LINKAGE (V" << fmt(version,3,1) << ") WITH" << fmt(nlocus,3) << "-POINT";
			if (sexlink)	std::cout << " SEXLINKED DATA" << std::endl;
			else			std::cout << " AUTOSOMAL DATA" << std::endl;
			std::cout << " ORDER OF LOCI: ";
			for( i=1; i <= nlocus; i ++)
			{
				thissystem=1;
				while (order[thissystem]!=i)  thissystem=thissystem+1;
				std::cout << fmt(thissystem,3);
			}
			std::cout << std::endl;
		}
	}
	
	void inputdata(std::istream& datafile, std::istream& ipedfile, std::istream& speedfile)
	{
		array1d<1,maxlocus,int> nfactor;
		readloci(datafile,nfactor);	// read datafile.dat
		readped(ipedfile,nfactor);	// read ipedfile.dat
		if (usespeed) readspeed(speedfile); // read speedfile.dat
	}
	
	void likelihood(int thisped,ind proband);
	
	
	static void cleanup(ind& p)
	{
		{
			thisperson& with = *p;
			delete with.gen;
			with.gen=nil;
		}
	}
	
	static void getvect(ind p);
	
	
	static real quanfun(phenpoint phen,locuspoint thislocus,int i,int j,thesemeans& mean)
	{
		real val;
		int it,jt;
		
		real quanfun_result;
		{
			phenotype& with = *phen;
			{
				locusvalues& with1 = *thislocus;
				val=1.0;
				if (! with.s_quantitative.missing)
				{
					val=0;
					for( it=1; it <= with1.s_quantitative.ntrait; it ++)
						for( jt=1; jt <= with1.s_quantitative.ntrait; jt ++)
							if (i==j)
								val=val+with1.s_quantitative.vmat[it][jt]*(with.s_quantitative.x[jt]-mean[jt])*(with.s_quantitative.x[it]-mean[it]);
							else
								val=val+with1.s_quantitative.conmat*with1.s_quantitative.vmat[it][jt]*(with.s_quantitative.x[jt]-mean[jt])*(with.s_quantitative.x[it]-mean[it]);
					val=with1.s_quantitative.det*exp(-val*0.5);
					if (i!=j)  val=val*with1.s_quantitative.contrait;
				}
				quanfun_result=val;
			}
		}
		return quanfun_result;
	}
	
	
	static void getval(int syste,int i,int j,real& val, ind& p)
	{
		{
			locusvalues& with = *thislocus[syste];
			switch (with.which) {
				case quantitative:val=val*quanfun(p->phen[syste],thislocus[syste],i,j,with.s_quantitative.pm[i][j]); break;
				case affection:{ phenotype& with1 = *p->phen[syste];  val=val*with.s_affection.pen[i][j][with1.s_affection.aff][with1.s_affection.liability];}
					break;
				case null: exit_error("thislocus[system]->which = null (fbj)");
				case last_locustype: exit_error("thislocus[system]->which = last_locustype (fbj)");
				case binary: exit_error("thislocus[system]->which = binary (fbj)");
			}
		}
		//		std::cerr<<p->ped<<' '<<p->id<<' '<<syste<<' '<<i<<' '<<j;
		//		if (p->phen[2]!=nil) std::cerr<<' '<<(p->phen[2]->phenf.has(1)?2:1)<<' '<<(p->phen[2]->phenf.has(2)?2:1);
		//		if (p->gen!=nil) for (int i=1;i<=10;i++) std::cerr<<' '<<p->gen->genarray[i];
		//		std::cerr<<std::endl;
	}
	
	static void setval(real val, ind& p, hapvector& hap1, hapvector& hap2);
	
	
	static void prior(ind& p, real& val, hapvector& hap1, hapvector& hap2, int& nhap1, int& nhap2)
	{
		int i;
		
		if (! disequi)
			if (sexlink && p->male)
				for( i=1; i <= nlocus; i ++)
				{	locusvalues& with = *thislocus[i];
					val=val*with.freq[hap1[i]]; }
			else
			{
				for( i=1; i <= nlocus; i ++)
				{	locusvalues& with = *thislocus[i];
					val=val*with.freq[hap1[i]]*with.freq[hap2[i]]; }
				if (nhap1!=nhap2)  val=2.0*val;
			}
			else
			{
				thisarray& with = *hapfreq;
				if (sexlink && p->male)
					val=val*with.genarray[nhap1];
				else {
					val=val*with.genarray[nhap1]*with.genarray[nhap2];
					if (nhap1!=nhap2)  val=2.0*val;
				}
			}
		val=val*segscale;
	}
	
	
	static void setval(real val, ind& p, hapvector& hap1, hapvector& hap2)
	{
		int here,count,nhap1,nhap2,i;
		
		count=1;
		nhap1=1;
		nhap2=1;
		if (sexlink && p->male)
		{
			for( i=1; i <= nlocus; i ++) nhap1=nhap1+increment[i]*(hap1[i]-1);
			here=nhap1;
		}
		else
		{
			gennurec& with = *gennustruct;
			for( i=1; i <= nlocus; i ++)
			{
				nhap1=nhap1+increment[i]*(hap1[i]-1);
				nhap2=nhap2+increment[i]*(hap2[i]-1);
				if (hap1[i]!=hap2[i])  count=2*count;
				here=with.genenumber[nhap1][nhap2];
			}
		}
		if (p->pa==nil)  prior(p, val, hap1, hap2, nhap1, nhap2);
		{
			thisarray& with = *p->gen;
			if (! disequi)
			{
				if (count!=1)  count=count / 2;
				for( i=1; i <= count; i ++)
				{
					with.genarray[here]=val;
					here=here+1;
				}
			}
			else with.genarray[here]=val;
		}
	}
	
	static void getgene(int syste,real val, ind& p, hapvector& hap1, hapvector& hap2);
	
	
	static void facmale(ind& p, int& syste, hapvector& hap1, real& val, hapvector& hap2)
	{
		int j;
		
		/*facmale*/
		{
			thisperson& with = *p;
			{
				locusvalues& with1 = *thislocus[syste];
				for( j=1; j <= with1.nallele; j ++)
					if ((with.phen[syste]->phenf==with1.allele[j]) || (with.phen[syste]->phenf==set::of(eos)))
					{
						hap1[syste]=j;
						if (syste!=nlocus) getgene(syste+1,val, p, hap1, hap2);
						else setval(val, p, hap1, hap2);
					}
			}
		}
	}       /*facmale*/
	
	
	static void affmale(int& syste, real& newval, real& val, ind& p, hapvector& hap1, hapvector& hap2)
	{
		int j;
		
		/*affmale*/
		{
			locusvalues& with = *thislocus[syste];
			for( j=1; j <= with.nallele; j ++)
			{
				newval=val;
				getval(syste,0,j,newval, p);
				hap1[syste]=j;
				if (newval!=0.0)
				{	if (syste!=nlocus) getgene(syste+1,newval, p, hap1, hap2);
				else setval(newval, p, hap1, hap2);
				}
			}
		}
	}       /*affmale*/
	
	
	static void quanmale(ind& p, int& syste, real& newval, real& val, hapvector& hap1, hapvector& hap2)
	{
		int j;
		
		/*quanmale*/
		{
			thisperson& with = *p;
			{
				locusvalues& with1 = *thislocus[syste];
				if ((with.phen[syste]->s_affection.aff==affall) || (with.phen[syste]->s_affection.aff==missaff))
				{
					newval=val;
					hap1[syste]=affall;
					if (newval!=0.0)
					{	if (syste!=nlocus) getgene(syste+1,newval, p, hap1, hap2);
					else setval(newval, p, hap1, hap2);
					}
				}
				if ((with.phen[syste]->s_affection.aff!=affall) || (with.phen[syste]->s_affection.aff==missaff))
					for( j=1; j <= with1.nallele; j ++)
						if (j!=affall)
						{
							newval=val;
							hap1[syste]=j;
							if (newval!=0.0)
							{
								if (syste!=nlocus) getgene(syste+1,newval, p, hap1, hap2);
								else setval(newval, p, hap1, hap2);
							}
						}
			}
		}
	}       /*quanmale*/
	
	
	static void fac(ind& p, int& syste, hapvector& hap1, hapvector& hap2, real& val)
	{
		int i,j;
		
		/*fac*/
		{
			thisperson& with = *p;
			{
				locusvalues& with1 = *thislocus[syste];
				for( i=1; i <= with1.nallele; i ++)
				{
					hap1[syste]=i;
					for( j=i; j <= with1.nallele; j ++)
						if ((with.phen[syste]->phenf==(with1.allele[i]+with1.allele[j])) || (with.phen[syste]->phenf==set::of(eos)))
						{
							hap2[syste]=j;
							if (syste!=nlocus)
								getgene(syste+1,val, p, hap1, hap2);
							else setval(val, p, hap1, hap2);
						}
				}
			}
		}
		if (disequi)
		{
			thisperson& with = *p;
			{
				locusvalues& with1 = *thislocus[syste];
				for( i=1; i <= with1.nallele; i ++)
				{
					hap2[syste]=i;
					for( j=i; j <= with1.nallele; j ++)
						if ((with.phen[syste]->phenf==(with1.allele[i]+with1.allele[j])) || (with.phen[syste]->phenf==set::of(eos)))
						{
							hap1[syste]=j;
							if (syste!=nlocus)
								getgene(syste+1,val, p, hap1, hap2);
							else setval(val, p, hap1, hap2);
						}
				}
			}
		}
	}       /*fac*/
	
	
	
	static void aff(int& syste, hapvector& hap1, hapvector& hap2, real& newval, real& val, ind& p)
	/*Used with an affection status phenotype or when
	 thislocus[syste]^which is null*/
	{
		int i,j;
		
		{
			locusvalues& with = *thislocus[syste];
			for( i=1; i <= with.nallele; i ++)
			{
				hap1[syste]=i;
				for( j=i; j <= with.nallele; j ++)
					if ((! nohom[syste]) || (i!=affall) || (j!=affall))
					{
						hap2[syste]=j;
						newval=val;
						getval(syste,i,j,newval, p);
						if (newval!=0.0)
						{	if (syste!=nlocus)  getgene(syste+1,newval, p, hap1, hap2);
						else setval(newval, p, hap1, hap2);
						}
					}
			}
		}
		if (disequi)
		{
			locusvalues& with = *thislocus[syste];
			for( i=1; i <= with.nallele; i ++)
			{
				hap2[syste]=i;
				for( j=i; j <= with.nallele; j ++)
					if ((! nohom[syste]) || (i!=affall) || (j!=affall))
					{
						hap1[syste]=j;
						newval=val;
						getval(syste,i,j,newval, p);
						if (newval!=0.0)
						{
							if (syste!=nlocus)  getgene(syste+1,newval, p, hap1, hap2);
							else setval(newval, p, hap1, hap2);
						}
					}
			}
		}
	}
	
	
	
	
	static void quanval(int& syste, hapvector& hap1, hapvector& hap2, real& newval, real& val, ind& p)
	/*Uses this only when thislocus[syste]^.which is not null*/
	{
		int i,j;
		
		{
			locusvalues& with = *thislocus[syste];
			for( i=1; i <= with.nallele; i ++)
			{
				hap1[syste]=i;
				for( j=i; j <= with.nallele; j ++)
					if ((! nohom[syste]) || (i!=affall) || (j!=affall))
					{
						hap2[syste]=j;
						newval=val;
						getval(syste,i,j,newval, p);
						if (newval!=0.0)
						{
							if (syste!=nlocus)  getgene(syste+1,newval, p, hap1, hap2);
							else setval(newval, p, hap1, hap2);
						}
					}
			}
		}
		if (disequi)
		{
			locusvalues& with = *thislocus[syste];
			for( i=1; i <= with.nallele; i ++)
			{
				hap2[syste]=i;
				for( j=i; j <= with.nallele; j ++)
					if ((! nohom[syste]) || (i!=affall) || (j!=affall))
					{
						hap1[syste]=j;
						newval=val;
						getval(syste,i,j,newval, p);
						if (newval!=0.0)
						{
							if (syste!=nlocus)  getgene(syste+1,newval, p, hap1, hap2);
							else setval(newval, p, hap1, hap2);
						}
					}
			}
		}
	}
	
	
	static void getgene(int syste,real val, ind& p, hapvector& hap1, hapvector& hap2)
	{
		real newval;
		
		{
			locusvalues& with = *thislocus[syste];
			if (sexlink && p->male)
				switch (with.which) {
					case binary:facmale(p, syste, hap1, val, hap2); break;
					case affection:affmale(syste, newval, val, p, hap1, hap2); break;
					case quantitative:quanmale(p, syste, newval, val, hap1, hap2); break;
					case null:
						if (with.privlocus->which==affection)
							affmale(syste, newval, val, p, hap1, hap2);
						else quanmale(p, syste, newval, val, hap1, hap2);
						break;
					case last_locustype: exit_error("thislocus[system]->which = last_locustype (fbj)");
				}
			else
				switch (with.which) {
					case binary:fac(p, syste, hap1, hap2, val); break;
					case affection:aff(syste, hap1, hap2, newval, val, p); break;
					case quantitative:quanval(syste, hap1, hap2, newval, val, p); break;
					case null:aff(syste, hap1, hap2, newval, val, p); break;
					case last_locustype: exit_error("thislocus[system]->which = last_locustype (fbj)");
				}
		}
	}
	
	static void ugetgene(int syste,real val, ind& p, hapvector& hap1, hapvector& hap2);
	
	static void facmale1(ind& p, int& syste, hapvector& hap1, real& val, hapvector& hap2)
	{
		int j;
		
		/*facmale*/
		{
			thisperson& with = *p;
			{ information& with1 = *p->store_;
				{ locusvalues& with2 = *thislocus[syste];
					for( j=1; j <= with2.nallele; j ++)
						if ((with.phen[syste]->phenf==with2.allele[j]) || (with.phen[syste]->phenf==set::of(eos)))
							if (with1.possible[syste][1][j])
							{
								hap1[syste]=j;
								if (syste!=nlocus)
									ugetgene(syste+1,val, p, hap1, hap2);
								else setval(val, p, hap1, hap2);
							}
				}
			}
		}
	}       /*facmale*/
	
	
	static void affmale1(ind& p, int& syste, real& newval, real& val, hapvector& hap1, hapvector& hap2)
	{
		int j;
		
		{
			information& with = *p->store_;
			{
				locusvalues& with1 = *thislocus[syste];
				for( j=1; j <= with1.nallele; j ++)
					if (with.possible[syste][1][j])
					{
						newval=val;
						getval(syste,0,j,newval, p);
						hap1[syste]=j;
						if (newval!=0.0)
						{
							if (syste!=nlocus)  ugetgene(syste+1,newval, p, hap1, hap2);
							else setval(newval, p, hap1, hap2);
						}
					}
			}
		}
	}
	
	
	static void quanmale1(ind& p, int& syste, real& newval, real& val, hapvector& hap1, hapvector& hap2)
	{
		int j;
		
		
		{
			thisperson& with = *p;
			{
				information& with1 = *p->store_;
				{
					locusvalues& with2 = *thislocus[syste];
					
					if ((with.phen[syste]->s_affection.aff==affall) || (with.phen[syste]->s_affection.aff==missaff))
						
						if (with1.possible[syste][1][affall])
						{
							newval=val;
							hap1[syste]=affall;
							if (newval!=0.0)
							{
								if (syste!=nlocus)  ugetgene(syste+1,newval, p, hap1, hap2);
								else setval(newval, p, hap1, hap2);
							}
						}
					{
						information& with3 = *p->store_;
						if ((with.phen[syste]->s_affection.aff!=affall)
							|| (with.phen[syste]->s_affection.aff==missaff))
							
							for( j=1; j <= with2.nallele; j ++)
								if (j!=affall)
									
									if (with3.possible[syste][1][j])
									{
										newval=val;
										hap1[syste]=j;
										if (newval!=0.0)
										{
											if (syste!=nlocus)  ugetgene(syste+1,newval, p, hap1, hap2);
											else setval(newval, p, hap1, hap2);
										}
									}
					}
				}
			}
		}
	}
	
	
	static void fac1(ind& p, int& syste, hapvector& hap1, hapvector& hap2, real& val)
	{
		int i,j;
		
		
		{
			thisperson& with = *p;
			{
				information& with1 = *p->store_;
				{
					locusvalues& with2 = *thislocus[syste];
					for( i=1; i <= with2.nallele; i ++)
					{
						hap1[syste]=i;
						for( j=i; j <= with2.nallele; j ++)
							if ((with.phen[syste]->phenf==(with2.allele[i]+with2.allele[j]))
								|| (with.phen[syste]->phenf==set::of(eos)))
								
								if (with1.possible[syste][i][j])
								{
									hap2[syste]=j;
									if (syste!=nlocus)  ugetgene(syste+1,val, p, hap1, hap2);
									else setval(val, p, hap1, hap2);
								}
					}
				}
			}
		}
		if (disequi)
		{
			thisperson& with = *p;
			{
				information& with1 = *p->store_;
				{
					locusvalues& with2 = *thislocus[syste];
					for( i=1; i <= with2.nallele; i ++)
					{
						hap2[syste]=i;
						for( j=i; j <= with2.nallele; j ++)
							if ((with.phen[syste]->phenf==(with2.allele[i]+with2.allele[j]))
								|| (with.phen[syste]->phenf==set::of(eos)))
								
								if (with1.possible[syste][i][j])
								{
									hap1[syste]=j;
									if (syste!=nlocus)  ugetgene(syste+1,val, p, hap1, hap2);
									else setval(val, p, hap1, hap2);
								}
					}
				}
			}
		}
	}
	
	
	static void aff1(ind& p, int& syste, hapvector& hap1, hapvector& hap2, real& newval, real& val)
	/*Used with an affection status phenotype or when
	 thislocus[syste]^which is null*/
	{
		int i,j;
		
		{
			//fbj	thisperson& with = *p;
			{
				information& with1 = *p->store_;
				{
					locusvalues& with2 = *thislocus[syste];
					for( i=1; i <= with2.nallele; i ++)
					{
						hap1[syste]=i;
						for( j=i; j <= with2.nallele; j ++)
							if ((! nohom[syste]) || (i!=affall) || (j!=affall))
								if (with1.possible[syste][i][j])
								{
									hap2[syste]=j;
									newval=val;
									getval(syste,i,j,newval, p);
									if (newval!=0.0)
									{	if (syste!=nlocus) ugetgene(syste+1,newval, p, hap1, hap2);
									else setval(newval, p, hap1, hap2);
									}
								}
					}
				}
			}
		}
		if (disequi)
		{
			//fbj	thisperson& with = *p;
			{
				information& with1 = *p->store_;
				{
					locusvalues& with2 = *thislocus[syste];
					for( i=1; i <= with2.nallele; i ++)
					{
						hap2[syste]=i;
						for( j=i; j <= with2.nallele; j ++)
							if ((! nohom[syste]) || (i!=affall) || (j!=affall))
								if (with1.possible[syste][i][j])
								{
									hap1[syste]=j;
									newval=val;
									getval(syste,i,j,newval, p);
									if (newval!=0.0)
									{
										if (syste!=nlocus) ugetgene(syste+1,newval, p, hap1, hap2);
										else setval(newval, p, hap1, hap2);
									}
								}
					}
				}
			}
		}
	}
	
	
	
	
	static void quanval1(ind& p, int& syste, hapvector& hap1, hapvector& hap2, real& newval, real& val)
	/*Uses this only when thislocus[syste]^.which is not null*/
	{
		int i,j;
		
		
		{
			//fbj	thisperson& with = *p;
			{
				information& with1 = *p->store_;
				{
					locusvalues& with2 = *thislocus[syste];
					for( i=1; i <= with2.nallele; i ++)
					{
						hap1[syste]=i;
						for( j=i; j <= with2.nallele; j ++)
							if ((! nohom[syste]) || (i!=affall) || (j!=affall))
								if (with1.possible[syste][i][j])
								{
									hap2[syste]=j;
									newval=val;
									getval(syste,i,j,newval, p);
									if (newval!=0.0)
									{
										if (syste!=nlocus) ugetgene(syste+1,newval, p, hap1, hap2);
										else setval(newval, p, hap1, hap2);
									}
								}
					}
				}
			}
		}
		if (disequi)
		{
			//fbj	thisperson& with = *p;
			{
				information& with1 = *p->store_;
				{
					locusvalues& with2 = *thislocus[syste];
					for( i=1; i <= with2.nallele; i ++)
					{
						hap2[syste]=i;
						for( j=i; j <= with2.nallele; j ++)
							if ((! nohom[syste]) || (i!=affall) || (j!=affall))
								if (with1.possible[syste][i][j])
								{
									hap1[syste]=j;
									newval=val;
									getval(syste,i,j,newval, p);
									if (newval!=0.0)
									{
										if (syste!=nlocus) ugetgene(syste+1,newval, p, hap1, hap2);
										else setval(newval, p, hap1, hap2);
									}
								}
					}
				}
			}
		}
	}
	
	
	
	static void ugetgene(int syste,
						 real val, ind& p, hapvector& hap1, hapvector& hap2)
	{
		real newval;
		
		
		{
			locusvalues& with = *thislocus[syste];
			if (sexlink && p->male)
				
				switch (with.which) {
					case binary:facmale1(p, syste, hap1, val, hap2); break;
					case affection:affmale1(p, syste, newval, val, hap1, hap2); break;
					case quantitative:quanmale1(p, syste, newval, val, hap1, hap2); break;
					case null:
						if (with.privlocus->which==affection)
							affmale1(p, syste, newval, val, hap1, hap2);
						else quanmale1(p, syste, newval, val, hap1, hap2);
						break;
					case last_locustype: exit_error("thislocus[system]->which = last_locustype (fbj)");
				}
			else
				switch (with.which) {
					case binary:fac1(p, syste, hap1, hap2, val); break;
					case affection:aff1(p, syste, hap1, hap2, newval, val); break;
					case quantitative:quanval1(p, syste, hap1, hap2, newval, val); break;
					case null:aff1(p, syste, hap1, hap2, newval, val); break;
					case last_locustype: exit_error("thislocus[system]->which = last_locustype (fbj)");
				}
		}
	}
	
	
	static void getvect(ind p)
	{
		hapvector hap1,hap2;
		
		if (p->unknown)  ugetgene(1,1.0, p, hap1, hap2);
		else getgene(1,1.0, p, hap1, hap2);
	}
	
	
	static void prob(ind& p, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale)
	
	{
		int i;
		
		/*prob*/
		{
			thisperson& with = *p;
			if (with.gen==nil)
			{
				with.gen = new thisarray;
				{
					thisarray& with1 = *with.gen;
					for( i=1; i <= fgeno; i ++) with1.genarray[i]=0.0;
				}
				if (with.inloop>0)
					if (looppers[thisped][with.inloop][1]==p)
						with.gen->genarray[loopgen[with.inloop]]=holdpoint[with.inloop]->genarray[loopgen[with.inloop]];
					else
						with.gen->genarray[loopgen[with.inloop]]=1.0;
					else
					{
						getvect(p);
						if (p->pa==nil)  nuscale=nuscale+1;
					}
			}
		}
	}       /*prob*/
	
	static void seg(ind& p,ind& q,ind& r,direction peel, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale);
	
	
	static void getapprox(ind& p, int& thisped)
	
	{
		int first;
		real maxval;
		
		/*getapprox*/
		maxval=p->gen->genarray[1];
		{
			thisarray& with = *p->gen;
			for( first=1; first <= fgeno; first ++)
				if (with.genarray[first]>maxval)  maxval=with.genarray[first];}
		{
			approxrec& with = *approxstruct;
			{
				thisarray& with1 = *p->gen;
				for( first=1; first <= fgeno; first ++)
					with.approxarray[thisped][first]=with1.genarray[first]>maxval*epsilon;}}
		if (! lasttime)
		{
			approxrec& with = *approxstruct;
			{
				thisarray& with1 = *p->gen;
				for( first=1; first <= fgeno; first ++)
					if (! with.approxarray[thisped][first])  with1.genarray[first]=0.0;
			}
		}
	}       /*getapprox*/
	
	
	static real msegsex(ind& p, int& fseg, int& secondseg, int& sseg, thetapoint& secondsex, int& sstart, int& send, real& ps, real& pf, int& firstseg, thetapoint& firstsex, int& fstart, int& fend)
	
	{
		int g1,g2,g3,g4,g5,g6,g7,g8,ms,mf,ms1,ms2,mf1,mf2,j,k,f1,f2,s1,s2;
		real val,temp2;
		array1d<1,maxchild,real> temp;
		
		/*msegsex*/
		real msegsex_result;
		for( k=1; k <= nchild; k ++) temp[k]=0.0;
		{
			gennurec& with = *gennustruct;
			if (p->male)
			{
				mf=muthap[fseg];
				secondseg=sseg;
				{
					thetavalues& with1 = *secondsex;
					for( j=sstart; j <= send; j ++)
						if (with1.segprob[j]==0.0)
							secondseg=secondseg+1;
						else
						{
							temp2=with1.segprob[j];
							s1=seghap1[secondseg];
							s2=seghap2[secondseg];
							ms1=muthap[s1];
							ms2=muthap[s2];
							if (s1!=s2)
							{
								g1=with.genenumber[fseg][s1];
								g2=with.genenumber[fseg][s2];
								g3=with.genenumber[fseg][ms1];
								g4=with.genenumber[fseg][ms2];
								g5=with.genenumber[mf][s1];
								g6=with.genenumber[mf][s2];
								g7=with.genenumber[mf][ms1];
								g8=with.genenumber[mf][ms2];
								for( k=1; k <= nchild; k ++)
								{ thisarray& with2 = *thischild[k];
									if (malechild[k])
										temp[k]=temp[k]+temp2*((1-ps)*(with2.genarray[s1]+with2.genarray[s2])
															   +ps*(with2.genarray[ms1]+with2.genarray[ms2]));
									else temp[k]=temp[k]+temp2*((1-pf)*(1-ps)*(with2.genarray[g1]+with2.genarray[g2])
																+(1-pf)*ps*(with2.genarray[g3]+with2.genarray[g4])+pf*(1-ps)*(with2.genarray[g5]+with2.genarray[g6])
																+pf*ps*(with2.genarray[g7]+with2.genarray[g8]));}
							}
							else {
								g1=with.genenumber[fseg][s1];
								g3=with.genenumber[fseg][ms1];
								g5=with.genenumber[mf][s1];
								g7=with.genenumber[mf][ms1];
								for( k=1; k <= nchild; k ++)
								{ thisarray& with2 = *thischild[k];
									if (malechild[k])
										temp[k]=temp[k]+2.0*temp2*((1-ps)*with2.genarray[s1]+ps*with2.genarray[ms1]);
									else temp[k]=temp[k]+2.0*temp2*((1-pf)*(1-ps)*with2.genarray[g1]
																	+ps*(1-pf)*with2.genarray[g3]+pf*(1-ps)*with2.genarray[g5]+pf*ps*with2.genarray[g7]);}
							}
							secondseg=secondseg+1;
						}
				}
			}
			else
			{
				firstseg=fseg;
				ms=muthap[sseg];
				{
					thetavalues& with1 = *firstsex;
					for( j=fstart; j <= fend; j ++)
						if (with1.segprob[j]==0.0)
							firstseg=firstseg+1;
						else
						{
							temp2=with1.segprob[j];
							f1=seghap1[firstseg];
							f2=seghap2[firstseg];
							mf1=muthap[f1];
							mf2=muthap[f2];
							if (f1!=f2)
							{
								g1=with.genenumber[sseg][f1];
								g2=with.genenumber[sseg][f2];
								g3=with.genenumber[sseg][mf1];
								g4=with.genenumber[sseg][mf2];
								g5=with.genenumber[ms][f1];
								g6=with.genenumber[ms][f2];
								g7=with.genenumber[ms][mf1];
								g8=with.genenumber[ms][mf2];
								for( k=1; k <= nchild; k ++)
								{ thisarray& with2 = *thischild[k];
									if (malechild[k])
										temp[k]=temp[k]+temp2*((1-pf)*(with2.genarray[f1]+with2.genarray[f2])
															   +pf*(with2.genarray[mf1]+with2.genarray[mf2]));
									else
										temp[k]=temp[k]+temp2*((1-pf)*(1-ps)*(with2.genarray[g1]+with2.genarray[g2])
															   +(1-ps)*pf*(with2.genarray[g3]+with2.genarray[g4])+ps*(1-pf)*(with2.genarray[g5]+with2.genarray[g6])
															   +pf*ps*(with2.genarray[g7]+with2.genarray[g8]));}
							}
							else {
								g1=with.genenumber[sseg][f1];
								g3=with.genenumber[sseg][mf1];
								g5=with.genenumber[ms][f1];
								g7=with.genenumber[ms][mf1];
								for( k=1; k <= nchild; k ++)
								{ thisarray& with2 = *thischild[k];
									if (malechild[k])
										temp[k]=temp[k]+2.0*temp2*((1-pf)*with2.genarray[f1]+pf*with2.genarray[mf1]);
									else
										temp[k]=temp[k]+2.0*temp2*((1-pf)*(1-ps)*with2.genarray[g1]
																   +pf*(1-ps)*with2.genarray[g3]+ps*(1-pf)*with2.genarray[g5]+ps*pf*with2.genarray[g7]);}
							}
							firstseg=firstseg+1;
						}
				}
			}
		}
		val=1.0;
		for( k=1; k <= nchild; k ++) val=val*temp[k];
		msegsex_result=val;
		return msegsex_result;
	}       /*msegsex*/
	
	
	static real msegsexf(int& fseg, int& send, int& sstart, int& secondseg, int& sseg, thetapoint& secondsex, real& ps, real& pf)
	
	{
		int g1,g2,g3,g4,g5,g6,g7,g8,mf,ms1,ms2,j,k,l,s1,s2,
		slength;
		real val,temp1;
		array2d<1,maxchild,1,maxseg,real> temp;
		array1d<1,maxseg,real> temp2;
		
		/*msegsexf*/
		real msegsexf_result;
		mf=muthap[fseg];
		slength=send-sstart+1;
		for( k=1; k <= nchild; k ++)
			for( l=1; l <= slength; l ++)
				temp[k][l]=0.0;
		secondseg=sseg;
		{
			gennurec& with = *gennustruct;
			{
				thetavalues& with1 = *secondsex;
				for( j=sstart; j <= send; j ++)
				{
					for( l=1; l <= slength; l ++)
						temp2[l]=with1.segprob[j+(l-1)*slength];
					s1=seghap1[secondseg];
					s2=seghap2[secondseg];
					ms1=muthap[s1];
					ms2=muthap[s2];
					if (s1!=s2)
					{
						g1=with.genenumber[fseg][s1];
						g2=with.genenumber[fseg][s2];
						g3=with.genenumber[fseg][ms1];
						g4=with.genenumber[fseg][ms2];
						g5=with.genenumber[mf][s1];
						g6=with.genenumber[mf][s2];
						g7=with.genenumber[mf][ms1];
						g8=with.genenumber[mf][ms2];
						for( k=1; k <= nchild; k ++)
						{
							{ thisarray& with2 = *thischild[k];
								if (malechild[k])
									val=(1-ps)*(with2.genarray[s1]+with2.genarray[s2])+ps*(with2.genarray[ms1]+with2.genarray[ms2]);
								else val=(1-pf)*(1-ps)*(with2.genarray[g1]+with2.genarray[g2])
									+(1-pf)*ps*(with2.genarray[g3]+with2.genarray[g4])+pf*(1-ps)*(with2.genarray[g5]+with2.genarray[g6])
									+pf*ps*(with2.genarray[g7]+with2.genarray[g8]);}
							for( l=1; l <= slength; l ++)
								temp[k][l]=temp[k][l]+temp2[l]*val;
						}
					}
					else {
						g1=with.genenumber[fseg][s1];
						g3=with.genenumber[fseg][ms1];
						g5=with.genenumber[mf][s1];
						g7=with.genenumber[mf][ms1];
						for( k=1; k <= nchild; k ++)
						{
							{ thisarray& with2 = *thischild[k];
								if (malechild[k])
									val=2.0*((1-ps)*with2.genarray[s1]+ps*with2.genarray[ms1]);
								else val=2.0*((1-pf)*(1-ps)*with2.genarray[g1]+ps*(1-pf)*with2.genarray[g3]
											  +pf*(1-ps)*with2.genarray[g5]+ps*pf*with2.genarray[g7]);}
							for( l=1; l <= slength; l ++)
								temp[k][l]=temp[k][l]+temp2[l]*val;
						}
					}
					secondseg=secondseg+1;
				}
			}
		}
		temp1=0.0;
		for( l=1; l <= slength; l ++)
		{
			val=1.0;
			for( k=1; k <= nchild; k ++) val=val*temp[k][l];
			temp1=temp1+val;
		}
		msegsexf_result=temp1;
		return msegsexf_result;
	}       /*msegsexf*/
	
	
	
	
	static real segsex(ind& p, int& secondseg, int& sseg, thetapoint& secondsex, int& sstart, int& send, int& fseg, int& firstseg, thetapoint& firstsex, int& fstart, int& fend)
	
	{
		int g1,g2,j,k,f1,f2,s1,s2;
		real val,temp2;
		array1d<1,maxchild,real> temp;
		
		/*segsex*/
		real segsex_result;
		for( k=1; k <= nchild; k ++) temp[k]=0.0;
		{
			gennurec& with = *gennustruct;
			if (p->male)
			{
				secondseg=sseg;
				{
					thetavalues& with1 = *secondsex;
					for( j=sstart; j <= send; j ++)
						if (with1.segprob[j]==0.0)
							secondseg=secondseg+1;
						else
						{
							temp2=with1.segprob[j];
							s1=seghap1[secondseg];
							s2=seghap2[secondseg];
							if (s1!=s2)
							{
								g1=with.genenumber[fseg][s1];
								g2=with.genenumber[fseg][s2];
								for( k=1; k <= nchild; k ++)
								{ thisarray& with2 = *thischild[k];
									if (malechild[k])
										temp[k]=temp[k]+temp2*(with2.genarray[s1]+with2.genarray[s2]);
									else temp[k]=temp[k]+temp2*(with2.genarray[g1]+with2.genarray[g2]);}
							}
							else {
								g1=with.genenumber[fseg][s1];
								for( k=1; k <= nchild; k ++)
								{ thisarray& with2 = *thischild[k];
									if (malechild[k])
										temp[k]=temp[k]+2.0*temp2*with2.genarray[s1];
									else temp[k]=temp[k]+2.0*temp2*with2.genarray[g1];}
							}
							secondseg=secondseg+1;
						}
				}
			}
			else {
				firstseg=fseg;
				{
					thetavalues& with1 = *firstsex;
					for( j=fstart; j <= fend; j ++)
						if (with1.segprob[j]==0.0)
							firstseg=firstseg+1;
						else
						{
							temp2=with1.segprob[j];
							f1=seghap1[firstseg];
							f2=seghap2[firstseg];
							if (f1!=f2)
							{
								g1=with.genenumber[sseg][f1];
								g2=with.genenumber[sseg][f2];
								for( k=1; k <= nchild; k ++)
								{ thisarray& with2 = *thischild[k];
									if (malechild[k])
										temp[k]=temp[k]+temp2*(with2.genarray[f1]+with2.genarray[f2]);
									else temp[k]=temp[k]+temp2*(with2.genarray[g1]+with2.genarray[g2]);}
							}
							else {
								g1=with.genenumber[sseg][f1];
								for( k=1; k <= nchild; k ++)
								{ thisarray& with2 = *thischild[k];
									if (malechild[k])
										temp[k]=temp[k]+2.0*temp2*with2.genarray[f1];
									else temp[k]=temp[k]+2.0*temp2*with2.genarray[g1];}
							}
							firstseg=firstseg+1;
						}
				}
			}
		}
		val=1.0;
		for( k=1; k <= nchild; k ++) val=val*temp[k];
		segsex_result=val;
		return segsex_result;
	}       /*segsex*/
	
	
	
	static real segsexf(int& send, int& sstart, int& secondseg, int& sseg, thetapoint& secondsex, int& fseg)
	{
		int g1,g2,j,k,l,s1,s2,slength;
		real val,temp1;
		array2d<1,maxchild,1,maxseg,real> temp;
		array1d<1,maxseg,real> temp2;
		
		
		real segsexf_result;
		slength=send-sstart+1;
		for( k=1; k <= nchild; k ++)
			for( l=1; l <= slength; l ++) temp[k][l]=0.0;
		secondseg=sseg;
		{
			gennurec& with = *gennustruct;
			{
				thetavalues& with1 = *secondsex;
				for( j=sstart; j <= send; j ++)
				{
					for( l=1; l <= slength; l ++)
						temp2[l]=with1.segprob[j+(l-1)*slength];
					s1=seghap1[secondseg];
					s2=seghap2[secondseg];
					if (s1!=s2)
					{
						g1=with.genenumber[fseg][s1];
						g2=with.genenumber[fseg][s2];
						for( k=1; k <= nchild; k ++)
						{
							thisarray& with2 = *thischild[k];
							
							if (malechild[k])
								val=with2.genarray[s1]+with2.genarray[s2];
							else
								val=with2.genarray[g1]+with2.genarray[g2];
							if (val!=0.0)
								for( l=1; l <= slength; l ++)
									temp[k][l]=temp[k][l]+temp2[l]*val;
						}
					}
					else
					{
						g1=with.genenumber[fseg][s1];
						for( k=1; k <= nchild; k ++)
						{
							thisarray& with2 = *thischild[k];
							
							if (malechild[k])
								val=2.0*with2.genarray[s1];
							else
								val=2.0*with2.genarray[g1];
							if (val!=0.0)
								for( l=1; l <= slength; l ++)
									temp[k][l]=temp[k][l]+temp2[l]*val;
						}
					}
					secondseg=secondseg+1;
				}
			}
		}
		temp1=0.0;
		for( l=1; l <= slength; l ++)
		{
			val=1.0;
			for( k=1; k <= nchild; k ++) val=val*temp[k][l];
			temp1=temp1+val;
		}
		segsexf_result=temp1;
		return segsexf_result;
	}
	
	
	
	static real segfun(int& firstseg, int& fseg, thetapoint& firstsex, int& fstart, int& fend, int& secondseg, int& sseg, thetapoint& secondsex, int& sstart, int& send)
	{
		int g1,g2,g3,g4,i,j,k,f1,f2,s1,s2;
		real val,temp1,temp2;
		array1d<1,maxchild,real> temp;
		
		
		real segfun_result;
		firstseg=fseg;
		for( k=1; k <= nchild; k ++) temp[k]=0.0;
		{
			gennurec& with = *gennustruct;
			{
				thetavalues& with1 = *firstsex;
				for( i=fstart; i <= fend; i ++)
					if (with1.segprob[i]==0.0)
						firstseg=firstseg+1;
					else
					{
						temp1=with1.segprob[i];
						f1=seghap1[firstseg];
						f2=seghap2[firstseg];
						secondseg=sseg;
						{
							thetavalues& with2 = *secondsex;
							if (f1!=f2)
							{
								for( j=sstart; j <= send; j ++)
									if (with2.segprob[j]==0.0)
										secondseg=secondseg+1;
									else
									{
										temp2=temp1*with2.segprob[j];
										s1=seghap1[secondseg];
										s2=seghap2[secondseg];
										if (s1!=s2)
										{
											g1=with.genenumber[f1][s1];
											g2=with.genenumber[f1][s2];
											g3=with.genenumber[f2][s1];
											g4=with.genenumber[f2][s2];
											for( k=1; k <= nchild; k ++)
											{ thisarray& with3 = *thischild[k];
												temp[k]=temp[k]+temp2*(
																	   with3.genarray[g1]+
																	   with3.genarray[g2]+
																	   with3.genarray[g3]+
																	   with3.genarray[g4]);}
										}
										else
										{
											g1=with.genenumber[f1][s1];
											g3=with.genenumber[f2][s1];
											for( k=1; k <= nchild; k ++)
											{ thisarray& with3 = *thischild[k];
												temp[k]=temp[k]+2*temp2*(
																		 with3.genarray[g1]+
																		 with3.genarray[g3]);}
										}
										secondseg=secondseg+1;
									}
							}
							else
								for( j=sstart; j <= send; j ++)
									if (with2.segprob[j]==0.0)
										secondseg=secondseg+1;
									else
									{
										temp2=temp1*with2.segprob[j];
										s1=seghap1[secondseg];
										s2=seghap2[secondseg];
										if (s1!=s2)
										{
											g1=with.genenumber[f1][s1];
											g2=with.genenumber[f1][s2];
											for( k=1; k <= nchild; k ++)
											{ thisarray& with3 = *thischild[k];
												temp[k]=temp[k]+2*temp2*(
																		 with3.genarray[g1]+
																		 with3.genarray[g2]);}
										}
										else
										{
											g1=with.genenumber[f1][s1];
											for( k=1; k <= nchild; k ++)
											{ thisarray& with3 = *thischild[k];
												temp[k]=temp[k]+4*temp2*(
																		 with3.genarray[g1]);}
										}
										secondseg=secondseg+1;
									}
						}
						firstseg=firstseg+1;
					}
			}
		}
		val=1.0;
		for( k=1; k <= nchild; k ++) val=val*temp[k];
		segfun_result=val;
		return segfun_result;
	}
	
	
	
	
	static real msegfast(int& firstseg, int& fseg, int& send, int& sstart, thetapoint& firstsex, int& fstart, int& fend, int& secondseg, int& sseg, thetapoint& secondsex, real& ps, real& pf)
	
	{
		int g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,
		i,j,k,l,f1,f2,s1,s2,ms1,ms2,mf1,mf2,slength;
		real val,temp1;
		array2d<1,maxchild,1,maxseg,real> temp;
		array1d<1,maxseg,real> temp2;
		
		/*msegfast*/
		real msegfast_result;
		firstseg=fseg;
		slength=send-sstart+1;
		for( k=1; k <= nchild; k ++)
			for( l=1; l <= slength; l ++)
				temp[k][l]=0.0;
		{
			gennurec& with = *gennustruct;
			{
				thetavalues& with1 = *firstsex;
				for( i=fstart; i <= fend; i ++)
					if (with1.segprob[i]==0.0)
						firstseg=firstseg+1;
					else
					{
						temp1=with1.segprob[i];
						f1=seghap1[firstseg];
						f2=seghap2[firstseg];
						mf1=muthap[f1];
						mf2=muthap[f2];
						secondseg=sseg;
						{ thetavalues& with2 = *secondsex;
							if (f1!=f2)
								for( j=sstart; j <= send; j ++)
								{
									for( l=1; l <= slength; l ++)
										temp2[l]=temp1*with2.segprob[j+(l-1)*slength];
									s1=seghap1[secondseg];
									s2=seghap2[secondseg];
									ms1=muthap[s1];
									ms2=muthap[s2];
									if (s1!=s2)
									{
										g1=with.genenumber[f1][s1];
										g2=with.genenumber[f1][s2];
										g3=with.genenumber[f2][s1];
										g4=with.genenumber[f2][s2];
										g5=with.genenumber[f1][ms1];
										g6=with.genenumber[f1][ms2];
										g7=with.genenumber[f2][ms1];
										g8=with.genenumber[f2][ms2];
										g9=with.genenumber[mf1][s1];
										g10=with.genenumber[mf1][s2];
										g11=with.genenumber[mf2][s1];
										g12=with.genenumber[mf2][s2];
										g13=with.genenumber[mf1][ms1];
										g14=with.genenumber[mf1][ms2];
										g15=with.genenumber[mf2][ms1];
										g16=with.genenumber[mf2][ms2];
										for( k=1; k <= nchild; k ++)
										{
											thisarray& with3 = *thischild[k];
											
											val=(1-ps)*(1-pf)*(with3.genarray[g1]+with3.genarray[g2]
															   +with3.genarray[g3]+with3.genarray[g4])+ps*(1-pf)*(with3.genarray[g5]+with3.genarray[g6]
																												  +with3.genarray[g7]+with3.genarray[g8])+pf*(1-ps)*(with3.genarray[g9]+with3.genarray[g10]
																																									 +with3.genarray[g11]+with3.genarray[g12])+pf*ps*(with3.genarray[g13]+with3.genarray[g14]
																																																					  +with3.genarray[g15]+with3.genarray[g16]);
											if (val!=0.0)
												for( l=1; l <= slength; l ++)
													temp[k][l]=temp[k][l]+temp2[l]*val;
										}
									}
									else {
										g1=with.genenumber[f1][s1];
										g3=with.genenumber[f2][s1];
										g5=with.genenumber[f1][ms1];
										g7=with.genenumber[f2][ms1];
										g9=with.genenumber[mf1][s1];
										g11=with.genenumber[mf2][s1];
										g13=with.genenumber[mf1][ms1];
										g15=with.genenumber[mf2][ms1];
										for( k=1; k <= nchild; k ++)
										{
											thisarray& with3 = *thischild[k];
											
											val=2*((1-ps)*(1-pf)*(with3.genarray[g1]+with3.genarray[g3])
												   +ps*(1-pf)*(with3.genarray[g5]+with3.genarray[g7])+pf*(1-ps)*(with3.genarray[g9]
																												 +with3.genarray[g11])+pf*ps*(with3.genarray[g13]+with3.genarray[g15]));
											for( l=1; l <= slength; l ++)
												temp[k][l]=temp[k][l]+temp2[l]*val;
										}
									}
									secondseg=secondseg+1;
								}
							else for( j=sstart; j <= send; j ++)
							{
								for( l=1; l <= slength; l ++)
									temp2[l]=temp1*with2.segprob[j+(l-1)*slength];
								s1=seghap1[secondseg];
								s2=seghap2[secondseg];
								ms1=muthap[s1];
								ms2=muthap[s2];
								if (s1!=s2)
								{
									g1=with.genenumber[f1][s1];
									g2=with.genenumber[f1][s2];
									g5=with.genenumber[f1][ms1];
									g6=with.genenumber[f1][ms2];
									g9=with.genenumber[mf1][s1];
									g10=with.genenumber[mf1][s2];
									g13=with.genenumber[mf1][ms1];
									g14=with.genenumber[mf1][ms2];
									for( k=1; k <= nchild; k ++)
									{
										thisarray& with3 = *thischild[k];
										
										val=2*((1-ps)*(1-pf)*(with3.genarray[g1]+with3.genarray[g2])
											   +ps*(1-pf)*(with3.genarray[g5]+with3.genarray[g6])+pf*(1-ps)*(with3.genarray[g9]
																											 +with3.genarray[g10])+pf*ps*(with3.genarray[g13]+with3.genarray[g14]));
										if (val!=0.0)
											for( l=1; l <= slength; l ++)
												temp[k][l]=temp[k][l]+temp2[l]*val;
									}
								}
								else {
									g1=with.genenumber[f1][s1];
									g5=with.genenumber[f1][ms1];
									g9=with.genenumber[mf1][s1];
									g13=with.genenumber[mf1][ms1];
									for( k=1; k <= nchild; k ++)
									{
										thisarray& with3 = *thischild[k];
										
										val=4*((1-ps)*(1-pf)*with3.genarray[g1]+ps*(1-pf)*with3.genarray[g5]
											   +pf*(1-ps)*with3.genarray[g9]+pf*ps*with3.genarray[g13]);
										if (val!=0.0)
											for( l=1; l <= slength; l ++)
												temp[k][l]=temp[k][l]+temp2[l]*val;
									}
								}
								secondseg=secondseg+1;
							}
						}
						firstseg=firstseg+1;
					}
			}
		}
		temp1=0.0;
		for( l=1; l <= slength; l ++)
		{
			val=1.0;
			for( k=1; k <= nchild; k ++) val=val*temp[k][l];
			temp1=temp1+val;
		}
		msegfast_result=temp1;
		return msegfast_result;
	}       /*msegfast*/
	
	
	static real msegfun(int& firstseg, int& fseg, thetapoint& firstsex, int& fstart, int& fend, int& secondseg, int& sseg, thetapoint& secondsex, int& sstart, int& send, real& ps, real& pf)
	{
		int g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15,g16,i,j,k,f1,f2,s1,s2,ms1,ms2,mf1,mf2;
		real val,temp1,temp2;
		array1d<1,maxchild,real> temp;
		
		/*msegfun*/
		real msegfun_result;
		firstseg=fseg;
		for( k=1; k <= nchild; k ++) temp[k]=0.0;
		{
			gennurec& with = *gennustruct;
			{
				thetavalues& with1 = *firstsex;
				for( i=fstart; i <= fend; i ++)
					if (with1.segprob[i]==0.0)
						firstseg=firstseg+1;
					else
					{
						temp1=with1.segprob[i];
						f1=seghap1[firstseg];
						f2=seghap2[firstseg];
						mf1=muthap[f1];
						mf2=muthap[f2];
						secondseg=sseg;
						{
							thetavalues& with2 = *secondsex;
							if (f1!=f2)
							{
								for( j=sstart; j <= send; j ++)
									if (with2.segprob[j]==0.0)
										secondseg=secondseg+1;
									else
									{
										temp2=temp1*with2.segprob[j];
										s1=seghap1[secondseg];
										s2=seghap2[secondseg];
										ms1=muthap[s1];
										ms2=muthap[s2];
										if (s1!=s2)
										{
											g1=with.genenumber[f1][s1];
											g2=with.genenumber[f1][s2];
											g3=with.genenumber[f2][s1];
											g4=with.genenumber[f2][s2];
											g5=with.genenumber[f1][ms1];
											g6=with.genenumber[f1][ms2];
											g7=with.genenumber[f2][ms1];
											g8=with.genenumber[f2][ms2];
											g9=with.genenumber[mf1][s1];
											g10=with.genenumber[mf1][s2];
											g11=with.genenumber[mf2][s1];
											g12=with.genenumber[mf2][s2];
											g13=with.genenumber[mf1][ms1];
											g14=with.genenumber[mf1][ms2];
											g15=with.genenumber[mf2][ms1];
											g16=with.genenumber[mf2][ms2];
											for( k=1; k <= nchild; k ++)
											{ thisarray& with3 = *thischild[k];
												temp[k]=temp[k]+temp2*((1-ps)*(1-pf)*(
																					  with3.genarray[g1]+with3.genarray[g2]+
																					  with3.genarray[g3]+with3.genarray[g4])+ps*(1-pf)*(
																																		with3.genarray[g5]+with3.genarray[g6]+
																																		with3.genarray[g7]+with3.genarray[g8])+pf*(1-ps)*(
																																														  with3.genarray[g9]+with3.genarray[g10]+
																																														  with3.genarray[g11]+with3.genarray[g12])+pf*ps*(
																																																										  with3.genarray[g13]+with3.genarray[g14]+
																																																										  with3.genarray[g15]+with3.genarray[g16]));}
										}
										else
										{
											g1=with.genenumber[f1][s1];
											g3=with.genenumber[f2][s1];
											g5=with.genenumber[f1][ms1];
											g7=with.genenumber[f2][ms1];
											g9=with.genenumber[mf1][s1];
											g11=with.genenumber[mf2][s1];
											g13=with.genenumber[mf1][ms1];
											g15=with.genenumber[mf2][ms1];
											for( k=1; k <= nchild; k ++)
											{ thisarray& with3 = *thischild[k];
												temp[k]=temp[k]+2*temp2*((1-ps)*(1-pf)*(
																						with3.genarray[g1]+with3.genarray[g3])+ps*(1-pf)*(
																																		  with3.genarray[g5]+with3.genarray[g7])+pf*(1-ps)*(
																																															with3.genarray[g9]+with3.genarray[g11])+pf*ps*(
																																																										   with3.genarray[g13]+with3.genarray[g15]));}
										}
										secondseg=secondseg+1;
									}
							}
							else
								for( j=sstart; j <= send; j ++)
									if (with2.segprob[j]==0.0)
										secondseg=secondseg+1;
									else
									{
										temp2=temp1*with2.segprob[j];
										s1=seghap1[secondseg];
										s2=seghap2[secondseg];
										ms1=muthap[s1];
										ms2=muthap[s2];
										if (s1!=s2)
										{
											g1=with.genenumber[f1][s1];
											g2=with.genenumber[f1][s2];
											g5=with.genenumber[f1][ms1];
											g6=with.genenumber[f1][ms2];
											g9=with.genenumber[mf1][s1];
											g10=with.genenumber[mf1][s2];
											g13=with.genenumber[mf1][ms1];
											g14=with.genenumber[mf1][ms2];
											for( k=1; k <= nchild; k ++)
											{	thisarray& with3 = *thischild[k];
												temp[k]=temp[k]+2*temp2*((1-ps)*(1-pf)*(
																						with3.genarray[g1]+with3.genarray[g2])+ps*(1-pf)*(
																																		  with3.genarray[g5]+with3.genarray[g6])+pf*(1-ps)*(
																																															with3.genarray[g9]+with3.genarray[g10])+pf*ps*(
																																																										   with3.genarray[g13]+with3.genarray[g14]));}
										}
										else
										{
											g1=with.genenumber[f1][s1];
											g5=with.genenumber[f1][ms1];
											g9=with.genenumber[mf1][s1];
											g13=with.genenumber[mf1][ms1];
											for( k=1; k <= nchild; k ++)
											{	thisarray& with3 = *thischild[k];
												temp[k]=temp[k]+4*temp2*(
																		 (1-ps)*(1-pf)*with3.genarray[g1]+
																		 ps*(1-pf)*with3.genarray[g5]+pf*(1-ps)*with3.genarray[g9]
																		 +pf*ps*with3.genarray[g13]);}
										}
										secondseg=secondseg+1;
									}
						}
						firstseg=firstseg+1;
					}
			}
		}
		val=1.0;
		for( k=1; k <= nchild; k ++) val=val*temp[k];
		msegfun_result=val;
		return msegfun_result;
	}       /*msegfun*/
	
	
	static real segfast(int& firstseg, int& fseg, int& send, int& sstart, thetapoint& firstsex, int& fstart, int& fend, int& secondseg, int& sseg, thetapoint& secondsex)
	{
		int g1,g2,g3,g4,i,j,k,l,f1,f2,s1,s2,slength;
		real val,temp1;
		array2d<1,maxchild,1,maxseg,real> temp;
		array1d<1,maxseg,real> temp2;
		
		/*segfast*/
		real segfast_result;
		firstseg=fseg;
		slength=send-sstart+1;
		for( k=1; k <= nchild; k ++)
			for( l=1; l <= slength; l ++) temp[k][l]=0.0;
		{
			gennurec& with = *gennustruct;
			{
				thetavalues& with1 = *firstsex;
				for( i=fstart; i <= fend; i ++)
					if (with1.segprob[i]==0.0)
						firstseg=firstseg+1;
					else
					{
						temp1=with1.segprob[i];
						f1=seghap1[firstseg];
						f2=seghap2[firstseg];
						secondseg=sseg;
						{
							thetavalues& with2 = *secondsex;
							if (f1!=f2)
								for( j=sstart; j <= send; j ++)
								{
									for( l=1; l <= slength; l ++)
										temp2[l]=temp1*with2.segprob[j+(l-1)*slength];
									s1=seghap1[secondseg];
									s2=seghap2[secondseg];
									if (s1!=s2)
									{
										g1=with.genenumber[f1][s1];
										g2=with.genenumber[f1][s2];
										g3=with.genenumber[f2][s1];
										g4=with.genenumber[f2][s2];
										for( k=1; k <= nchild; k ++)
										{
											thisarray& with3 = *thischild[k];
											
											val=with3.genarray[g1]+with3.genarray[g2]+with3.genarray[g3]+with3.genarray[g4];
											if (val!=0.0)
												for( l=1; l <= slength; l ++)
													temp[k][l]=temp[k][l]+temp2[l]*val;
										}
									}
									else {
										g1=with.genenumber[f1][s1];
										g3=with.genenumber[f2][s1];
										for( k=1; k <= nchild; k ++)
										{
											thisarray& with3 = *thischild[k];
											
											val=2.0*(with3.genarray[g1]+with3.genarray[g3]);
											if (val!=0.0)
												for( l=1; l <= slength; l ++)
													temp[k][l]=temp[k][l]+temp2[l]*val;
										}
									}
									secondseg=secondseg+1;
								}
							else for( j=sstart; j <= send; j ++)
							{
								for( l=1; l <= slength; l ++) temp2[l]=temp1*with2.segprob[j+(l-1)*slength];
								s1=seghap1[secondseg];
								s2=seghap2[secondseg];
								if (s1!=s2)
								{
									g1=with.genenumber[f1][s1];
									g2=with.genenumber[f1][s2];
									for( k=1; k <= nchild; k ++)
									{
										thisarray& with3 = *thischild[k];
										
										val=2.0*(with3.genarray[g1]+with3.genarray[g2]);
										if (val!=0.0)
											for( l=1; l <= slength; l ++)
												temp[k][l]=temp[k][l]+temp2[l]*val;
									}
								}
								else {
									g1=with.genenumber[f1][s1];
									for( k=1; k <= nchild; k ++)
									{
										thisarray& with3 = *thischild[k];
										
										val=4.0*with3.genarray[g1];
										if (val!=0.0)
											for( l=1; l <= slength; l ++)
												temp[k][l]=temp[k][l]+temp2[l]*val;
									}
								}
								secondseg=secondseg+1;
							}
						}
						firstseg=firstseg+1;
					}
			}
		}
		temp1=0.0;
		for( l=1; l <= slength; l ++)
		{
			val=1.0;
			for( k=1; k <= nchild; k ++) val=val*temp[k][l];
			temp1=temp1+val;
		}
		segfast_result=temp1;
		return segfast_result;
	}       /*segfast*/
	
	
	static void initseg(ind& p, int& nfirst, int& nsecond, thetapoint& firstsex, thetapoint& secondsex, real& pf, real& ps, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale, ind& father, ind& mother, ind& child)
	{     /*initseg*/
		if (p->male)
		{
			nfirst=mgeno;
			nsecond=fgeno;
			firstsex=maletheta;
			secondsex=femaletheta;
			pf=mutmale;
			ps=mutfemale;
		}
		else {
			nfirst=fgeno;
			nsecond=mgeno;
			firstsex=femaletheta;
			secondsex=maletheta;
			pf=mutfemale;
			ps=mutmale;
		}
		prob(father, thisped, loopgen, holdpoint, nuscale);
		prob(mother, thisped, loopgen, holdpoint, nuscale);
		child=father->foff;
		nchild=0;
		do {
			prob(child, thisped, loopgen, holdpoint, nuscale);
			if ((child->ma==mother) && (! child->up))
			{
				nchild=nchild+1;
				if (nchild > maxchild) exit_error("There are too many sibs in one of the pedigrees.");
				thischild[nchild]=child->gen;
				malechild[nchild]=child->male;
			}
			child=child->nextpa;
		} while (!(child==nil));
	}       /*initseg*/
	
	
	static void exitseg(int& nuscale, ind& child, ind& father, ind& mother)
	
	{    /*exitseg*/
		nuscale=nuscale+1;
		child=father->foff;
		do {
			if ((child->ma==mother) && ! child->up)
				cleanup(child);
			child=child->nextpa;
		} while (!(child==nil));
	}       /*exitseg*/
	
	
	static void segsexctop(ind& p, int& nfirst, int& nsecond, thetapoint& firstsex, thetapoint& secondsex, real& pf, real& ps, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale, ind& father, ind& mother, ind& child, int& fseg, ind& q, int& sstart, int& send, int& sseg, int& secondseg, int& fstart, int& fend, int& firstseg)
	{
		real segval;
		int first,second,sstop;
		boolean thiscensor,thisrare;
		
		/*segsexctop*/
		initseg(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child);
		{
			censorrec& with = *censorstruct;
			if (p->male)
			{
				{
					thisarray& with1 = *p->gen;
					for( first=1; first <= nfirst; first ++)
						if (with1.genarray[first]!=0.0)
						{
							segval=0.0;
							fseg=first;
							thisrare=rare[first];
							second=1;
							if (thisrare)
							{
								thisarray& with2 = *q->gen;
								do {
									sstart=probstart[second];
									send=probend[second];
									sseg=segstart[second];
									sstop=second+(send-sstart)+1;
									if (thisc<maxcensor)
									{
										thisc=thisc+1;
										thiscensor=with.censor[thisc];
									}
									else thiscensor=(with2.genarray[second]==0.0)
										|| rare[second];
									if (! thiscensor)
									{	if (mutsys!=0)
										segval=segval+with2.genarray[second]*msegsexf(fseg, send, sstart, secondseg, sseg, secondsex, ps, pf);
									else segval=segval+with2.genarray[second]*segsexf(send, sstart, secondseg, sseg, secondsex, fseg);
									}
									second=sstop;
								} while (!(second>nsecond));}
							else
							{
								thisarray& with2 = *q->gen;
								do {
									sstart=probstart[second];
									send=probend[second];
									sseg=segstart[second];
									sstop=second+(send-sstart)+1;
									if (thisc<maxcensor)
									{
										thisc=thisc+1;
										thiscensor=with.censor[thisc];
									}
									else thiscensor=(with2.genarray[second]==0.0);
									if (! thiscensor)
									{	if (mutsys!=0)
										segval=segval+with2.genarray[second]*msegsexf(fseg, send, sstart, secondseg, sseg, secondsex, ps, pf);
									else segval=segval+with2.genarray[second]*segsexf(send, sstart, secondseg, sseg, secondsex, fseg);
									}
									second=sstop;
								} while (!(second>nsecond));}
							with1.genarray[first]=with1.genarray[first]*segval*segscale;
						}
				}
			}
			else { thisarray& with1 = *p->gen;
				for( first=1; first <= nfirst; first ++)
					if (with1.genarray[first]!=0.0)
					{
						segval=0.0;
						fstart=probstart[first];
						fend=probend[first];
						fseg=segstart[first];
						thisrare=rare[first];
						second=1;
						if (thisrare)
						{
							thisarray& with2 = *q->gen;
							do {
								sseg=second;
								if (thisc<maxcensor)
								{
									thisc=thisc+1;
									thiscensor=with.censor[thisc];
								}
								else thiscensor=(with2.genarray[second]==0.0)
									|| rare[second];
								if (! thiscensor)
								{	if (mutsys!=0)
									segval=segval+with2.genarray[second]*msegsex(p, fseg, secondseg, sseg, secondsex, sstart, send, ps, pf, firstseg, firstsex, fstart, fend);
								else segval=segval+with2.genarray[second]*segsex(p, secondseg, sseg, secondsex, sstart, send, fseg, firstseg, firstsex, fstart, fend);
								}
								second=second+1;
							} while (!(second>nsecond));}
						else
						{
							thisarray& with2 = *q->gen;
							do {
								sseg=second;
								if (thisc<maxcensor)
								{
									thisc=thisc+1;
									thiscensor=with.censor[thisc];
								}
								else thiscensor=(with2.genarray[second]==0.0);
								if (! thiscensor)
								{	if (mutsys!=0)
									segval=segval+with2.genarray[second]*msegsex(p, fseg, secondseg, sseg, secondsex, sstart, send, ps, pf, firstseg, firstsex, fstart, fend);
								else segval=segval+with2.genarray[second]*segsex(p, secondseg, sseg, secondsex, sstart, send, fseg, firstseg, firstsex, fstart, fend);
								}
								second=second+1;
							} while (!(second>nsecond));}
						with1.genarray[first]=with1.genarray[first]*segval*segscale;
					}
			}
		}
		cleanup(q);
		exitseg(nuscale, child, father, mother);
	}       /*segsexctop*/
	
	
	static void segsextop(ind& p, int& nfirst, int& nsecond, thetapoint& firstsex, thetapoint& secondsex, real& pf, real& ps, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale, ind& father, ind& mother, ind& child, int& fseg, ind& q, int& sstart, int& send, int& sseg, int& secondseg, int& fstart, int& fend, int& firstseg)
	
	{
		real segval,val;
		int first,second,sstop;
		
		/*segsextop*/
		initseg(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child);
		{
			censorrec& with = *censorstruct;
			if (p->male)
			{
				{
					thisarray& with1 = *p->gen;
					for( first=1; first <= nfirst; first ++)
						if (with1.genarray[first]!=0.0)
						{
							segval=0.0;
							fseg=first;
							second=1;
							{ thisarray& with2 = *q->gen;
								do {
									sstart=probstart[second];
									send=probend[second];
									sseg=segstart[second];
									sstop=second+(send-sstart)+1;
									if (thisc<=maxcensor)  thisc=thisc+1;
									if (with2.genarray[second]==0.0)
									{
										second=sstop;
										if (thisc<=maxcensor)  with.censor[thisc]=true;
									}
									else {
										if (mutsys!=0)
											val=msegsexf(fseg, send, sstart, secondseg, sseg, secondsex, ps, pf);
										else val=segsexf(send, sstart, secondseg, sseg, secondsex, fseg);
										if (thisc<=maxcensor)  with.censor[thisc]=(val==0.0);
										segval=segval+with2.genarray[second]*val;
										second=sstop;
									}
								} while (!(second>nsecond));}
							with1.genarray[first]=with1.genarray[first]*segval*segscale;
						}
				}
			}
			else { thisarray& with1 = *p->gen;
				for( first=1; first <= nfirst; first ++)
					if (with1.genarray[first]!=0.0)
					{
						segval=0.0;
						fstart=probstart[first];
						fend=probend[first];
						fseg=segstart[first];
						second=1;
						{ thisarray& with2 = *q->gen;
							do {
								sseg=second;
								if (thisc<=maxcensor)  thisc=thisc+1;
								if (with2.genarray[second]==0.0)
								{
									second=second+1;
									if (thisc<=maxcensor)  with.censor[thisc]=true;
								}
								else {
									if (mutsys!=0)
										val=msegsex(p, fseg, secondseg, sseg, secondsex, sstart, send, ps, pf, firstseg, firstsex, fstart, fend);
									else val=segsex(p, secondseg, sseg, secondsex, sstart, send, fseg, firstseg, firstsex, fstart, fend);
									if (thisc<=maxcensor)  with.censor[thisc]=(val==0.0);
									segval=segval+with2.genarray[second]*val;
									second=second+1;
								}
							} while (!(second>nsecond));}
						with1.genarray[first]=with1.genarray[first]*segval*segscale;
					}
			}
		}
		cleanup(q);
		exitseg(nuscale, child, father, mother);
	}       /*segsextop*/
	
	
	static void segsexup(ind& p, int& nfirst, int& nsecond, thetapoint& firstsex, thetapoint& secondsex, real& pf, real& ps, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale, ind& father, ind& mother, ind& child, int& fseg, ind& q, int& sstart, int& send, int& sseg, int& secondseg, int& firstseg, int& fstart, int& fend)
	{
		real segval;
		int first,second;
		
		/*segsexup*/
		initseg(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child);
		{
			//fbj	censorrec& with = *censorstruct;
			if (p->male)
			{
				{
					thisarray& with1 = *p->gen;
					for( first=1; first <= nfirst; first ++)
						if (with1.genarray[first]!=0.0)
						{
							segval=0.0;
							fseg=first;
							second=1;
							{ thisarray& with2 = *q->gen;
								for( second=1; second <= nsecond; second ++)
									if (with2.genarray[second]!=0.0)
									{
										sstart=probstart[second];
										send=probend[second];
										sseg=segstart[second];
										if (mutsys!=0)
											segval=segval+with2.genarray[second]*msegsex(p, fseg, secondseg, sseg, secondsex, sstart, send, ps, pf, firstseg, firstsex, fstart, fend);
										else segval=segval+with2.genarray[second]*segsex(p, secondseg, sseg, secondsex, sstart, send, fseg, firstseg, firstsex, fstart, fend);
									}
							}
							with1.genarray[first]=with1.genarray[first]*segval*segscale;
						}
				}
			}
			else { thisarray& with1 = *p->gen;
				for( first=1; first <= nfirst; first ++)
					if (with1.genarray[first]!=0.0)
					{
						segval=0.0;
						fstart=probstart[first];
						fend=probend[first];
						fseg=segstart[first];
						{ thisarray& with2 = *q->gen;
							for( second=1; second <= nsecond; second ++)
								if (with2.genarray[second]!=0.0)
								{
									sseg=second;
									if (mutsys!=0)
										segval=segval+with2.genarray[second]*msegsex(p, fseg, secondseg, sseg, secondsex, sstart, send, ps, pf, firstseg, firstsex, fstart, fend);
									else segval=segval+with2.genarray[second]*segsex(p, secondseg, sseg, secondsex, sstart, send, fseg, firstseg, firstsex, fstart, fend);
								}
						}
						with1.genarray[first]=with1.genarray[first]*segval*segscale;
					}
			}
		}
		cleanup(q);
		exitseg(nuscale, child, father, mother);
	}       /*segsexup*/
	
	
	static void segsexdown(ind& p, int& nfirst, int& nsecond, thetapoint& firstsex, thetapoint& secondsex, real& pf, real& ps, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale, ind& father, ind& mother, ind& child, int& fseg, ind& q, int& sstart, int& send, int& sseg, int& secondseg, int& firstseg, int& fstart, int& fend, ind& r)
	{
		unsigned short here;
		genotype gene;
		real val,temp2;
		int j,first,second;
		
		/*segsexdown*/
		initseg(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child);
		for( here=1; here <= fgeno; here ++) gene[here]=0.0;
		{
			gennurec& with = *gennustruct;
			{
				//fbj		censorrec& with1 = *censorstruct;
				{
					thisarray& with2 = *p->gen;
					for( first=1; first <= nfirst; first ++)
						if (with2.genarray[first]!=0.0)
						{
							fseg=first;
							second=1;
							{
								thisarray& with3 = *q->gen;
								for( second=1; second <= nsecond; second ++)
									if (with3.genarray[second]!=0.0)
									{
										val=with3.genarray[second]*p->gen->genarray[first]*segscale;
										sstart=probstart[second];
										send=probend[second];
										sseg=segstart[second];
										if (nchild!=0)  val=val*segsex(p, secondseg, sseg, secondsex, sstart, send, fseg, firstseg, firstsex, fstart, fend);
										if (val!=0.0)
										{
											secondseg=sseg;
											{
												thetavalues& with4 = *femaletheta;
												for( j=sstart; j <= send; j ++)
												{
													temp2=with4.segprob[j];
													if (temp2!=0.0)
													{	if (r->male)
													{
														here=seghap1[secondseg];
														gene[here]=gene[here]+temp2*val;
														here=seghap2[secondseg];
														gene[here]=gene[here]+temp2*val;
													}
													else {
														here=with.genenumber[seghap1[secondseg]][first];
														gene[here]=gene[here]+temp2*val;
														here=with.genenumber[seghap2[secondseg]][first];
														gene[here]=gene[here]+temp2*val;
													}
													}
													secondseg=secondseg+1;
												}
											}
										}
									}
							} /*second*/
						}
				}
			}
		}/*first*/
		p->gen->genarray=r->gen->genarray;
		r->gen->genarray=gene;
		{ thisarray& with = *r->gen;
			for( first=1; first <= fgeno; first ++)
				with.genarray[first]=with.genarray[first]*p->gen->genarray[first];}
		cleanup(p);
		cleanup(q);
		exitseg(nuscale, child, father, mother);
	}       /*segsexdown*/
	
	
	static void segctop(ind& p, int& nfirst, int& nsecond, thetapoint& firstsex, thetapoint& secondsex, real& pf, real& ps, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale, ind& father, ind& mother, ind& child, int& fstart, int& fend, int& fseg, ind& q, int& sstart, int& send, int& sseg, int& firstseg, int& secondseg)
	{
		real segval;
		int first,second,sstop;
		boolean thiscensor,thisrare;
		
		/*segctop*/
		initseg(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child);
		{
			censorrec& with = *censorstruct;
			{
				thisarray& with1 = *p->gen;
				for( first=1; first <= fgeno; first ++)
					if (with1.genarray[first]!=0.0)
					{
						segval=0.0;
						fstart=probstart[first];
						fend=probend[first];
						fseg=segstart[first];
						thisrare=rare[first];
						second=1;
						if (thisrare)
						{ thisarray& with2 = *q->gen;
							do {
								sstart=probstart[second];
								send=probend[second];
								sseg=segstart[second];
								sstop=second+(send-sstart)+1;
								if (thisc<maxcensor)
								{
									thisc=thisc+1;
									thiscensor=with.censor[thisc];
								}
								else thiscensor=(with2.genarray[second]==0.0)
									|| rare[second];
								if (thiscensor)
									second=sstop;
								else {
									if (mutsys!=0)
										segval=segval+with2.genarray[second]*msegfast(firstseg, fseg, send, sstart, firstsex, fstart, fend, secondseg, sseg, secondsex, ps, pf);
									else segval=segval+with2.genarray[second]*segfast(firstseg, fseg, send, sstart, firstsex, fstart, fend, secondseg, sseg, secondsex);
									second=sstop;
								}
							} while (!(second>fgeno));}
						else
						{ thisarray& with2 = *q->gen;
							do {
								sstart=probstart[second];
								send=probend[second];
								sseg=segstart[second];
								sstop=second+(send-sstart)+1;
								if (thisc<maxcensor)
								{
									thisc=thisc+1;
									thiscensor=with.censor[thisc];
								}
								else thiscensor=(with2.genarray[second]==0.0);
								if (thiscensor)
									second=sstop;
								else {
									if (mutsys!=0)
										segval=segval+with2.genarray[second]*msegfast(firstseg, fseg, send, sstart, firstsex, fstart, fend, secondseg, sseg, secondsex, ps, pf);
									else segval=segval+with2.genarray[second]*segfast(firstseg, fseg, send, sstart, firstsex, fstart, fend, secondseg, sseg, secondsex);
									second=sstop;
								}
							} while (!(second>fgeno));}
						with1.genarray[first]=with1.genarray[first]*segval*segscale;
					}
			}
		}
		cleanup(q);
		exitseg(nuscale, child, father, mother);
		if (approximate && firstapprox && (p->pa==nil))
			getapprox(p, thisped);
	}       /*segctop*/
	
	
	static void segtop(ind& p, int& nfirst, int& nsecond, thetapoint& firstsex, thetapoint& secondsex, real& pf, real& ps, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale, ind& father, ind& mother, ind& child, int& fstart, int& fend, int& fseg, ind& q, int& sstart, int& send, int& sseg, int& firstseg, int& secondseg)
	{
		real segval,val;
		int first,second,sstop;
		//	boolean thatrare;
		boolean thisrare;
		
		/*segtop*/
		initseg(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child);
		{
			censorrec& with = *censorstruct;
			{
				thisarray& with1 = *p->gen;
				for( first=1; first <= fgeno; first ++)
					if (with1.genarray[first]!=0.0)
					{
						thisrare=rare[first];
						segval=0.0;
						fstart=probstart[first];
						fend=probend[first];
						fseg=segstart[first];
						second=1;
						if (thisrare)
						{ thisarray& with2 = *q->gen;
							do {
								sstart=probstart[second];
								send=probend[second];
								sseg=segstart[second];
								sstop=second+(send-sstart)+1;
								if (thisc<=maxcensor)  thisc=thisc+1;
								if ((with2.genarray[second]==0.0) || rare[second])
								{
									second=sstop;
									if (thisc<=maxcensor)  with.censor[thisc]=true;
								}
								else {
									if (mutsys!=0)
										val=msegfast(firstseg, fseg, send, sstart, firstsex, fstart, fend, secondseg, sseg, secondsex, ps, pf);
									else val=segfast(firstseg, fseg, send, sstart, firstsex, fstart, fend, secondseg, sseg, secondsex);
									if (thisc<=maxcensor)  with.censor[thisc]=(val==0.0);
									segval=segval+with2.genarray[second]*val;
									second=sstop;
								}
							} while (!(second>fgeno));}
						else
						{
							thisarray& with2 = *q->gen;
							do {
								sstart=probstart[second];
								send=probend[second];
								sseg=segstart[second];
								sstop=second+(send-sstart)+1;
								if (thisc<=maxcensor)  thisc=thisc+1;
								if (with2.genarray[second]==0.0)
								{
									second=sstop;
									if (thisc<=maxcensor)  with.censor[thisc]=true;
								}
								else {
									if (mutsys!=0)
										val=msegfast(firstseg, fseg, send, sstart, firstsex, fstart, fend, secondseg, sseg, secondsex, ps, pf);
									else val=segfast(firstseg, fseg, send, sstart, firstsex, fstart, fend, secondseg, sseg, secondsex);
									if (thisc<=maxcensor)  with.censor[thisc]=(val==0.0);
									segval=segval+with2.genarray[second]*val;
									second=sstop;
								}
							} while (!(second>fgeno));}
						with1.genarray[first]=with1.genarray[first]*segval*segscale;
					}
			}
		}
		cleanup(q);
		exitseg(nuscale, child, father, mother);
		if (approximate && firstapprox && (p->pa==nil))
			getapprox(p, thisped);
	}       /*segtop*/
	
	
	static void segcapprox(ind& p, int& nfirst, int& nsecond, thetapoint& firstsex, thetapoint& secondsex, real& pf, real& ps, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale, ind& father, ind& mother, ind& child, int& fstart, int& fend, int& fseg, ind& q, int& sstart, int& send, int& sseg, int& firstseg, int& secondseg)
	{
		real segval;
		int first,second,sstop;
		boolean thiscensor;
		
		/*segcapprox*/
		initseg(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child);
		{
			censorrec& with = *censorstruct;
			{
				approxrec& with1 = *approxstruct;
				{
					thisarray& with2 = *p->gen;
					for( first=1; first <= fgeno; first ++)
						if (with2.genarray[first]!=0.0)
						{
							segval=0.0;
							fstart=probstart[first];
							fend=probend[first];
							fseg=segstart[first];
							second=1;
							{
								thisarray& with3 = *q->gen;
								do {
									sstart=probstart[second];
									send=probend[second];
									sseg=segstart[second];
									sstop=second+(send-sstart)+1;
									if (thisc<maxcensor)
									{
										thisc=thisc+1;
										thiscensor=with.censor[thisc];
									}
									else thiscensor=(with3.genarray[second]==0.0);
									if (thiscensor || ! with1.approxarray[thisped][first])
										second=sstop;
									else {
										if (mutsys!=0)
											segval=segval+with3.genarray[second]*msegfast(firstseg, fseg, send, sstart, firstsex, fstart, fend, secondseg, sseg, secondsex, ps, pf);
										else segval=segval+with3.genarray[second]*segfast(firstseg, fseg, send, sstart, firstsex, fstart, fend, secondseg, sseg, secondsex);
										second=sstop;
									}
								} while (!(second>fgeno));}
							with2.genarray[first]=with2.genarray[first]*segval*segscale;
						}
				}
			}
		}
		cleanup(q);
		exitseg(nuscale, child, father, mother);
	}       /*segcapprox*/
	
	
	static void segup(ind& p, int& nfirst, int& nsecond, thetapoint& firstsex, thetapoint& secondsex, real& pf, real& ps, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale, ind& father, ind& mother, ind& child, int& fstart, int& fend, int& fseg, ind& q, int& sstart, int& send, int& sseg, int& firstseg, int& secondseg)
	{
		real segval;
		int first,second;
		
		/*segup*/
		initseg(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child);
		{
			//fbj	censorrec& with = *censorstruct;
			{
				thisarray& with1 = *p->gen;
				for( first=1; first <= fgeno; first ++)
					if (with1.genarray[first]!=0.0)
					{
						segval=0.0;
						fstart=probstart[first];
						fend=probend[first];
						fseg=segstart[first];
						{ thisarray& with2 = *q->gen;
							for( second=1; second <= fgeno; second ++)
								if (with2.genarray[second]!=0.0)
								{
									sstart=probstart[second];
									send=probend[second];
									sseg=segstart[second];
									if (mutsys!=0)
										segval=segval+with2.genarray[second]*msegfun(firstseg, fseg, firstsex, fstart, fend, secondseg, sseg, secondsex, sstart, send, ps, pf);
									else segval=segval+with2.genarray[second]*segfun(firstseg, fseg, firstsex, fstart, fend, secondseg, sseg, secondsex, sstart, send);
								}
						}
						with1.genarray[first]=with1.genarray[first]*segval*segscale;
					}
			}
		}
		cleanup(q);
		exitseg(nuscale, child, father, mother);
	}       /*segup*/
	
	
	static void segdown(ind& p, int& nfirst, int& nsecond, thetapoint& firstsex, thetapoint& secondsex, real& pf, real& ps, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale, ind& father, ind& mother, ind& child, int& fstart, int& fend, int& fseg, ind& q, int& sstart, int& send, int& sseg, int& firstseg, int& secondseg, ind& r)
	{
		unsigned short here;
		genotype gene;
		real val,temp1,temp2;
		int f1,f2,s1,s2,i,j,first,second;
		
		
		initseg(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child);
		for( here=1; here <= fgeno; here ++) gene[here]=0.0;
		{
			gennurec& with = *gennustruct;
			{
				//fbj		censorrec& with1 = *censorstruct;
				{
					thisarray& with2 = *p->gen;
					for( first=1; first <= mgeno; first ++)
						if (with2.genarray[first]!=0.0)
						{
							fstart=probstart[first];
							fend=probend[first];
							fseg=segstart[first];
							{
								thisarray& with3 = *q->gen;
								for( second=1; second <= fgeno; second ++)
									if (with3.genarray[second]!=0.0)
									{
										sstart=probstart[second];
										send=probend[second];
										sseg=segstart[second];
										val=segscale*with3.genarray[second]*p->gen->genarray[first];
										if (nchild!=0)  val=segfun(firstseg, fseg, firstsex, fstart, fend, secondseg, sseg, secondsex, sstart, send)*val;
										if (val!=0.0)
										{
											firstseg=fseg;
											{ thetavalues& with4 = *maletheta;
												for( i=fstart; i <= fend; i ++)
												{
													temp1=with4.segprob[i];
													secondseg=sseg;
													if (temp1!=0.0)
													{ thetavalues& with5 = *femaletheta;
														for( j=sstart; j <= send; j ++)
															if (with5.segprob[j]==0.0)
																secondseg=secondseg+1;
															else
															{
																temp2=with5.segprob[j]*temp1*val;
																f1=seghap1[firstseg];
																f2=seghap2[firstseg];
																s1=seghap1[secondseg];
																s2=seghap2[secondseg];
																here=with.genenumber[s1][f1];
																gene[here]=gene[here]+temp2;
																here=with.genenumber[s1][f2];
																gene[here]=gene[here]+temp2;
																here=with.genenumber[s2][f1];
																gene[here]=gene[here]+temp2;
																here=with.genenumber[s2][f2];
																gene[here]=gene[here]+temp2;
																secondseg=secondseg+1;
															}}
													firstseg=firstseg+1;
												}
											}
										}
									}
							} /*second*/
						}
				}
			}
		}/*first*/
		p->gen->genarray=r->gen->genarray;
		r->gen->genarray=gene;
		{ thisarray& with = *r->gen;
			for( first=1; first <= fgeno; first ++)
				with.genarray[first]=with.genarray[first]*p->gen->genarray[first];}
		cleanup(p);
		cleanup(q);
		exitseg(nuscale, child, father, mother);
	}       /*segdown*/
	
	
	static void msegsexdown(ind& p, int& nfirst, int& nsecond, thetapoint& firstsex, thetapoint& secondsex, real& pf, real& ps, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale, ind& father, ind& mother, ind& child, int& fseg, ind& q, int& sstart, int& send, int& sseg, int& secondseg, int& firstseg, int& fstart, int& fend, ind& r)
	{
		unsigned short here;
		genotype gene;
		real val,temp2;
		int ms2,ms1,mf,j,first,second;
		
		/*msegsexdown*/
		initseg(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child);
		for( here=1; here <= fgeno; here ++) gene[here]=0.0;
		{
			gennurec& with = *gennustruct;
			{
				//fbj		censorrec& with1 = *censorstruct;
				{
					thisarray& with2 = *p->gen;
					for( first=1; first <= nfirst; first ++)
						if (with2.genarray[first]!=0.0)
						{
							fseg=first;
							second=1;
							{
								thisarray& with3 = *q->gen;
								for( second=1; second <= nsecond; second ++)
									if (with3.genarray[second]!=0.0)
									{
										val=with3.genarray[second]*p->gen->genarray[first]*segscale;
										sstart=probstart[second];
										send=probend[second];
										sseg=segstart[second];
										if (nchild!=0)  val=val*msegsex(p, fseg, secondseg, sseg, secondsex, sstart, send, ps, pf, firstseg, firstsex, fstart, fend);
										if (val!=0.0)
										{
											mf=muthap[seghap1[fseg]];
											secondseg=sseg;
											{ thetavalues& with4 = *femaletheta;
												for( j=sstart; j <= send; j ++)
												{
													ms1=muthap[seghap1[secondseg]];
													ms2=muthap[seghap2[secondseg]];
													temp2=with4.segprob[j];
													if (temp2!=0.0)
													{	if (r->male)
													{
														here=seghap1[secondseg];
														gene[here]=gene[here]+(1-ps)*temp2*val;
														here=seghap2[secondseg];
														gene[here]=gene[here]+(1-ps)*temp2*val;
														here=ms1;
														gene[here]=gene[here]+ps*temp2*val;
														here=ms2;
														gene[here]=gene[here]+ps*temp2*val;
													}
													else {
														here=with.genenumber[seghap1[secondseg]][first];
														gene[here]=gene[here]+(1-pf)*(1-ps)*temp2*val;
														here=with.genenumber[seghap2[secondseg]][first];
														gene[here]=gene[here]+(1-pf)*(1-ps)*temp2*val;
														here=with.genenumber[seghap1[secondseg]][mf];
														gene[here]=gene[here]+(pf)*(1-ps)*temp2*val;
														here=with.genenumber[seghap2[secondseg]][mf];
														gene[here]=gene[here]+(pf)*(1-ps)*temp2*val;
														here=with.genenumber[ms1][first];
														gene[here]=gene[here]+(1-pf)*ps*temp2*val;
														here=with.genenumber[ms2][first];
														gene[here]=gene[here]+(1-pf)*ps*temp2*val;
														here=with.genenumber[ms1][mf];
														gene[here]=gene[here]+pf*ps*temp2*val;
														here=with.genenumber[ms2][mf];
														gene[here]=gene[here]+pf*ps*temp2*val;
													}
													}
													secondseg=secondseg+1;
												}
											}
										}
									}
							} /*second*/
						}
				}
			}
		}/*first*/
		p->gen->genarray=r->gen->genarray;
		r->gen->genarray=gene;
		{ thisarray& with = *r->gen;
			for( first=1; first <= fgeno; first ++)
				with.genarray[first]=with.genarray[first]*p->gen->genarray[first];}
		cleanup(p);
		cleanup(q);
		exitseg(nuscale, child, father, mother);
	}       /*msegsexdown*/
	
	
	static void msegdown(ind& p, int& nfirst, int& nsecond, thetapoint& firstsex, thetapoint& secondsex, real& pf, real& ps, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale, ind& father, ind& mother, ind& child, int& fstart, int& fend, int& fseg, ind& q, int& sstart, int& send, int& sseg, int& firstseg, int& secondseg, ind& r)
	{
		unsigned short here;
		genotype gene;
		real val,temp,temp1,temp2;
		int i,j,first,second,f1,f2,s1,s2,mf1,mf2,ms1,ms2;
		
		/*msegdown*/
		initseg(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child);
		for( here=1; here <= fgeno; here ++) gene[here]=0.0;
		{
			gennurec& with = *gennustruct;
			{
				//fbj		censorrec& with1 = *censorstruct;
				{
					thisarray& with2 = *p->gen;
					for( first=1; first <= fgeno; first ++)
						if (with2.genarray[first]!=0.0)
						{
							fstart=probstart[first];
							fend=probend[first];
							fseg=segstart[first];
							{ thisarray& with3 = *q->gen;
								for( second=1; second <= fgeno; second ++)
									if (with3.genarray[second]!=0.0)
									{
										sstart=probstart[second];
										send=probend[second];
										sseg=segstart[second];
										val=with3.genarray[second]*p->gen->genarray[first]*segscale;
										if (nchild!=0)  val=msegfun(firstseg, fseg, firstsex, fstart, fend, secondseg, sseg, secondsex, sstart, send, ps, pf)*val;
										if (val!=0.0)
										{
											firstseg=fseg;
											{ thetavalues& with4 = *maletheta;
												for( i=fstart; i <= fend; i ++)
												{
													temp1=with4.segprob[i];
													secondseg=sseg;
													if (temp1!=0.0)
													{ thetavalues& with5 = *femaletheta;
														for( j=sstart; j <= send; j ++)
															if (with5.segprob[j]==0.0)
																secondseg=secondseg+1;
															else
																
															{
																temp2=with5.segprob[j];
																f1=seghap1[firstseg];
																f2=seghap2[firstseg];
																s1=seghap1[secondseg];
																s2=seghap2[secondseg];
																temp=(1-pf)*(1-ps)*temp1*temp2*val;
																here=with.genenumber[s1][f1];
																gene[here]=gene[here]+temp;
																here=with.genenumber[s1][f2];
																gene[here]=gene[here]+temp;
																here=with.genenumber[s2][f1];
																gene[here]=gene[here]+temp;
																here=with.genenumber[s2][f2];
																gene[here]=gene[here]+temp;
																ms1=muthap[s1];
																ms2=muthap[s2];
																temp=(1-pf)*ps*temp1*temp2*val;
																here=with.genenumber[ms1][f1];
																gene[here]=gene[here]+temp;
																here=with.genenumber[ms1][f2];
																gene[here]=gene[here]+temp;
																here=with.genenumber[ms2][f1];
																gene[here]=gene[here]+temp;
																here=with.genenumber[ms2][f2];
																gene[here]=gene[here]+temp;
																mf1=muthap[f1];
																mf2=muthap[f2];
																temp=pf*(1-ps)*temp1*temp2*val;
																here=with.genenumber[mf1][s1];
																gene[here]=gene[here]+temp;
																here=with.genenumber[mf1][s2];
																gene[here]=gene[here]+temp;
																here=with.genenumber[mf2][s1];
																gene[here]=gene[here]+temp;
																here=with.genenumber[mf2][s2];
																gene[here]=gene[here]+temp;
																temp=pf*ps*temp1*temp2*val;
																here=with.genenumber[mf1][ms1];
																gene[here]=gene[here]+temp;
																here=with.genenumber[mf1][ms2];
																gene[here]=gene[here]+temp;
																here=with.genenumber[mf2][ms1];
																gene[here]=gene[here]+temp;
																here=with.genenumber[mf2][ms2];
																gene[here]=gene[here]+temp;
																secondseg=secondseg+1;
															}
													}
													firstseg=firstseg+1;
												}
											}
										}
									}
							} /*second*/
						}
				}
			}
		}/*first*/
		p->gen->genarray=r->gen->genarray;
		r->gen->genarray=gene;
		{
			thisarray& with = *r->gen;
			for( first=1; first <= fgeno; first ++)
				with.genarray[first]=with.genarray[first]*p->gen->genarray[first];
		}
		cleanup(p);
		cleanup(q);
		exitseg(nuscale, child, father, mother);
	}       /*msegdown*/
	
	
	
	
	static void seg(ind& p,ind& q,ind& r,direction peel, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale)
	
	{
		ind child,father,mother;
		boolean phaseunkn;
		int fseg,sseg,sstart,send,fstart,fend,nfirst,nsecond,firstseg,secondseg;
		real pf,ps;
		thetapoint firstsex,secondsex;
		
		
		/*seg*/
		phaseunkn=((q->pa==nil) && (q->firstpass) && (q->inloop==0) && (! disequi));
		if (p->male)
		{
			father=p;
			mother=q;
		}
		else {
			father=q;
			mother=p;
		}
		if (peel==peelup)
			if ((p->pa==nil) && phaseunkn && approximate)
				if (firstapprox && firsttime)
					switch (thispath) {
						case auto_:segtop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg); break;
						case mauto:segtop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg); break;
						case sex:segsextop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, fstart, fend, firstseg); break;
						case msex:segsextop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, fstart, fend, firstseg); break;
						case last_pathway: exit_error("last_pathway occured (fbj)");
					}
				else/*first approximate not first time*/
					if (firstapprox)
						switch (thispath) {
							case auto_:segctop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg); break;
							case mauto:segctop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg); break;
							case sex:segsexctop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, fstart, fend, firstseg); break;
							case msex:segsexctop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, fstart, fend, firstseg); break;
							case last_pathway: exit_error("last_pathway occured (fbj)");
						}
					else/*approximate*/
						switch (thispath) {
							case auto_:segcapprox(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg); break;
							case mauto:segcapprox(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg); break;
							case sex:segsexctop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, fstart, fend, firstseg); break;
							case msex:segsexctop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, fstart, fend, firstseg); break;
							case last_pathway: exit_error("last_pathway occured (fbj)");
						}
					else/*do not approximate*/
						if (phaseunkn)
							if (firsttime)
								switch (thispath) {
									case auto_:segtop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg); break;
									case mauto:segtop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg); break;
									case sex:segsextop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, fstart, fend, firstseg); break;
									case msex:segsextop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, fstart, fend, firstseg); break;
									case last_pathway: exit_error("last_pathway occured (fbj)");
								}
							else/*not firsttime*/
								switch (thispath) {
									case auto_:segctop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg); break;
									case mauto:segctop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg); break;
									case sex:segsexctop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, fstart, fend, firstseg); break;
									case msex:segsexctop(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, fstart, fend, firstseg); break;
									case last_pathway: exit_error("last_pathway occured (fbj)");
								}
							else/*phaseinfo*/
								switch (thispath) {
									case auto_:segup(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg); break;
									case mauto:segup(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg); break;
									case sex:segsexup(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, firstseg, fstart, fend); break;
									case msex:segsexup(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, firstseg, fstart, fend); break;
									case last_pathway: exit_error("last_pathway occured (fbj)");
								}
							else/*not peelup*/
								switch (thispath) {
									case auto_:segdown(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg, r); break;
									case mauto:msegdown(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fstart, fend, fseg, q, sstart, send, sseg, firstseg, secondseg, r); break;
									case sex:segsexdown(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, firstseg, fstart, fend, r); break;
									case msex:msegsexdown(p, nfirst, nsecond, firstsex, secondsex, pf, ps, thisped, loopgen, holdpoint, nuscale, father, mother, child, fseg, q, sstart, send, sseg, secondseg, firstseg, fstart, fend, r); break;
									case last_pathway: exit_error("last_pathway occured (fbj)");
								}
		q->firstpass=false;
		p->firstpass=false;
	}       /*seg*/
	
	static void collapsedown(ind p, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale);
	static void collapseup(ind p, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale);
	
	static void collapseup(ind p, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale)
	{
		ind q,child,nextchild;
		boolean down;
		
		/*collapseup*/
		p->done=true;
		if (p->foff!=nil)
		{
			down=false;
			child=p->foff;
			while (child!=nil)
			{
				down=false;
				if (p->male)
					q=child->ma;
				else q=child->pa;
				if (! q->done)
				{
					collapsedown(q, thisped, loopgen, holdpoint, nuscale);
					nextchild=child;
					while (nextchild!=nil)
					{
						if ((nextchild->pa==q) || (nextchild->ma==q))
						{
							if (!nextchild->up)
								collapseup(nextchild, thisped, loopgen, holdpoint, nuscale);
							else
								down=true;
						}
						if (p->male)	nextchild=nextchild->nextpa;
						else			nextchild=nextchild->nextma;
					}
					if (q->multi)	collapseup(q, thisped, loopgen, holdpoint, nuscale);
					if (!down)		seg(p,q,child,peelup, thisped, loopgen, holdpoint, nuscale);
					else			collapsedown(p, thisped, loopgen, holdpoint, nuscale);
				}
				if (p->male)	child=child->nextpa;
				else			child=child->nextma;
			}
		}
	}       /*collapseup*/
	
	static void collapsedown(ind p, int& thisped, array1d<1,maxloop,int>& loopgen, array1d<1,maxloop,genpoint>& holdpoint, int& nuscale)
	{     /*collapsedown*/
		if (p->pa!=nil)
		{
			p->up=true;
			collapseup(p->pa, thisped, loopgen, holdpoint, nuscale);
			seg(p->pa,p->ma,p,peeldown, thisped, loopgen, holdpoint, nuscale);
		}
	}       /*collapsedown*/
	
	static void riskcumul(ind& proband, real& hetero, real& homo)
	
	{
		int i;
		
		/*riskcumul*/
		{
			thisarray& with = *proband->gen;
			if (sexlink && proband->male)
			{
				for( i=1; i <= mgeno; i ++)
					if (riskmale[i])  hetero=hetero+with.genarray[i];
			}
			else
			{
				for( i=1; i <= fgeno; i ++)
					if (risk2[i])
						homo=homo+with.genarray[i];
					else if (risk1[i])
						hetero=hetero+with.genarray[i];
			}
		}
	}       /*riskcumul*/
	
	
	static void riskcalc(real& homo, real& hetero, ind& proband)
	{
		real normal;
		
		/*riskcalc*/
		homo=homo/like;
		hetero=hetero/like;
		normal=1-homo-hetero;
		if (dooutput)
		{
			std::cout << "RISK FOR PERSON " << fmt(proband->id,6) << " IN PEDIGREE " << fmt(proband->ped,7) << std::endl;
			if ((! proband->male) || (! sexlink))
				std::cout << "HOMOZYGOTE CARRIER   : " << fmt(homo,8,5) << std::endl;
			if ((! proband->male) || (! sexlink))
				std::cout << "HETEROZYGOTE CARRIER : " << fmt(hetero,8,5) << std::endl;
			else std::cout << "MALE CARRIER         : " << fmt(hetero,8,5) << std::endl;
			std::cout << "NORMAL               : " << fmt(normal,8,5) << std::endl;
		}
	}       /*riskcalc*/
	
	void likelihood(int thisped,ind proband)
	{
		array1d<1,maxloop,int> loopmax,loopgen;
		real homo,hetero,tmplike;
		int i,j,nuscale,startscale;
		array1d<1,maxloop,genpoint> holdpoint;
		boolean gocalc,alldone;
		
		/*likelihood*/
		if (informative[thisped])
		{
			homo=0.0;
			hetero=0.0;
			tmplike=0.0;
			alldone=false;
			startscale=0;
			for( i=1; i <= maxloop; i ++)
			{
				loopgen[i]=1;
				loopmax[i]=1;
				holdpoint[i]=nil;
				if (looppers[thisped][i][1]!=nil)
				{
					thisperson& with = *looppers[thisped][i][1];
					
					with.gen = new thisarray;
					{ thisarray& with1 = *with.gen;
						for( j=1; j <= fgeno; j ++) with1.genarray[j]=0.0;}
					getvect(looppers[thisped][i][1]);
					if (looppers[thisped][i][1]->pa==nil)  startscale=startscale+1;
					holdpoint[i]=with.gen;
					with.gen=nil;
					if (with.male)
						loopmax[i]=mgeno;
					else loopmax[i]=fgeno;
				}
			}
			loopgen[1]=0;
			do {
				i=1;
				do {
					loopgen[i]=loopgen[i]+1;
					if (loopgen[i]>loopmax[i])
						loopgen[i]=1;
					else i=maxloop;
					i=i+1;
				} while (!(i>maxloop));
				gocalc=true;
				for( i=1; i <= maxloop; i ++)
				/*ML change*/
					if (holdpoint[i]!=nil)
						if (holdpoint[i]->genarray[loopgen[i]]==0.0)
							gocalc=false;
				if (gocalc)
				{
					for( i=1; i <= totperson; i ++)
					{
						thisperson& with = *person[i];
						
						with.gen=nil;
						with.done=false;
						with.up=false;
					}
					nuscale=startscale;
					collapseup(proband, thisped, loopgen, holdpoint, nuscale);
					collapsedown(proband, thisped, loopgen, holdpoint, nuscale);
					like=0.0;
					{ thisarray& with = *proband->gen;
						if (proband->male)
							for( i=1; i <= mgeno; i ++)
								like=like+with.genarray[i];
						else for( i=1; i <= fgeno; i ++)
							like=like+with.genarray[i];}
					tmplike=tmplike+like;
					if ((risk) && (like!=0.0))  riskcumul(proband, hetero, homo);
					for( i=1; i <= totperson; i ++)
						if (person[i]->gen!=nil)  cleanup(person[i]);
				}
				alldone=true;
				for( i=1; i <= maxloop; i ++) alldone=alldone && (loopgen[i]==loopmax[i]);
			} while (!alldone);
			like=tmplike;
			if (risk && (like!=0.0))  riskcalc(homo, hetero, proband);
			if (like==0.0)
				like=zerolike;
			else like=log(like)-nuscale*log(segscale);
			for( i=1; i <= maxloop; i ++)
				if (holdpoint[i]!=nil)
				{
					delete holdpoint[i];
					holdpoint[i]=nil;
				}
		}
		else like=0.0;
	}      /*likelihood*/
	
	
	void checkzero()
	{
		if (! firsttime)
		{
			{ thetavalues& with = *maletheta;
				for( i=1; i <= nlocus-1; i ++)
					if ((with.theta[i]!=0.0) && zeromale[i])  firsttime=true;}
			{ thetavalues& with = *femaletheta;
				for( i=1; i <= nlocus-1; i ++)
					if ((with.theta[i]!=0.0) && zerofemale[i])  firsttime=true;}
		}
		if (maletheta->theta[which]==0.0)  firsttime=true;
		if (femaletheta->theta[which]==0.0)  firsttime=true;
		if (firsttime)
		{
			{ thetavalues& with = *maletheta;
				for( i=1; i <= nlocus-1; i ++)
					zeromale[i]=with.theta[i]==0.0;}
			{ thetavalues& with = *femaletheta;
				for( i=1; i <= nlocus-1; i ++)
					zerofemale[i]=with.theta[i]==0.0;}
		}
	}
	
	struct HLOD_result {
		double LOD,alpha,HLOD;
		HLOD_result():LOD(0),alpha(0),HLOD(0){}
		HLOD_result(double l, double a, double h):LOD(l),alpha(a),HLOD(h){}
	};
	
	inline HLOD_result cal_HLOD(const std::vector<double>& lodvec)
	{
		double minAlpha=0, maxAlpha=1+std::numeric_limits<double>::epsilon()*20, LOD=0;
		size_t size=lodvec.size();
		for (size_t i=0;i<size;i++) LOD+=lodvec[i];
		double maxhlod=0;
		double maxa=0;
		for (;;)
		{
			double incremnt = (maxAlpha-minAlpha)/20;
			for (double a=minAlpha; a<=maxAlpha; a+=incremnt)
			{
				double hlod=0;
				for (unsigned i=0;i<size;i++) hlod += log10( a*pow(10,lodvec[i]) + 1-a );
				if (hlod>maxhlod) { maxhlod=hlod; maxa=a; }
			}
			if (incremnt<=0.001) break;
			minAlpha = std::max(minAlpha, maxa - incremnt);
			maxAlpha = std::min(maxAlpha, maxa + incremnt);
		}
		return HLOD_result(LOD,maxa,maxhlod);
	}
	
	inline double cal_PLOD (const std::vector<double>& lodvec)
	{
		double plod = 0;
		size_t size=lodvec.size();
		for (size_t i=0;i<size;i++) if (lodvec[i]>0) plod+=lodvec[i];
		return plod;
	}
	
	void iterpeds(double& log10L, double& lod, double& hlod, double& alpha, double& plod)
	{
		int i;
		unsigned char thisped;
		
		
		tlike=0.0;
		alike=0.0;
		for( i=1; i <= totperson; i ++)
			person[i]->done=false;
		for( i=1; i <= totperson; i ++)
			person[i]->firstpass=true;
		thisc=minint;
		recombination();
		checkzero();
		
		if (dooutput)
		{
			for( i=1; i <= 35; i ++) std::cout << '-';
			std::cout << std::endl;
			for( i=1; i <= 35; i ++) std::cout << '-';
			std::cout << std::endl;
			if (sexdif)
				std::cout << "MALE THETAS   ";
			else
				std::cout << "THETAS ";
			for( i=1; i <= nlocus-1; i ++) std::cout << fmt(maletheta->theta[i],6,3);
			if (interfer)
				std::cout << fmt(maletheta->theta[nlocus],6,3);
			std::cout << std::endl;
			if (sexdif)
			{
				std::cout << "FEMALE THETAS ";
				for( i=1; i <= nlocus-1; i ++) std::cout << fmt(femaletheta->theta[i],6,3);
				if (interfer)
					std::cout << fmt(femaletheta->theta[nlocus],6,3);
				std::cout << std::endl;
			}
			for( i=1; i <= 35; i ++) std::cout << '-';
			std::cout << std::endl;
			std::cout << "PEDIGREE |  LN LIKE  | LOG 10 LIKE" << std::endl;
			for( i=1; i <= 35; i ++) std::cout << '-';
			std::cout << std::endl;
		}
		std::vector<real> LodPerFam;
		for( thisped=1; thisped <= nuped; thisped ++)
		{
			likelihood(thisped,proband[thisped]);
			if (dooutput)
				std::cout << fmt(proband[thisped]->ped,9) << ' ' << fmt(like,12,6) << ' ';
			alike=alike+like;
			like=like/log_10;
			if (dooutput)
				std::cout << fmt(like,12,6) << std::endl;
			tlike=tlike+like;
			if (maletheta->theta[which]==0.5)	log10like_theta5[thisped]=like;
			else								LodPerFam.push_back(like-log10like_theta5[thisped]);
		}
		log10L=tlike;
		if (dooutput)
		{
			for( i=1; i <= 35; i ++) std::cout << '-';
			std::cout << std::endl;
			std::cout << "TOTALS    " << fmt(alike,12,6) << ' ' << fmt(tlike,12,6) << std::endl;
		}
		alike=-2*alike;
		if (dooutput) std::cout << "-2 LN(LIKE) = " << alike;
		if (score && ! risk)
		{
			if (nlocus==2)
			{
				if (maletheta->theta[which]==0.5)
					scorevalue=tlike;
				alike=tlike-scorevalue;
				lod=alike;
				if (dooutput)
					std::cout << " LOD SCORE = " << fmt(alike,12,6);
				if (maletheta->theta[which]!=0.5)
				{
					HLOD_result r = cal_HLOD(LodPerFam);
					alpha=r.alpha;
					hlod =r.HLOD;
					plod =cal_PLOD(LodPerFam);
					if (dooutput)
						std::cout << " HLOD = " << fmt(hlod,12,6) << " alpha = " << fmt(alpha,12,6) ;
				}
			}
			else
			{
				if (maletheta->theta[which]==0.5)
					scorevalue=alike;
				alike=scorevalue-alike;
				lod=alike;
				if (dooutput)
					std::cout << " LOG LIKE DIFFERENCE = " << fmt(alike,12,6);
			}
		}
		if (dooutput)
			std::cout << std::endl;
		if (firsttime)
		{
			if (thisc<maxcensor)
			{
				if (dooutput) std::cout << "Maxcensor can be reduced to " << thisc << std::endl;
			}
			else
			{	if (thisc>maxcensor)
				if (dooutput) std::cout << "you may gain efficiency by increasing maxcensor" << std::endl;
			}
		}
		firsttime=false;
	}
	
	
	void initmlink()
	{
		if (dooutput)
		{
			std::cout << "Copyright (C) CEPH, University of Utah, and Columbia University 1990" << std::endl;
			std::cout << "Program MLINK version" << fmt(version,6,2) << std::endl;
			std::cout << std::endl;
			std::cout << "The program constants are set to the following maxima:" << std::endl;
			std::cout << fmt(maxlocus,6) << " loci in mapping problem (maxlocus)" << std::endl;
			std::cout << fmt(maxall,6) << " alleles at a single locus (maxall)" << std::endl;
			std::cout << fmt(maxneed,6) << " recombination probabilities (maxneed)" << std::endl;
			std::cout << fmt(maxcensor,6) << " maximum of censoring array (maxcensor)" << std::endl;
			std::cout << fmt(maxhap,6) << " haplotypes = n1 x n2 x ... where ni = current # alleles at locus i" << std::endl;
			std::cout << fmt(maxfem,6) << " joint genotypes for a female" << std::endl;
			std::cout << fmt(maxmal,6) << " joint genotypes for a male" << std::endl;
			std::cout << fmt(maxind,6) << " individuals in all pedigrees combined (maxind)" << std::endl;
			std::cout << fmt(maxped,6) << " pedigrees (maxped)" << std::endl;
			std::cout << fmt(maxfact,6) << " binary codes at a single locus (maxfact)" << std::endl;
			std::cout << fmt(maxtrait,6) << " quantitative factor(s) at a single locus" << std::endl;
			std::cout << fmt(maxliab,6) << " liability classes" << std::endl;
			std::cout << fmt(maxfact,6) << " binary codes at a single locus" << std::endl;
			std::cout << fmt(scale,8,2) << " base scaling factor for likelihood (scale)" << std::endl;
			std::cout << fmt(scalemult,8,2) << " scale multiplier for each locus (scalemult)" << std::endl;
			std::cout << fmt(minfreq,8,5) << " frequency for elimination of heterozygotes (minfreq)" << std::endl;
			if (minfreq!=0.0)
			{
				std::cout << "IMPORTANT : RECOMPILE THIS PROGRAM WITH MINFREQ=0.0" << std::endl;
				std::cout << "FOR THE ANALYSIS OF RECESSIVE TRAITS" << std::endl;
			}
			std::cout << std::endl;
#if __cplusplus >= 201103L
			std::cout << "Length of real variables = " << sizeof(__typeof__(like)) << " bytes" << std::endl;
#elif defined __GNUG__
			std::cout << "Length of real variables = " << sizeof(typeof(like)) << " bytes" << std::endl;
#else
	#error This library requires either a GNU C++ compiler or C++0x/C++11 features.
#endif

			
		}
	}
	
} // end namespace

void mlink_common(std::istream& datafile, std::istream& ipedfile, std::istream& speedfile, double& log10L)
{
	initmlink();
	inputdata(datafile, ipedfile, speedfile);
	firsttime=true;
	lasttime=false;
	dolod=false;
	censorstruct = new censorrec;
	gennustruct = new gennurec;
	firstapprox=true;
	getlocations();
	
	scorevalue=0.0;
	log10like_theta5.assign(nuped+1, 0);
	double lod, hlod, alpha, plod;
	if (score && ! risk)
	{
		// likelihood for theta=0.5 is necessary for LOD calculation on other theta values.
		holdtheta=maletheta->theta[which];
		maletheta->theta[which]=0.5;
		iterpeds(log10L, lod, hlod, alpha, plod);
		maletheta->theta[which]=holdtheta;
	}
}

void mlink_calculate(std::istream& datafile, std::istream& ipedfile, std::istream& speedfile, double& theta, double& alpha, double& result)
{
	try
	{
		mlink_common(datafile, ipedfile, speedfile, result);
		if (linkage::defaultResultType==linkage::LOG10L && theta==0.5) return;
		double log10L, lod=0, hlod=0, hlod_alpha=1, plod=0;
		do {
			while (maletheta->theta[which]<=finish)
			{
				iterpeds(log10L, lod, hlod, hlod_alpha, plod);
				theta=maletheta->theta[which];
				alpha=hlod_alpha;
				if		(linkage::defaultResultType == linkage::LOD)	result=lod;
				else if (linkage::defaultResultType == linkage::HLOD)	result=hlod;
				else if (linkage::defaultResultType == linkage::PLOD)	result=plod;
				else if	(linkage::defaultResultType == linkage::LOG10L)	result=log10L;
				else exit(1);
				maletheta->theta[which]=maletheta->theta[which]+inc1;
			}
			
			if (! extra_lines.empty())
			{
				maletheta->theta[which] = extra_lines.front().t;
				inc1 = extra_lines.front().i;
				finish = extra_lines.front().f;
				extra_lines.pop();
			}
			else
				which=0;
			finish=finish+0.0001;
		} while (!(which==0));
	}
	catch (int e)
	{
		// There is an error/warning code e.
		theta = std::numeric_limits<double>::signaling_NaN();
		alpha = std::numeric_limits<double>::signaling_NaN();
		result = std::numeric_limits<double>::signaling_NaN();
	}
}

#ifdef LINKAGE_OPTIMIZE

#include <nlopt.hpp>

double stepwise_opt_max(double (*functionPtr)(const std::vector<double>&, std::vector<double>&, void*), std::vector<double>& x, double& lb, double& ub, double ms)
{
	double max_y = -DBL_MAX;
	double max_a = -DBL_MAX;
	std::vector<double> grad;
	void* f_data=NULL;
	for (;;)
	{
		double sl = (ub-lb)/50;
		for (double a=lb; a<=ub; )
		{
			x[0]=a;
			double y = (*functionPtr)(x, grad, f_data);
			if (y>max_y) { max_y=y; max_a=a; }
			if (a==ub) break;
			else if (a+sl<ub) a+=sl;
			else a=ub;
		}
		if (sl<=ms) break;
		lb = std::max(lb, max_a - sl);
		ub = std::min(ub, max_a + sl);
	}
	return max_y;
}

double _mlink_opt_theta(const std::vector<double> &x, std::vector<double> &grad, void* f_data)
{
	double log10L, lod, hlod, alpha, plod;
	maletheta->theta[which] = x[0];
	iterpeds(log10L, lod, hlod, alpha, plod);
	if		(linkage::defaultResultType == linkage::LOD)  return lod;
	else if (linkage::defaultResultType == linkage::HLOD) return hlod;
	else if (linkage::defaultResultType == linkage::PLOD) return plod;
	else if (linkage::defaultResultType == linkage::LOG10L) return log10L;
	else exit(1);
}

void mlink_opt_theta(std::istream& datafile, std::istream& ipedfile, std::istream& speedfile, double& theta, double& result)
{
	try
	{
		double log10L;
		mlink_common(datafile, ipedfile, speedfile, log10L);
		double rs;
		nlopt::opt oo(nlopt::LN_NELDERMEAD, 1); // LN_NELDERMEAD > LN_SBPLX > LN_BOBYQA >> LN_COBYLA LN_AUGLAG LN_PRAXIS GN_ISRES
		std::vector<double> x(1,0);
		std::vector<double> g;
		oo.set_max_objective(_mlink_opt_theta, NULL);
		oo.set_lower_bounds(0);
		oo.set_upper_bounds(0.5);
		nlopt::result nlopt_result=oo.optimize(x, rs);
		if (nlopt_result<0) { std::cerr<<"nlopt failed with code "<<nlopt_result<<std::endl; exit(1); }
		theta=x[0];
		result=rs;
	}
	catch (int e)
	{
		// There is an error/warning code e.
		theta = std::numeric_limits<double>::signaling_NaN();
		result = std::numeric_limits<double>::signaling_NaN();
	}
}

double _mlink_opt_scale(const std::vector<double> &x, std::vector<double> &grad, void* f_data)
{
	double log10L, lod, hlod, alpha, plod;
	ch_pen_dom( x[0] );
	holdtheta=maletheta->theta[which];
	maletheta->theta[which]=0.5;
	iterpeds(log10L, lod, hlod, alpha, plod);
	maletheta->theta[which]=holdtheta;
	iterpeds(log10L, lod, hlod, alpha, plod);
	if		(linkage::defaultResultType == linkage::LOD)  return lod;
	else if (linkage::defaultResultType == linkage::HLOD) return hlod;
	else if (linkage::defaultResultType == linkage::PLOD) return plod;
	else if (linkage::defaultResultType == linkage::LOG10L) return log10L;
	else exit(1);
}

void mlink_opt_scale(std::istream& datafile, std::istream& ipedfile, std::istream& speedfile, double& scale, double& result)
{
	try
	{
		double log10L;
		mlink_common(datafile, ipedfile, speedfile, log10L);
		double rs = std::numeric_limits<double>::signaling_NaN();
		double lb = 0.0001 / minpen;
		double ub = 0.9999 / maxpen;
		std::vector<double> x(1,1);
		
		//	method 1 = step-wise searching, good to find the global maximum but hard to pinpoint even when min_step=0.0001
		//	method 2 = Nelder Mead, may fall into local maximum.
		//	method 3 = combine these two: first step-wise search to find the region, then use Nelder Mead. May fall into local because step-wide is not global search
		//	method 4 = multiple nelder mead on equally divided regions. Never successfully ran. NLopt has bugs!
		
		// method 1
		rs = stepwise_opt_max(_mlink_opt_scale,x,lb,ub,0.001);
		
		// method 2 or 3
		nlopt::opt oo(nlopt::LN_NELDERMEAD, 1); // LN_NELDERMEAD > LN_SBPLX > LN_BOBYQA >> LN_COBYLA LN_AUGLAG LN_PRAXIS GN_ISRES
		oo.set_max_objective(_mlink_opt_scale, NULL);
		oo.set_lower_bounds(lb);
		oo.set_upper_bounds(ub);
		oo.optimize(x, rs);
		
		/*/ method 4
		 double len = (ub-lb)/5;
		 double max_y=-DBL_MAX;
		 double max_x=-DBL_MAX;
		 for (double bgn=lb; bgn<ub; bgn+=len)
		 {
		 double y;
		 double end=std::min(bgn+len,ub);
		 oo.set_lower_bounds(bgn);
		 oo.set_upper_bounds(end);
		 oo.optimize(x, y);
		 if (y>max_y) { max_y=y; max_x=x[0]; }
		 }
		 rs=max_y;
		 x[0]=max_x; //*/
		
		// report results;
		//	theta=maletheta->theta[which];
		scale=x[0];
		result=rs;
	}
	catch (int e)
	{
		// There is an error/warning code e.
		theta = std::numeric_limits<double>::signaling_NaN();
		result = std::numeric_limits<double>::signaling_NaN();
	}
}

double _mlink_opt_domin(const std::vector<double> &x, std::vector<double> &grad, void* f_data)
{
	double log10L, lod, hlod, alpha, plod;
	ch_pen_dom( x[0],x[1] );
	holdtheta=maletheta->theta[which];
	maletheta->theta[which]=0.5;
	iterpeds(log10L, lod, hlod, alpha, plod);
	maletheta->theta[which]=holdtheta;
	iterpeds(log10L, lod, hlod, alpha, plod);
	if		(linkage::defaultResultType == linkage::LOD)  return lod;
	else if (linkage::defaultResultType == linkage::HLOD) return hlod;
	else if (linkage::defaultResultType == linkage::PLOD) return plod;
	else if (linkage::defaultResultType == linkage::LOG10L) return log10L;
	else exit(1);
}

void mlink_opt_domin(std::istream& datafile, std::istream& ipedfile, std::istream& speedfile, double& phenocopy, double& penetrance, double& result)
{
	try
	{
		double log10L;
		mlink_common(datafile, ipedfile, speedfile, log10L);
		double rs = std::numeric_limits<double>::signaling_NaN();
		std::vector<double> x  = { 0.02,   0.2 };
		std::vector<double> lb = { 0.0001, 0.0001 };
		std::vector<double> ub = { 0.9999, 0.9999 };
		
		nlopt::opt oo(nlopt::LN_NELDERMEAD, 2); // LN_NELDERMEAD > LN_SBPLX > LN_BOBYQA >> LN_COBYLA LN_AUGLAG LN_PRAXIS GN_ISRES
		oo.set_max_objective(_mlink_opt_domin, NULL);
		oo.set_lower_bounds(lb);
		oo.set_upper_bounds(ub);
		oo.optimize(x, rs);
		
		phenocopy=x[0];
		penetrance=x[1];
		result=rs;
	}
	catch (int e)
	{
		// There is an error/warning code e.
		phenocopy = std::numeric_limits<double>::signaling_NaN();
		penetrance = std::numeric_limits<double>::signaling_NaN();
		result = std::numeric_limits<double>::signaling_NaN();
	}
}

#endif

#ifdef _MAKE_PROG

int main(int argc, const char * argv[])
{
	// pio_initialize(argc, argv); // removed since I don't have "output <<" anymore.
	// example data in /Users/fbj/work/projects/Ps/6.Eseq/2.linkage/try/
	std::fstream datafile("datafile.dat",std::ios_base::in);
	std::fstream ipedfile("ipedfile.dat",std::ios_base::in);
	std::fstream speedfile("speedfile.dat",std::ios_base::in);
	
	double theta, alpha, lod;
	mlink_calculate(datafile, ipedfile, speedfile, theta, alpha, lod);
	
	datafile.close();
	ipedfile.close();
	speedfile.close();
	return EXIT_SUCCESS;
}

#endif

int mlink_main()
{
	// pio_initialize(argc, argv); // removed since I don't have "output <<" anymore.
	// example data in /Users/fbj/work/projects/Ps/6.Eseq/2.linkage/try/
	std::fstream datafile("datafile.dat",std::ios_base::in);
	std::fstream ipedfile("ipedfile.dat",std::ios_base::in);
	std::fstream speedfile("speedfile.dat",std::ios_base::in);
	
	double theta, alpha, lod;
	mlink_calculate(datafile, ipedfile, speedfile, theta, alpha, lod);
	
	datafile.close();
	ipedfile.close();
	speedfile.close();
	return EXIT_SUCCESS;
}
