#include <ptoc/ptoc.h>
#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

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

	/*Unix Pascal*/
	/*Changes at Columbia U. indicated by "change"*/
	/*7 Feb 1995*/
	
	const char version[] = "5.23";	/*PRESENT VERSION OF LINKAGE*/
	const int maxlocus = 20;        /*MAXIMUM NUMBER OF LOCI*/
	const int maxall = 100;         /*MAX NUMBER OF ALLELES AT A SINGLE LOCUS 15*/
	const int maxgeno = maxall*(maxall+1) / 2; /*MAX NUMBER OF SINGLE LOCUS GENOTYPES*/
	/*Different definition than in analysis programs!*/
	const int maxind = 500000;		/*MAXIMUM NUMBER OF INDIVIDUALS IN A PEDIGREE previously 500 */
	const int maxmarriage = 3;		/*MAXIMUM NUMBER OF MARRIAGES FOR ONE MALE*/
	const int maxfact = maxall+2;	/*MAXIMUM NUMBER OF BINARY CODES AT A SINGLE LOCUS*/
	const int affall = 2;			/*DISEASE ALLELE FOR QUANT. TRAITS OR AFFECTION STATUS*/

	/* QUANTITATIVE TRAIT */
	const int maxtrait = 3;			/*MAXIMUM NUMBER OF QUANTITATIVE FACTORS AT A SINGLE LOCUS*/
	const real missval = 0.0;		/*MISSING VALUES FOR QUANTITATIVE TRAITS*/

	/* AFFECTION STATUS */
	const int missaff = 0;			/*MISSING VALUE FOR AFFECTION STATUS*/
	const int affval = 2;           /*CODE FOR AFFECTED INDIVIDUAL*/
	const int maxliab = 120;        /*MAXIMUM NUMBER OF LIABILITY CLASSES*/
	
	typedef array1d<1,maxgeno,boolean> genotype;
	enum direction {peelup,peeldown, last_direction};
	typedef set binset;
	typedef array1d<1,maxall,binset> phenarray;
	typedef array2d<1,maxall,1,maxall,boolean> possvect;
	typedef array1d<1,maxlocus,possvect> possarray;
	typedef struct locusvalues* locuspoint;
	typedef struct phenotype* phenpoint;
	enum locustype {affection,quantitative,binary, last_locustype};
	
	bool dooutput=false;
	
	struct locusvalues {
		int nallele;
		locustype which;
		union {
			int ntrait;
			struct {array4d<0,maxall,1,maxall,0,2,1,maxliab,real> pen;int nclass;} s_affection;
			struct {phenarray allele;int nfactor,format;} s_binary;
		};
	};
	
	struct phenotype {
		locustype which;
		union {
			struct {array1d<1,maxtrait,real> x;boolean missing;} s_quantitative;
			struct {int aff,liability;} s_affection;
			binset phenf;
		};
	};
	typedef struct thisperson* ind;
	typedef struct thisarray* genpoint;
	typedef array1d<1,maxlocus,phenpoint> indphen;
	
	struct thisarray {
		genotype genarray;
	};
	typedef struct information* infoptr;
	
	struct information {
		possarray possible;
	};
	
	struct thisperson {
		int id,paid,maid,offid,npaid,nmaid,sex,profield,oldped,nseq;
		ind pa,ma,foff,nextpa,nextma;
		genpoint gen;
		indphen phen;
		infoptr store_;
		array1d<1,maxlocus,boolean> thisunknown;
		boolean unknown,multi,done,up,male;
		boolean downvisit; // added 2015-05-17
	};
	enum haplotype {a,b, last_haplotype};
	typedef array1d<a,b,int> subhap;
	
	array1d<1,maxgeno,subhap> seghap;
	array1d<1,maxlocus,locuspoint> thislocus;
	array1d<0,maxind,ind> person;
	ind proband,loop1,loop2;
	array2d<1,maxgeno,1,maxgeno,int> genenumber;
	int risksys,mutsys,nsystem;
	unsigned fgeno,mgeno;
	int nsequence,newped,whichsys,totperson;
	boolean sexlink,risk,disequi;
	genotype gene;
	real one=1.00001; /*changed*/
	boolean makehomozygous;  /* Change - Added 7/8/93 */
	// below added on 2015-05-17
	int numind;  /*number of individuals in pedigree, for loop detection*/
	int depth;  /*depth of recursion*/
	const int DEPTH_MULTIPLE=3;
	const boolean DOWN_CHECK=true;
	
	void respond()
	/*Include file to LINKAGE programs*/
	{
		std::cerr << "  Press Enter key to continue or Ctrl-C to abort" << std::endl;
		input >> NL;
	}
	
	
	void inputerror(int nerror,int par1,int par2)
	{
		std::cerr << "Fatal error detected in procedure inputdata" << std::endl;
		switch (nerror) {
			case 0:std::cerr << "Number of loci " << fmt(par1,2) << " exceeds the constant maxlocus" << std::endl; break;
			case 1:std::cerr << "Number of loci read " << fmt(par1,2) << ". Less than minimum of 1" << std::endl; break;
			case 2:std::cerr << "Error detected reading loci order. Locus number " << fmt(par2,2) << " in position " << fmt(par1,2) << " exceeds number of loci" << std::endl; break;
			case 3:std::cerr << "Error detected reading loci order. Illegal locus number " << fmt(par2,2) << " in position " << fmt(par1,2) << std::endl; break;
			case 4:std::cerr << "Error detected reading loci order. Locus number repeated in positions " << fmt(par1,2) << " and " << fmt(par2,2) << std::endl; break;
			case 5:std::cerr << "Error detected reading locus description. Illegal locus type " << fmt(par2,2) << " for locus " << fmt(par1,2) << std::endl; break;
			case 6:std::cerr << "Error detected reading locus description for system " << fmt(par1,2) << ". Number of alleles  " << fmt(par1,2) << " exceeds maxall" << std::endl; break;
			case 7:std::cerr << "Error detected reading locus description for system " << fmt(par1,2) << ". Illegal number of alleles  " << fmt(par2,2) << std::endl; break;
			case 8:std::cerr << "Error detected reading locus description for system " << fmt(par1,2) << ". Number of factors  " << fmt(par2,2) << " exceeds maxfact" << std::endl; break;
			case 9:std::cerr << "Error detected reading locus description for system " << fmt(par1,2) << ". Illegal number of factors  " << fmt(par2,2) << std::endl; break;
			case 10:std::cerr << "Error detected reading locus description for system " << fmt(par1,2) << ". Alleles not codominant" << std::endl; break;
			case 11:std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Illegal code for sex " << fmt(par2,2) << std::endl; break;
			case 12:std::cerr << "Error detected reading pedigree record at pedigree" << fmt(par1,2) << ". Maximum number of pedigree records exceeded" << std::endl; break;
			case 13:std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Maximum number of individuals exceeded" << std::endl; break;
			case 14:std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Illegal binary factor code " << fmt(par2,2) << std::endl; break;
			case 15:std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". No allelic pair for genotype" << std::endl; break;
			case 16:std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Allele number " << fmt(par2,2) << " exceeds maxall" << std::endl; break;
			case 17:std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Illegal allele number " << fmt(par2,2) << std::endl; break;
			case 18:std::cerr << "Number of systems after factorization (" << fmt(par1,3) << ") exceeds maxsystem" << std::endl; break;
			case 19:std::cerr << "Number of systems after factorization (" << fmt(par1,3) << ") less than minimum of 1" << std::endl; break;
			case 20:std::cerr << "Number of recombination types (" << fmt(par1,3) << ") exceeds maxrectype" << std::endl; break;
			case 21:std::cerr << "Number of recombination types (" << fmt(par1,3) << ") less than minimum of 1" << std::endl; break;
			case 22:std::cerr << "End of file detected in tempdat by procedure readthg before all data found" << std::endl; break;
			case 23:std::cerr << "Error detected reading iterated locus in datafile. Value (" << fmt(par1,3) << ") greater than nlocus" << std::endl; break;
			case 24:std::cerr << "Error detected reading iterated locus in datafile. Illegal value (" << fmt(par1,3) << ')' << std::endl; break;
			case 25:std::cerr << "Number of iterated parameters greater then maxn" << std::endl; break;
			case 26:std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Liability class (" << fmt(par2,2) << ") exceeds nclass" << std::endl; break;
			case 27:std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". Illegal liability class (" << fmt(par2,2) << ')' << std::endl; break;
			case 28:std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Liability classes (" << fmt(par2,3) << ") exceed maxliab" << std::endl; break;
			case 29:std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Illegal number of liability classes (" << fmt(par2,3) << ')' << std::endl; break;
			case 30:std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Penetrance out of range" << std::endl; break;
			case 31:std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Number of traits (" << fmt(par2,3) << ") exceeds maxtrait" << std::endl; break;
			case 32:std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Number of traits out of range (" << fmt(par2,3) << ')' << std::endl; break;
			case 33:std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Variance must be positive" << std::endl; break;
			case 34:std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Variance multiplier must be positive" << std::endl; break;
			case 35:std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Risk allele " << fmt(par2,3) << ") exceeds nallele" << std::endl; break;
			case 36:std::cerr << "Error detected reading locus description for system" << fmt(par1,2) << ". Illegal risk allele (" << fmt(par2,3) << ')' << std::endl; break;
			case 37:std::cerr << "Error detected reading datafile. Risk locus " << fmt(par2,3) << ") exceeds nlocus" << std::endl; break;
			case 38:std::cerr << "Error detected reading datafile. Illegal value for risk locus " << fmt(par2,3) << ')' << std::endl; break;
			case 39:std::cerr << "Error detected reading datafile. Mutation locus " << fmt(par2,3) << ") exceeds nlocus" << std::endl; break;
			case 40:std::cerr << "Error detected reading datafile. Illegal value for mutation locus " << fmt(par2,3) << ')' << std::endl; break;
			case 41:std::cerr << "Error detected reading datafile. Linkage disequilibrium is not allowed with this program" << std::endl; break;
			case 42:std::cerr << "Locus " << fmt(par1,5) << " in lod score list exceeds nlocus " << fmt(par2,5) << std::endl; break;
			case 43:std::cerr << "Illegal locus number " << fmt(par1,5) << " in lod score list" << std::endl; break;
			case 44:std::cerr << "Error detected reading pedigree record " << fmt(par1,2) << ". One 0 allele" << std::endl; break;
		}
		exit(1);
		respond(); /*changed*/
	}
	
	
	void inputwarning(int nwarning,int par1,int par2)
	{
		std::cerr << "Warning number from procedure inputdata" << std::endl;
		switch (nwarning) {
			case 0:std::cerr << "Illegal sex difference parameter " << fmt(par1,2) << " Parameter should be 0, 1, or 2" << std::endl; break;
			case 1:std::cerr << "Illegal interference parameter " << fmt(par1,2) << " Lack of interference assumed" << std::endl; break;
			case 2:std::cerr << "Illegal sex difference parameter " << fmt(par1,2) << " Parameter must be 0 with sex-linked data" << std::endl; break;
			case 3:std::cerr << "Non-standard affection status" << fmt(par2,4) << /*" interpreted as normal*/" in pedigree record" << fmt(par1,5) << std::endl; break;
		}
		exit(1);
		respond(); /*changed*/
	}
	
	void writespeed(std::ostream& speedfile)
	{
		int i,j,a,b;
		
		/*writespeed*/
		for( i=1; i <= totperson; i ++)
		{ thisperson& with = *person[i];
			if (with.unknown && (with.foff!=nil))
			{
				information& with1 = *with.store_;
				
				speedfile << "id" << fmt(with.nseq,7) << std::endl;
				for( j=1; j <= nsystem; j ++)
					for( a=1; a <= thislocus[j]->nallele; a ++)
						for( b=1; b <= thislocus[j]->nallele; b ++)
							if (with1.possible[j][a][b])
								speedfile << fmt(j,3) << fmt(a,3) << fmt(b,3) << std::endl;
			}}
	}      /*writespeed*/
	
	
	void writeped(std::ostream& ipedfile)
	{
		int i,j,k,a,b;
		
		/*writeped*/
		for( i=1; i <= totperson; i ++)
		{
			thisperson& with = *person[i];
			
			ipedfile << fmt(with.oldped,6) << ' ' << fmt(with.id,4) << ' ' << fmt(with.paid,4) << ' ' << fmt(with.maid,4) << ' ' << fmt(with.offid,4) << ' ' << fmt(with.npaid,4) << ' ';
			ipedfile << fmt(with.nmaid,4) << ' ' << fmt(with.sex,1) << ' ' << fmt(with.profield,1) << ' ';
			for( j=1; j <= nsystem; j ++)
			{
				{
					phenotype& with1 = *with.phen[j];
					{
						locusvalues& with2 = *thislocus[j];
						if (with2.which==binary)
							if (with2.s_binary.format==2)
								for( k=1; k <= with2.s_binary.nfactor; k ++)
									if (with1.phenf.has(k))
										ipedfile << " 1";
									else ipedfile << " 0";
									else
									{
										a=0;
										b=0;
										for( k=1; k <= with2.nallele; k ++)
											if (with1.phenf.has(k))
											{
												if (a==0) a=k;
												else b=k;
											}
										if (b==0)  b=a;
										ipedfile << fmt(a,3) << fmt(b,3);
									}
									else if (with2.which==quantitative)
										if ((! sexlink) || (! with.male))
											for( k=1; k <= with2.ntrait; k ++)
												ipedfile << ' ' << fmt(with1.s_quantitative.x[k],9,4);
										else for( k=1; k <= with2.ntrait; k ++)
											ipedfile << ' ' << fmt(with1.s_affection.aff,9);
										else
										{
											ipedfile << fmt(with1.s_affection.aff,2);
											if (with2.s_affection.nclass!=1)
												ipedfile << fmt(with1.s_affection.liability,3);
										}
					}
				}
				if (j!=nsystem)
					ipedfile << ' ';
			}
			ipedfile << std::endl;
		}
	}      /*writeped*/
	
	
	void infer()
	{
		int i,j,k,l,kposs,lposs,count,pacount,macount;
		boolean someknown;
		
		/*infer*/
		for( i=1; i <= totperson; i ++)
			if (person[i]->unknown)
			{ thisperson& with = *person[i];
				{
					information& with1 = *with.store_;
					
					for( j=1; j <= nsystem; j ++)
						if (thislocus[j]->which==binary)
							if (with.phen[j]->phenf==set::of(eos))
							{
								locusvalues& with2 = *thislocus[j];
								
								count=0;
								for( k=1; k <= with2.nallele; k ++)
									for( l=k; l <= with2.nallele; l ++)
										if (with1.possible[j][k][l])
										{
											kposs=k;
											lposs=l;
											count=count+1;
										}
								if (count==1)
								{
									if (sexlink && with.male)
										with.phen[j]->phenf=with2.s_binary.allele[lposs];
									else with.phen[j]->phenf=with2.s_binary.allele[kposs]+with2.s_binary.allele[lposs];
								}
							}
					count=0;
					for( j=1; j <= nsystem; j ++)
						if (thislocus[j]->which!=binary)
							count=count+1;
						else if (with.phen[j]->phenf==set::of(eos))
							count=count+1;
					with.unknown=count!=0;
				}}
		/*Infer children when parents are homozygotes*/
		for( i=1; i <= totperson; i ++)
			if (person[i]->foff==nil)
			{ thisperson& with = *person[i];
				for( j=1; j <= nsystem; j ++)
				{ locusvalues& with1 = *thislocus[j];
					if (with.phen[j]->which==binary)
						if (with.phen[j]->phenf==set::of(eos))
							if (with.pa!=nil)
							{
								pacount=0;
								macount=0;
								for( k=1; k <= thislocus[j]->nallele; k ++)
									if (with1.s_binary.allele[k]<=with.pa->phen[j]->phenf)
									{
										kposs=k;
										pacount=pacount+1;
									}
								for( l=1; l <= thislocus[j]->nallele; l ++)
									if (with1.s_binary.allele[l]<=with.ma->phen[j]->phenf)
									{
										lposs=l;
										macount=macount+1;
									}
								if (((macount==1) && (pacount==1)) && !(with.male && sexlink))
								{
									with.phen[j]->phenf=with1.s_binary.allele[kposs]+with1.s_binary.allele[lposs];
								}
								else if ((macount==1) && (with.male && sexlink))
								{
									with.phen[j]->phenf=with1.s_binary.allele[lposs];
								}
							}}}
		/*Replace by homozygotes if all unknown in a pedigree*/
		if (makehomozygous)    /*change - added*/
			for( j=1; j <= nsystem; j ++)
			{ locusvalues& with = *thislocus[j];
				if ((with.which==binary) && (with.s_binary.format==3)) /*change - 'format=3' added*/
				{
					someknown=false;
					for( i=1; i <= totperson; i ++)
						if (person[i]->phen[j]->phenf!=set::of(eos))
							someknown=true;
					if (! someknown)
						for( i=1; i <= totperson; i ++)
							person[i]->phen[j]->phenf=with.s_binary.allele[1];
				}}
	}      /*infer*/
	
	
	void getunknown()
	{
		int i,j,n,ahap,bhap;
		
		/*getunknown*/
		for( i=1; i <= totperson; i ++)
			person[i]->unknown=false;
		for( i=1; i <= totperson; i ++)
			for( j=1; j <= nsystem; j ++)
				person[i]->thisunknown[j]=false;
		for( i=1; i <= totperson; i ++)
		{
			thisperson& with = *person[i];
			for( j=1; j <= nsystem; j ++)
			{
				if (thislocus[j]->which==binary)
				{
					if (with.phen[j]->phenf==set::of(eos))
						with.thisunknown[j]=true;
					else if (thislocus[j]->which==quantitative)
					{	if (with.phen[j]->s_quantitative.x[1]==missval)
						with.thisunknown[j]=true;
					else if (with.phen[j]->s_affection.aff==missaff)
						with.thisunknown[j]=true;
					}
				}
				if (with.thisunknown[j])
					with.unknown=true;
			}
		}
		for( i=1; i <= totperson; i ++)
		{ thisperson& with = *person[i];
			if (with.unknown)
			{
				with.store_ = new information;
				{ information& with1 = *with.store_;
					for( n=1; n <= nsystem; n ++)
					{ locusvalues& with2 = *thislocus[n];
						for( ahap=1; ahap <= with2.nallele; ahap ++)
							for( bhap=1; bhap <= with2.nallele; bhap ++)
								with1.possible[n][ahap][bhap]=true;}}
			}
		}
	}      /*getunknown*/
	
	
	void getlocation(locuspoint thislocus)
	
	{
		int ahap,bhap,here;
		
		/*getlocation*/
		here=0;
		{
			locusvalues& with = *thislocus;
			for( ahap=1; ahap <= with.nallele; ahap ++)
				for( bhap=ahap; bhap <= with.nallele; bhap ++)
				{
					here=here+1;
					genenumber[ahap][bhap]=here;
					genenumber[bhap][ahap]=here;
					seghap[here][a]=ahap;
					seghap[here][b]=bhap;
				}
		}
	}      /*getlocation*/
	
	
	static void readbin(std::istream& pedfile, phenpoint& phen,locuspoint thislocus, ind& p)
	{
		int i,j;
		
		{
			locusvalues& with = *thislocus;
			{
				phenotype& with1 = *phen;
				
				with1.which=binary;
				with1.phenf=set::of(eos);
				for( i=1; i <= with.s_binary.nfactor; i ++)
				{
					pedfile >> j;
					if ((j!=0) && (j!=1))  inputerror(14,p->id,j);
					if (j==1)
						with1.phenf=with1.phenf+set::of(i, eos);
				}
			}
		}
	}
	
	
	static void readnumber(std::istream& pedfile, phenpoint& phen, ind& p)
	{
		int j,k;
		
		{
			phenotype& with = *phen;
			
			with.which=binary;
			with.phenf=set::of(eos);
			pedfile >> j >> k;
			if (j>maxall)
				inputerror(16,p->id,j);
			if (j<0)
				inputerror(17,p->id,j);
			if (k>maxall)
				inputerror(16,p->id,k);
			if (k<0)
				inputerror(17,p->id,k);
			if (((j==0) || (k==0)) && (j!=k))
				inputerror(44,p->id,j);
			else
			{
				if (j!=0)  with.phenf=with.phenf+set::of(j, eos);
				if (k!=0)  with.phenf=with.phenf+set::of(k, eos);
			}
		}
	}
	
	
	
	static void readaff(std::istream& pedfile, phenpoint& phen,locuspoint thislocus, ind& p)
	{
		int thisval;
		
		
		{
			phenotype& with = *phen;
			
			with.which=affection;
			pedfile >> thisval;
			if (thisval==missaff)
				
				with.s_affection.aff=0;
			else
				if (thisval==affval)
					with.s_affection.aff=2;
				else
				{
					if (thisval!=1)  inputwarning(3,p->id,thisval);
					with.s_affection.aff=1;
				}
			if (thislocus->s_affection.nclass==1)
				with.s_affection.liability=1;
			else pedfile >> with.s_affection.liability;
			if (with.s_affection.liability>thislocus->s_affection.nclass)  inputerror(26,p->id,with.s_affection.liability);
			if (with.s_affection.liability<=0)  inputerror(27,p->id,with.s_affection.liability);
		}
	}
	
	
	
	static void readquan(std::istream& pedfile, phenpoint& phen,locuspoint thislocus, ind& p)
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
					
					for( i=1; i <= with1.ntrait; i ++) pedfile >> with.s_quantitative.x[i];
					with.s_quantitative.missing=true;
					for( i=1; i <= with1.ntrait; i ++)
						if (with.s_quantitative.x[i]!=missval)
							with.s_quantitative.missing=false;
				}
			}
			else
			{
				with.which=affection;
				pedfile >> xval;
				{
					locusvalues& with1 = *thislocus;
					
					if (xval==missval)
						
						with.s_affection.aff=missaff;
					else
						if (xval==affall)
							with.s_affection.aff=affall;
						else with.s_affection.aff=-11;
					with.s_affection.liability=1;
					for( i=2; i <= with1.ntrait; i ++) pedfile >> xval;
				}
			}
		}
	}
	
	
	
	static void getphenotype(std::istream& pedfile, ind& p)
	{
		int thisread,system;
		{
			thisperson& with = *p;
			for( thisread=1; thisread <= nsystem; thisread ++)
			{
				system=thisread;
				with.phen[system]=nil;
				with.phen[system] = new phenotype;
				switch (thislocus[system]->which) {
					case quantitative:readquan(pedfile,with.phen[system],thislocus[system], p); break;
					case affection:readaff(pedfile,with.phen[system],thislocus[system], p); break;
					case binary:
						if (thislocus[system]->s_binary.format==3) readnumber(pedfile,with.phen[system], p);
						else readbin(pedfile,with.phen[system],thislocus[system], p);
						break;
					case last_locustype: exit_error("thislocus[system]->which = last_locustype (fbj)");
				}
			}
		}
	}
	
	
	static void getind(std::istream& pedfile, int& id,int& seq)
	{     /*getind*/
		id=0;
		pedfile >> seq;
		if (seq!=0)
		{
			id=seq;
			if (id>maxind)
				inputerror(13,id,id);
			if (person[id]==nil)
			{
				numind++;
				person[id] = new thisperson;
				{
					thisperson& with = *person[id];
					with.gen = new thisarray;
					with.nseq=seq+nsequence;
				}
			}
		}
	}       /*getind*/
	
	
	static void multimarriage(ind& p)
	{
		ind q,child;
		
		/*multimarriage*/
		if (p->foff!=nil)
		{
			thisperson& with = *p;
			
			if (with.male)
				q=with.foff->ma;
			else q=with.foff->pa;
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
	}       /*multimarriage*/
	
	
	
	
	static void getrest(std::istream& pedfile)
	{
		char whichchr;
		
		/*getrest*/
		whichchr=' ';
		while (!(is_EndOfLine(pedfile.peek()) || (whichchr==chr(12))))
			pedfile >> whichchr;
		pedfile >> GEOL;
	}       /*getrest*/
	
	void readped(std::istream& pedfile)
	{
		numind=0;
		depth=0;
		int i,newid,thisone,thisped,tempid;
		
		for( i=0; i <= maxind;
			i ++) person[i]=nil;
		totperson=0;
		loop1=nil;
		loop2=nil;
		proband=nil;
		thisped=newped;
		if (dooutput) std::cerr << "Ped. " << fmt(thisped,2) << std::endl; /*change - added*/
		while (! is_blank_row(pedfile) && (thisped==newped))
		{
			totperson=totperson+1;
			getind(pedfile,thisone,tempid);
			if (proband==nil)
				proband=person[thisone];
			{
				thisperson& with = *person[thisone];
				
				with.id=tempid;
				with.oldped=newped;
				getind(pedfile,newid,with.paid);
				with.pa=person[newid];
				getind(pedfile,newid,with.maid);
				with.ma=person[newid];
				getind(pedfile,newid,with.offid);
				with.foff=person[newid];
				getind(pedfile,newid,with.npaid);
				with.nextpa=person[newid];
				getind(pedfile,newid,with.nmaid);
				with.nextma=person[newid];
				with.store_=nil;
				with.up=false;
				pedfile >> with.sex;
				if ((with.sex!=1) && (with.sex!=2))
					inputerror(11,with.id,with.sex);
				with.male=(with.sex==1);
				pedfile >> with.profield;
				if (with.profield==1)
					proband=person[thisone];
				else if (with.profield==2)
				{	if (loop2==nil)
					loop2=person[thisone];
				else loop1=person[thisone];
				}
			}
			getphenotype(pedfile,person[thisone]);
			getrest(pedfile);
			if (! is_blank_row(pedfile))
				pedfile >> newped;
		}
		nsequence=totperson+nsequence;
		if ((loop2!=nil) && (loop1==nil))
			loop1=proband;
		for( thisone=1; thisone <= totperson;
			thisone ++) multimarriage(person[thisone]);
	}
	
	
	
	static void getquan(std::istream& datafile, locuspoint& locus, int& system)
	{
		int i;
		
		/*getquan*/
		{
			locusvalues& with = *locus;
			
			datafile >> with.ntrait >> GEOL;
			if (with.ntrait>maxtrait)  inputerror(31,system,with.ntrait);
			if (with.ntrait<=0)  inputerror(32,system,with.s_affection.nclass);
			for( i=1; i <= with.ntrait; i ++) /*changed*/
				datafile >> GEOL;   /*skip genotype means*/
			for( i=1; i <= 2; i ++) datafile >> GEOL;
		}
	}       /*getquan*/
	
	
	static void getpen(std::istream& datafile, locuspoint& locus, int& system)
	{
		int i,j,k,l;
		
		/*getpen*/
		{
			locusvalues& with = *locus;
			
			datafile >> with.s_affection.nclass >> GEOL;
			if (with.s_affection.nclass>maxliab)  inputerror(28,system,with.s_affection.nclass);
			if (with.s_affection.nclass<=0)  inputerror(29,system,with.s_affection.nclass);
			for( l=1; l <= with.s_affection.nclass; l ++)
			{
				for( i=1; i <= with.nallele; i ++)
					for( j=i; j <= with.nallele; j ++)
					{
						datafile >> with.s_affection.pen[i][j][2][l];
						if ((with.s_affection.pen[i][j][2][l]<0) || (with.s_affection.pen[i][j][2][l]>one))  inputerror(30,system,system);
						with.s_affection.pen[i][j][1][l]=1-with.s_affection.pen[i][j][2][l];
						with.s_affection.pen[i][j][0][l]=1.0;
						for( k=0; k <= 2; k ++)
							with.s_affection.pen[j][i][k][l]=with.s_affection.pen[i][j][k][l];
					}
				datafile >> GEOL;
				for( i=1; i <= with.nallele; i ++)
					with.s_affection.pen[0][i][0][l]=1.0;
				if (sexlink)
				{
					for( i=1; i <= with.nallele; i ++)
						datafile >> with.s_affection.pen[0][i][2][l];
					if ((with.s_affection.pen[0][j][2][l]<0) || (with.s_affection.pen[0][j][2][l]>one))  inputerror(30,system,system);
					for( i=1; i <= with.nallele; i ++)
						with.s_affection.pen[0][i][1][l]=1.0-with.s_affection.pen[0][i][2][l];
					datafile >> GEOL;
				}
			}
		}
	}       /*getpen*/
	
	
	static void getbin(std::istream& datafile, locuspoint& locus, int& system)
	{
		int i,j,k;
		
		/*getbin*/
		{
			locusvalues& with = *locus;
			
			datafile >> with.s_binary.nfactor >> GEOL;
			if (with.s_binary.nfactor>maxfact)  inputerror(8,system,with.s_binary.nfactor);
			if (with.s_binary.nfactor<=0)  inputerror(9,system,with.s_binary.nfactor);
			for( i=1; i <= with.nallele; i ++)
				with.s_binary.allele[i]=set::of(eos);
			for( i=1; i <= with.nallele; i ++)
				for( j=1; j <= with.s_binary.nfactor; j ++)
				{
					datafile >> k;
					if (k==1)
						with.s_binary.allele[i]=with.s_binary.allele[i]+set::of(j, eos);
				}
		}
		datafile >> GEOL;
	}       /*getbin*/
	
	
	static void getnumber(locuspoint& locus)
	{
		int i;
		
		/*getnumber*/
		{ locusvalues& with = *locus;
			for( i=1; i <= with.nallele; i ++)
				with.s_binary.allele[i]=set::of(i, eos);}
	}       /*getnumber*/
	
	
	static void getlocus(std::istream& datafile, int system, int& whichtype)
	{     /*getlocus*/
		thislocus[system] = new locusvalues;
		{
			locusvalues& with = *thislocus[system];
			
			datafile >> whichtype >> with.nallele;
			switch (whichtype) {
				case 0:with.which=quantitative; break;
				case 1:with.which=affection; break;
				case 2:{
					with.which=binary;
					with.s_binary.format=2;
				}
					break;
				case 3:{
					with.which=binary;
					with.s_binary.format=3;
				}
					break;
			}
			datafile >> GEOL;
			if (! disequi)
				datafile >> GEOL;
			switch (with.which) {
				case quantitative:getquan(datafile,thislocus[system], system); break;
				case affection:getpen(datafile,thislocus[system], system); break;
				case binary:
					if (with.s_binary.format==2) getbin(datafile,thislocus[system], system);
					else getnumber(thislocus[system]);
					break;
				case last_locustype: exit_error("thislocus[system]->which = last_locustype (fbj)");
			}
			if (risk && (system==risksys))
				datafile >> GEOL;
		}
	}       /*getlocus*/
	
	void readloci(std::istream& datafile)
	{
		int i,coupling,autosomal,whichtype;
		real mu;
		
		/*readloci*/
		datafile >> nsystem >> risksys >> autosomal >> GEOL;
		if (nsystem>maxlocus)  inputerror(0,nsystem,nsystem);
		if (nsystem<=0)  inputerror(1,nsystem,nsystem);
		risk=(risksys!=0);
		sexlink=(autosomal==1);
		if (dooutput)
		{
			std::cout << "YOU ARE USING LINKAGE (V" << version << ") WITH" << fmt(nsystem,3) << "-POINT";
			if (sexlink) std::cout << " SEXLINKED DATA" << std::endl;
			else std::cout << " AUTOSOMAL DATA" << std::endl;
		}
		datafile >> mutsys >> mu >> mu >> coupling >> GEOL;
		disequi=(coupling==1);
		datafile >> GEOL;
		for( i=1; i <= nsystem; i ++) getlocus(datafile, i, whichtype);
	}      /*readloci*/
	
	
	void cleanup(ind& p)
	{
		unsigned i,j; // prv int
		
		/*cleanup*/
		{ thisperson& with = *p;
			if (with.unknown)
			{ information& with1 = *with.store_;
				{ thisarray& with2 = *with.gen;
					if (sexlink && with.male)
					{
						for( i=1; i <= mgeno; i ++)
							if (! with2.genarray[i])
								with1.possible[whichsys][1][i]=false;
						for( i=2; i <= mgeno; i ++)
							for( j=1; j <= mgeno; j ++)
								with1.possible[whichsys][i][j]=false;
					}
					else for( i=1; i <= fgeno; i ++)
					{
						if (! with2.genarray[i])
							with1.possible[whichsys][seghap[i][a]][seghap[i][b]]=false;
						with1.possible[whichsys][seghap[i][b]][seghap[i][a]]=with1.possible[whichsys][seghap[i][a]][seghap[i][b]];
					}}}}
	}      /*cleanup*/
	
	
	void getgene(int system,ind p,indphen phen)
	
	{
		int here,i,j;
		
		/*getgene*/
		here=0;
		{
			locusvalues& with = *thislocus[system];
			if (sexlink && p->male)
				for( i=1; i <= with.nallele; i ++)
				{
					here=here+1;
					switch (with.which) {
						case quantitative:{ phenotype& with1 = *phen[system];
							if (i==affall)
								p->gen->genarray[here]=(with1.s_affection.aff==affall) || (with1.s_affection.aff==missaff);
							else p->gen->genarray[here]=(with1.s_affection.aff!=affall) || (with1.s_affection.aff==missaff);}
							break;
						case affection:{ phenotype& with1 = *phen[system];  p->gen->genarray[here]=(with.s_affection.pen[0][i][with1.s_affection.aff][with1.s_affection.liability]>0.0);}
							break;
						case binary:{ phenotype& with1 = *phen[system];  p->gen->genarray[here]=(with1.phenf==with.s_binary.allele[i]) || (with1.phenf==set::of(eos));}
							break;
						case last_locustype: exit_error("thislocus[system]->which = last_locustype (fbj)");
					}
				}
			else for( i=1; i <= with.nallele; i ++)
				for( j=i; j <= with.nallele; j ++)
				{
					here=here+1;
					switch (with.which) {
						case quantitative:p->gen->genarray[here]=true; break;
						case affection:{ phenotype& with1 = *phen[system];  p->gen->genarray[here]=(with.s_affection.pen[i][j][with1.s_affection.aff][with1.s_affection.liability]>0.0);}
							break;
						case binary:{ phenotype& with1 = *phen[system];  p->gen->genarray[here]=(with1.phenf==(with.s_binary.allele[i]+with.s_binary.allele[j])) || (with1.phenf==set::of(eos));}
							break;
						case last_locustype: exit_error("thislocus[system]->which = last_locustype (fbj)");
					}
				}
		}
	}      /*getgene*/
	
	
	void likelihood(ind& proband);
	static void seg(ind& p,ind& q,ind& r,direction peel);
	
	
	
	static boolean segfun(ind& child,int first,int second, subhap& secondhap, subhap& firsthap, ind& p)
	
	{
		boolean temp;
		haplotype thishap,thathap;
		
		/*segfun*/
		boolean segfun_result;
		temp=false;
		{ thisarray& with = *child->gen;
			if (! sexlink)
				for( thishap=a; thishap <= b; thishap = succ(haplotype,thishap))
					for( thathap=a; thathap <= b; thathap = succ(haplotype,thathap))
						temp=temp || with.genarray[genenumber[secondhap[thishap]][firsthap[thathap]]];
			else if (child->male)
				if (p->male)
					for( thathap=a; thathap <= b; thathap = succ(haplotype,thathap))
						temp=temp || with.genarray[secondhap[thathap]];
				else for( thathap=a; thathap <= b; thathap = succ(haplotype,thathap))
					temp=temp || with.genarray[firsthap[thathap]];
				else if (p->male)
					for( thathap=a; thathap <= b; thathap = succ(haplotype,thathap))
						temp=temp || with.genarray[genenumber[secondhap[thathap]][first]];
				else for( thathap=a; thathap <= b; thathap = succ(haplotype,thathap))
					temp=temp || with.genarray[genenumber[firsthap[thathap]][second]];}
		segfun_result=temp;
		return segfun_result;
	}       /*segfun*/
	
	
	
	
	static void segup(ind& p, int& nfirst, int& nsecond, subhap& firsthap, ind& q, subhap& secondhap, ind& child, ind& father, ind& mother)
	
	{
		int first,second;
		boolean segval,val;
		
		/*segup*/
		if (p->male)
		{
			nfirst=mgeno;
			nsecond=fgeno;
		}
		else {
			nfirst=fgeno;
			nsecond=mgeno;
		}
		{ thisarray& with = *p->gen;
			for( first=1; first <= nfirst; first ++)
				if (with.genarray[first])
				{
					segval=false;
					firsthap=seghap[first];
					{ thisarray& with1 = *q->gen;
						for( second=1; second <= nsecond; second ++)
							if (with1.genarray[second])
							{
								secondhap=seghap[second];
								val=true;
								child=father->foff;
								do {
									if (mother==child->ma)
										val=segfun(child,first,second, secondhap, firsthap, p);
									child=child->nextpa;
								} while (!((! val) || (child==nil)));
								segval=val || segval;
							}}
					with.genarray[first]=segval;
				}}
		cleanup(q);
		child=father->foff;
		do {
			if (child->ma==mother)
				cleanup(child);
			child=child->nextpa;
		} while (!(child==nil));
	}       /*segup*/
	
	
	
	
	static void segdown(ind& p, subhap& firsthap, ind& q, subhap& secondhap, ind& child, ind& father, ind& mother, ind& r)
	
	{
		unsigned first,second; // prv int
		int here;
		boolean val;
		haplotype thishap,thathap;
		
		/*segdown*/
		for( first=1; first <= fgeno; first ++) gene[first]=false;
		{ thisarray& with = *p->gen;
			for( first=1; first <= mgeno; first ++)
				if (with.genarray[first])
				{
					firsthap=seghap[first];
					{ thisarray& with1 = *q->gen;
						for( second=1; second <= fgeno; second ++)
							if (with1.genarray[second])
							{
								secondhap=seghap[second];
								val=with1.genarray[second] && p->gen->genarray[first];
								child=father->foff;
								do {
									if (child->ma==mother)
										if (! child->up)
											val=segfun(child,first,second, secondhap, firsthap, p);
									child=child->nextpa;
								} while (!((! val) || (child==nil)));
								if (val)
								{
									if (! sexlink)
										for( thishap=a; thishap <= b; thishap = succ(haplotype,thishap))
											for( thathap=a; thathap <= b; thathap = succ(haplotype,thathap))
											{
												here=genenumber[secondhap[thishap]][firsthap[thathap]];
												gene[here]=gene[here] || val;
											}
									else if (r->male)
										for( thathap=a; thathap <= b; thathap = succ(haplotype,thathap))
										{
											here=secondhap[thathap];
											gene[here]=gene[here] || val;
										}
									else for( thathap=a; thathap <= b; thathap = succ(haplotype,thathap))
									{
										here=genenumber[secondhap[thathap]][first];
										gene[here]=gene[here] || val;
									}
								}
							}}
				}}
		{ thisarray& with = *r->gen;
			for( first=1; first <= fgeno; first ++)
				with.genarray[first]=with.genarray[first] && gene[first];}
		cleanup(p);
		cleanup(q);
		child=father->foff;
		do {
			if (child->ma==mother)
				if (! child->up)
					cleanup(child);
			child=child->nextpa;
		} while (!(child==nil));
	}       /*segdown*/
	
	
	
	
	static void seg(ind& p,ind& q,ind& r,direction peel)
	
	{
		ind child,father,mother;
		subhap firsthap,secondhap;
		int nfirst,nsecond;
		
		
		/*seg*/
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
			segup(p, nfirst, nsecond, firsthap, q, secondhap, child, father, mother);
		else segdown(p, firsthap, q, secondhap, child, father, mother, r);
	}       /*seg*/
	
	
	static void collapsedown(ind p);
	static void collapseup(ind p);
	
	static void collapseup(ind p)
	
	{
		ind q,child,nextchild;
		boolean down;
		
		if (++depth > (DEPTH_MULTIPLE*numind)) exit_error("There may be loops in one of the pedigrees.");

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
					collapsedown(q);
					nextchild=child;
					while (nextchild!=nil)
					{
						if ((nextchild->pa==q) || (nextchild->ma==q))
						{	if (! nextchild->up)
							collapseup(nextchild);
						else down=true;
						}
						if (p->male)
							nextchild=nextchild->nextpa;
						else nextchild=nextchild->nextma;
					}
					if (q->multi)
						collapseup(q);
					if (! down)
						seg(p,q,child,peelup);
					else collapsedown(p);
				}
				if (p->male)
					child=child->nextpa;
				else child=child->nextma;
			}
		}
		--depth;
	}       /*collapseup*/
	
	
	static void collapsedown(ind p)
	{
		if (++depth > (DEPTH_MULTIPLE * numind))
			exit_error("There may be loops in one of the pedigrees.");
		if (DOWN_CHECK && p->downvisit && (p != proband))
			exit_error("There may be loops in one of the pedigrees.");
		else
			p->downvisit = true;
		
		/*collapsedown*/
		if (p->pa!=nil)
		{
			p->up=true;
			collapseup(p->pa);
			seg(p->pa,p->ma,p,peeldown);
		}
		--depth;
	}       /*collapsedown*/
	
	void likelihood(ind& proband)
	
	{
		int i,j;
		
		
		void collapsedown(ind p);
		
		
		/*likelihood*/
		collapsedown(proband);
		collapseup(proband);
		if (proband->thisunknown[whichsys])
		{
			{ thisperson& with = *proband;
				{ information& with1 = *with.store_;
					{ thisarray& with2 = *with.gen;
						{ locusvalues& with3 = *thislocus[whichsys];
							if (sexlink && with.male)
								for( j=1; j <= with3.nallele; j ++) with1.possible[whichsys][1][j]=with2.genarray[j];
							else for( i=1; i <= with3.nallele; i ++)
								for( j=i; j <= with3.nallele; j ++)
								{
									with1.possible[whichsys][i][j]=with2.genarray[genenumber[i][j]];
									with1.possible[whichsys][j][i]=with1.possible[whichsys][i][j];
								}}}}}
			cleanup(proband);
		}
	}      /*likelihood*/
	
	int iterpeds_errors=0;
	
	void iterpeds()
	{
		int i,j;
		boolean compattest,compat;
		
		/*iterpeds*/
		if ((loop1==nil) && (loop2==nil))
		/* This means that this part of unknown is not active for pedigrees with loops! */
		{
			for( i=1; i <= totperson; i ++)
				getgene(whichsys,person[i],person[i]->phen);
			compattest=false;
			compat=false;
			for( i=1; i <= totperson; i ++)
				if ((! compattest) || (person[i]->thisunknown[whichsys] && compat))
				{
					for( j=1; j <= totperson; j ++) person[j]->done=false;
					for( j=1; j <= totperson; j ++) person[j]->up=false;
					for( j=1; j <= totperson; j ++) person[j]->downvisit=false;
					likelihood(person[i]);
					if (! compattest)
					{
						compattest=true;
						for( j=1; j <= (int)fgeno; j ++)
							compat=compat || person[i]->gen->genarray[j];
						if (! compat)  {
//							std::cerr << "ERROR: Incompatibility detected in this family for locus " << whichsys << std::endl;
//							respond();
							++iterpeds_errors;
						}
					}
				}
			for( i=1; i <= totperson; i ++)
				if (person[i]->unknown)
					cleanup(person[i]);
		}
	}      /*iterpeds*/
	
	
	void reinit()
	{
		int i,j;
		
		/*reinit*/
		for( i=1; i <= totperson; i ++)
			for( j=1; j <= nsystem; j ++)
				delete person[i]->phen[j];
		for( i=1; i <= totperson; i ++)
			if (person[i]->store_!=nil)
				delete person[i]->store_;
		for( i=1; i <= totperson; i ++)
		{
			delete person[i]->gen;
			delete person[i];
		}
		for( i=1; i <= totperson; i ++)
			person[i]=nil;
	}      /*reinit*/
	
	
	static boolean testhets(std::istream& datafile)  /*Change: Function added by Joe Terwilliger 7/8/93*/
	{
		real a,b,c;
		int prog,d,numl,lc,nall,sexl,nqv,i,j,sexd,int_;
		boolean tmp;
		char fff;
		
		boolean testhets_result;
		tmp=true;
		datafile >> numl >> b >> sexl >> prog >> GEOL;
		datafile >> a >> b >> c >> d >> GEOL;
		if (d==1)  { tmp=false; goto L10; }
		if ((prog!=1) && (prog!=3))  goto L10;
		datafile >> a >> GEOL;
		for( j=1; j <= numl; j ++) {
			datafile >> lc;    /*locus type*/
			switch (lc) {
				case 0:{
					datafile >> nall >> GEOL;
					datafile >> GEOL;   /*skip gene frequencies*/
					datafile >> nqv >> GEOL; /*number of traits*/
					for( i=1; i <= nqv; i ++) datafile >> a >> GEOL; /*genotype means*/
					datafile >> a >> GEOL; /*var cov matrix*/
					datafile >> a >> GEOL;   /*multiplier*/
				}
					break;
					
				case 1:{
					datafile >> nall >> GEOL;
					datafile >> a >> GEOL;
					datafile >> nall >> GEOL;
					if (sexl==0)  for( i=1; i <= nall; i ++) datafile >> a >> GEOL; else
						for( i=1; i <= (nall+nall); i ++) datafile >> a >> GEOL;
				}
					break;
					
				case 2:{
					datafile >> nall >> GEOL;
					for( i=1; i <= (nall+2); i ++) datafile >> a >> GEOL;
				}
					break;
				case 3:{
					datafile >> nall >> GEOL;
					datafile >> a >> GEOL;
				}
					break;
			}
		}
		datafile >> sexd >> int_ >> GEOL;
		if (sexd!=0)
		{	if (numl==1)  /*changed 10 Aug 93*/
			datafile >> GEOL;
		else datafile >> a >> GEOL; /*fem/male ratio, or female recombination values*/
		}
		if (numl==1)  /*changed 10 Aug 93*/
			datafile >> GEOL;
		else datafile >> a >> GEOL; /*recombination values*/
		numl=numl-1;
		if ((numl==2) && (int_==1))  numl=3;
		if ((sexd==1) && (numl > 0))  numl=numl+1;     /*changed 10 Aug 93*/
		if (sexd==2)  numl=numl+numl;
		datafile >> a >> GEOL;
		for( i=1; i <= numl; i ++) datafile >> d;
		fff=' ';
		while ((! is_EndOfLine(datafile.peek())) && (fff!='1'))  datafile >> fff;
		if (fff=='1')  tmp=false;
	L10: testhets_result=tmp;
		return testhets_result;
	}
	
	void initunknown()
	{
		if (dooutput)
		{
			std::cout << "Program UNKNOWN version " << version << std::endl;
			std::cout << "The following maximum values are in effect:" << std::endl;
			std::cout << fmt(maxlocus,8) << " loci" << std::endl;
			std::cout << fmt(maxgeno,8) << " single locus genotypes" << std::endl;
			std::cout << fmt(maxall,8) << " alleles at a single locus" << std::endl;
			std::cout << fmt(maxind,8) << " individuals in one pedigree" << std::endl;
			std::cout << fmt(maxmarriage,8) << " marriage(s) for one male" << std::endl;
			std::cout << fmt(maxtrait,8) << " quantitative factor(s) at a single locus" << std::endl;
			std::cout << fmt(maxliab,8) << " liability classes" << std::endl;
			std::cout << fmt(maxfact,8) << " binary codes at a single locus" << std::endl;
		}
	}
	
} // end namespace

void unknown_program(std::istream& datafile, std::istream& pedfile, std::ostream& ipedfile, std::ostream& speedfile, int& errors)
{
	iterpeds_errors=0;
	initunknown();

	makehomozygous=testhets(datafile); /*Change - lines added by Joe Terwilliger 7/8/93*/
	datafile.seekg(0,datafile.beg);
	readloci(datafile);

	nsequence=0; // used in reading pedfile.dat
	if (! is_blank_row(pedfile)) pedfile >> newped;
	while (! is_blank_row(pedfile))
	{
		readped(pedfile);
		getunknown();
		for( whichsys=1; whichsys <= nsystem; whichsys ++)
			if (mutsys!=whichsys)
			{
				{
					locusvalues& with = *thislocus[whichsys];
					fgeno=with.nallele*(with.nallele+1) / 2;
					if (sexlink) mgeno=with.nallele;
					else mgeno=fgeno;
				}
				getlocation(thislocus[whichsys]);
				iterpeds();
			}
		infer();
		writeped(ipedfile);
		writespeed(speedfile);
		reinit();
	}
	errors=iterpeds_errors;
}

#ifdef _MAKE_PROG

int main(int argc, const char* argv[])
{
	//	pio_initialize(argc, argv); // removed since I don't have "output <<" anymore.

	std::fstream datafile("datafile.dat",std::ios_base::in);
	std::fstream pedfile("pedfile.dat",std::ios_base::in);
	std::fstream ipedfile("ipedfile.dat",std::ios_base::out);
	std::fstream speedfile("speedfile.dat",std::ios_base::out);
	int unknown_errors=0;
	unknown_program(datafile,pedfile,ipedfile,speedfile,unknown_errors);
	datafile.close();
	pedfile.close();
	ipedfile.close();
	speedfile.close();
	if (unknown_errors) { std::cerr<<"There're "<<unknown_errors<<" incompatibilities.\n"; return 1; }
	return EXIT_SUCCESS;
}

#endif

int unknown_main()
{
	//	pio_initialize(argc, argv); // removed since I don't have "output <<" anymore.
	
	std::fstream datafile("datafile.dat",std::ios_base::in);
	std::fstream pedfile("pedfile.dat",std::ios_base::in);
	std::fstream ipedfile("ipedfile.dat",std::ios_base::out);
	std::fstream speedfile("speedfile.dat",std::ios_base::out);
	int unknown_errors=0;
	unknown_program(datafile,pedfile,ipedfile,speedfile,unknown_errors);
	datafile.close();
	pedfile.close();
	ipedfile.close();
	speedfile.close();
	if (unknown_errors) { std::cerr<<"There're "<<unknown_errors<<" incompatibilities.\n"; return 1; }
	return EXIT_SUCCESS;
}
