// link libs: -lnlopt_cxx

/*
 Features:
 1) MAF is not by MLE.
 2) Truely weighting, and can weight by Deleteriousness+Quality+GuiltByAssociation+Segregation.
 3) Use LINKAGE to detect Mendelian error rather than just trios.
 4) Calculating probability of damaging of a haplotype rather than counting pathogenic variants, so it's robust to LD.
 5) Provide different risk models (polygenic dominant recessive additive).
 6) Can handle counpound heterozygosity.
 7) Can do logistic regression and linear regression, so can be used for QTs.
 8) regression uses sample weights to control for familial correlations.
 9) allow user-defined R codes, so it can do cox regression or other analysis.
 10) write involve_MHC, so user can adjust for HLA gene depending on whether it's in the MHC region.
 11) All methods can include covariates, including FET.
 12) Use multithreading to do permutation.
 13) genome-wide permutation p-value, less conservative than Bonferroni and robust to overlapping variant sets (overlapping gene or pathway).
 14) Do not filter MAF in controls, but in all samples, which is unbiased.
 15) one-sided test is less noisy.
 16) provide options (--flip-aff --filt-MaxAF --filt-FdrAF --filt-del) to search for protective genes.
 17) jointly consider large deletions while analyzing small variants.
 18) gene set analysis
 
 Input: study data contain cases and controls.
 Input: must contain Func_Gene & is sorted. Use them to identify variant groups.
 Input: must contain Chr,Start,End,Ref,Alt. Use them to identify variant duplicates.
 If there're duplicated variants within a group, only the first is used (deleteriousness better be descending).
 It writes only rows whose case genotype is not empty, a useful feature for combined data (eg, Ps cs + BrCa ct).
 
 Caveats:
 1) Taking the first variant among duplicates requires input sorted by BayesDel.
 2) SSU includes matrix operations, which slow down very fast when the number of variants increases. Previously, I used cap (--cap=0.1) to solve it, but it doesn't work for GSA (BayesDel not sorted).
 
 Hidden options:
 --wt-del B         Weight variants by deleteriousness {_Default_wt_del}
 --wt-vqs B         Weight variants by call quality {_Default_wt_vqs}
 --biol FILE        Weight variants by biological relevance scores in FILE {_Default_biol}
 --neg-biol B       Weight variants by biological relevance even when the score is negative {_Default_neg_biol}
 --no-del           Turn off both --filt-del and --wt-del
 --add-info         Output to the INFO column instead of adding a new column
 --araf D           Aggregated frequency of risk alleles in a causal gene {_Default_araf}
 --cap D            Set a maximum number of variants per gene as NumberOfCases*D {_Default_cap}
 --top B            Set a maximum number of variants per person (1 for dominent, 2 for recessive/additive) {_Default_top}
 --opt-pen          Do HLR with optimized penetrance (the input penetrance will become just the starting point)
 --SSUc             Do Sum Squared U test, customized weights
 --permute I        Do permutation test with I repeats {_Default_permute}
 --one-file         The input Genotype File contains all variants. If set, --permute will calculate family-wise p-values, otherwise output maximum statistics to stdout. {_Default_one_file}
 --perm-data Fs     Read maximum statistics from Fs {_Default_perm_data}
 --out-pv B         Output p-value instead of log10 Bayes factor (for --logistic --linear --SSUx) {_Default_out_pv}
 --af Ss            Ordered allele frequency fields in INFO for polygenic odds ratio calculation (--por) {_Default_af}
 --prs B            using Polygenic Risk Score (requires PRS_beta and PRS_allele in INFO) {_Default_prs}
 --por B            using Polygenic Odds Ratio reference to general population (also requires PRS_freq in INFO or --af) {_Default_por}
 --ok-no-match B    skip the variants that have no match in the Genotype File {_Default_ok_no_match}
 
 Removed options:
 --recessive B      Do recessive model in FET {_Default_recessive}
 --polygenic        Do HLR using a polygenic model
 */

#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_ldt.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_mtjobs.hpp>
#include <tft/libfbj_fet.hpp>
#include <tft/libfbj_regress.hpp>
#include <nlopt.hpp>
#include "victor_par.hpp"
#include "vAAA_ssu.hpp"
#include "vAAA_GLR.hpp"
#include "vAAA_par.hpp"
#include "vAAA_if.hpp"
#include "vAAA_sv.hpp"
#include <mutex> // requires c++11. g++47 needs it! clang++4.0/g++44 doesn't.

using namespace std;
typedef genepi::genotype GTP;

// ----------- functions for doing jobs -----------

// These data must be set before any job and be constant afterward.
bool						showID = false;		// show SeqID for --detail
bool						showGT = false;		// show genotype for --detail
bool						showIF = false;		// show incidental findings for --detail
bool						AdjPrp = false;		// adjusted multinomial proportions. Def=false because a) 0 obs is not a problem, b) result is weired, see file5 file6
bool						XAF_fp = false;		// calcualte population frequency from MaxAF (not robust to LD), otherwise from data (not robust to small study)
bool						h_test = true;		// trend test, otherwise model test based on penetrance (could be recessive/additive/dominant)
bool						flr_lr = true;		// Firth's logistic regression use likelihood ratio test instead of wald test
bool						one_sided = true;	// do one-sided test by logistic / linear regression
bool						out_lp = false;		// output -log10(p)  for logistic / linear regression, otherwise convert p to LOD and output LOD
bool						out_pv = false;		// output p-value    for logistic / linear regression, otherwise convert p to LOD and output LOD
bool						out_or = false;		// output odds ratio for logistic only
bool						usePRS = false;		// output polygenic risk score
bool						usePOR = false;		// output polygenic odds ratio
bool						out_se = false;		// output beta and se of beta from FLR
bool						immedC = true;		// immediate continuity correction in Mantel-Haenszel OR calcualtion
bool						delayC = false;		// delay continuity correction in Mantel-Haenszel OR calcualtion
bool						constC = false;		// constant continuity correction in Mantel-Haenszel OR calcualtion (even for tables w/o a zero cell)
bool						FETrec = false;		// FET do recessive model
int							FETmin = 1;			// FET minimum cohort (case/control) size
int							FETiwt = false;		// FET use individual weight
int							RNKmin = 20;		// RNK minimum sample size (case+control)
int							denovo = 0;			// at least how many offspring has a de novo mutation, 0 = no filter
double						defMAF = 0.00001;	// min MaxAF if MaxAF=0. ExAC has 60706 samples, so 1/60706/2=0.00000823641815.
int							NumCSs;				// number of cases
int							NumCTs;				// number of controls
vector<double>				DepVar;				// Dependent variable  : DepVar[sampleID]
vector< vector<double> >	CovVar;				// covariates for REG  : CovVar[sampleID][variable]
vector<string>				CovNID;				// covariate new ID    : CovVar[variable]
vector< vector<int> >		PmtPtr;				// Permuted pointer    : IndPtr[repeatNo][sampleID] (repeatNo=0 is original pointer, repeatNo>0 is permutation pointer, size=1+permte)
vector<int>					OriPtr;				// Original pointer    : OriPtr[sampleID]=sampleID (not used anymore)
vector<string>				strata;				// strata def by cov   : strata[sampleID]
vector<string>				record;				// SequencedID         : record[sampleID]
vector<double>				sqrtWt;				// sqrt sample weight  : sqrtWt[sampleID]
vector<double>				origWt;				// effective No. chrom : origWt[sampleID]
map<string,int>				seq2co;				// record coordinate   : seq2co[SequenceID]=sampleID (coordinate in record[])
string						prData;				// to print X,y data to prData.x prData.y
string						MisStr = "NA";		// missing value string
string						MyCode;				// for MyR.
int							permte = 0;			// do within job permutations # repeats to get p-values
double CI = qnorms(0.975);

// About analysis jobs
enum AnalysisMethod_t { CL1, CL2, FET, RNK, DET, HLR, POL, REG, SSU, GLR, OPN, OUT, FLR, MyR } AnalysisMethod=FET; // previously HLR, REG, FET
struct JobData {
	int						pID;	// permutation ID. 0 means no permutation, 1+ for permutation.
	vector<string>			gtM;	// input gtM[sampleID][variant] = genotype matrix.
	vector<string>			gtC;	// input gtC[sampleID][variant] = copy number 0-9.
	vector<string>			qcl;	// input qcl[row_num] = quality control log.
	vector<string>			idx;	// input idx[variant] = variant index.
	vector<double>			wts;	// input wts[variant] = variant weights.
	vector<double>			maf;	// input maf[variant] = minor allele frequency.
	vector<double>			del;	// input del[variant] = BayesDel score
	vector<double>			vqs;	// input vqs[variant] = VQSLOD score
	vector<double>			xiU;	// input xiU[variant] = number of alt allele
	vector<double>			niU;	// input niU[variant] = number of non-missing chromosomes
	vector<double>			xiA;	// input xiA[variant] = number of alt allele
	vector<double>			niA;	// input niA[variant] = number of non-missing chromosomes
	vector<int>				ndn;	// input ndn[variant] = number of trios have de novo mutations
	vector<string>			gnm;	// input gnm[variant] = gene name
	vector<int>				chr;	// input chr[variant] = chr_num, used by GLR
	vector<int>				pos;	// input pos[variant] = 1-based bp
	vector<string>			ref;	// input ref[variant] = REF
	vector<string>			alt;	// input alt[variant] = ALT
	vector<char>			rsk;	// input rsk[variant] = risk allele, R=REF A=ALT
	string					grp;	// input this_group ID
	bool					MHC;	// input this_group has variant(s) in the MHC region
	string					out;	// output string, used only by CL2
	string					xtr;	// output @ INFO, used by most AnalysisMethod
	double					res;	// output result, used by most AnalysisMethod
	JobData():pID(0),MHC(false),out(MisStr),res(std::numeric_limits<double>::signaling_NaN()) { }
	void clear_result() {
		out=MisStr;
		xtr.clear();
		res=std::numeric_limits<double>::signaling_NaN();
	}
	void setup(int SS, const string& this_group) {
		grp=this_group;
		MHC=false;
		gtM.assign(SS,string());
		gtC.assign(SS,string());
		qcl.clear();
		idx.clear();
		wts.clear();
		maf.clear();
		del.clear();
		vqs.clear();
		xiU.clear();
		niU.clear();
		xiA.clear();
		niA.clear();
		ndn.clear();
		gnm.clear();
		chr.clear();
		pos.clear();
		ref.clear();
		alt.clear();
		rsk.clear();
		clear_result();
	}
	string chr_str() {
		string output;
		string previous_gene="_dummy_";
		for (size_t i=0;i<gnm.size();++i)
			if (gnm[i]!=previous_gene)
			{
				previous_gene=gnm[i];
				if (output.empty()) output =genepi::convert_chr_num(chr[i]);
				else { output+=','; output+=genepi::convert_chr_num(chr[i]); }
			}
		return output;
	}
} DefaultJob;

double truncated(double pen)
{
	if (pen>=1) return 0.9999999;
	if (pen<=0) return 0.0000001;
	return pen;
}

double penetrance_by_prs(const string& gt, const string& cp, const JobData& job, const vector<double>& pen) // only useful for unphased data, doesn't consider ploidy!=2
{
	double RR=1;
	double OR=pen[1]/pen[0];
	for (size_t l=0;l<gt.size();++l)
	{
		int ploidy = GTP::ploidy(gt[l]);
		if (ploidy==2)
		{
			double p = job.maf[l];
			double q = 1-p;
			double d = q*q + 2*p*q*pow(OR,job.wts[l]) + p*p*pow(OR,job.wts[l]*2);
			int	var = GTP::num_alt(gt[l]);
			if 		(cp[l]=='2')
			{
				RR = RR * pow(OR,var*job.wts[l]) / d;
			}
			else if (cp[l]=='1')
			{
				if 		(var==0) 		 { RR = RR * OR / d; }
				else if (var==2||var==1) { RR = RR * pow(OR,1+job.wts[l]) / d; }
				// else if (var==1) { lns<<showw<<"heterozygous genotype on a CNV region that the copy numbre is 1."<<flush_logger; }
			}
			else if (cp[l]=='0')
			{
				RR = RR * pow(OR,2) / d;
			}
		}
		else if (ploidy==1)
		{
			double p = job.maf[l];
			double q = 1-p;
			double d = q + p*pow(OR,job.wts[l]);
			int	var = GTP::num_alt(gt[l]);
			if 		(cp[l]=='2')
			{
				RR = RR * pow(OR,var*job.wts[l]) / d;
			}
			else if (cp[l]=='0')
			{
				RR = RR * OR / d;
			}
		}
	}
	return truncated(perch::preval*RR);
}

double polygenic_score(const string& gt, const string& cp, const JobData& job)
{
	double PRS=0; // polygenic risk score
	double GPS=0; // general population score
	for (size_t l=0;l<gt.size();++l)
	{
		int ploidy = GTP::ploidy(gt[l]);
		int	var = GTP::num_alt(gt[l]);
		if (ploidy==2)
		{
			if 		(cp[l]=='2') ;
			else if (cp[l]=='1') { if (var<2) ++var; }
			else if (cp[l]=='0') { var=2; }
		}
		else if (ploidy==1)
		{
			if 		(cp[l]=='2') ;
			else if (cp[l]=='0') var=1;
		}
		if (job.rsk[l]=='A') PRS += var*job.wts[l];
		else				 PRS += (ploidy-var)*job.wts[l];
		if (usePOR)
		{
			double p = job.maf[l];
			GPS += 2*p*job.wts[l]; // q*q*beta*0 + 2*p*q*beta*1 + p*p*beta*2 = 2*p*beta*(q+p) = 2*p*beta
		}
	}
	if (usePOR) return exp(PRS-GPS);
	else		return PRS;
}

void prob_damaging(const string& g, const string& cp, size_t g_beg, size_t g_end, const JobData& p, double& p_0D, double& p_1D, double& p_2D) // robust to partically-phased or unphased data
{
	if (g.empty()) { p_0D=1; p_1D=0; p_2D=0; return; }
	
	int ploidy = GTP::ploidy(g[g_beg]);
	if (ploidy==2)
	{
		vector< tuple<double,double,double> > probs = { make_tuple(1,1,1) }; // p(genotype),p(hap1_not_damaged),p(hap2_not_damaged)
		for (size_t l=g_beg;l<g_end;++l)
		{
			if (cp[l]=='2')
			{
				double p_1 = GTP::prob_1(g[l]);
				double p1_ = GTP::prob1_(g[l]);
				size_t s = probs.size();
				if (p_1>0 && p_1<1) // unphased. p_1+p1_=1. They mean the probability of a heterozygous variant being on hap1 and hap2, respectively.
				{
					for (size_t i=0;i<s;++i)
					{
						double v0 = std::get<0>(probs[i]);
						double v1 = std::get<1>(probs[i]);
						double v2 = std::get<2>(probs[i]);
						probs[i] =		make_tuple(v0*p_1, v1*(1-p.wts[l]), v2);
						probs.push_back(make_tuple(v0*p1_, v1, v2*(1-p.wts[l])));
					}
				}
				else // phased or homozygous
				{
					for (size_t i=0;i<s;++i)
					{
						double v0 = std::get<0>(probs[i]);
						double v1 = std::get<1>(probs[i]);
						double v2 = std::get<2>(probs[i]);
						probs[i] = make_tuple(v0, v1*(1-p_1*p.wts[l]), v2*(1-p1_*p.wts[l]));
					}
				}
			}
			else if (cp[l]=='1')
			{
				int a=GTP::num_alt(g[l]);
				if		(a==2 || a==1)
				{
					size_t s = probs.size();
					for (size_t i=0;i<s;++i)
					{
						double v0 = std::get<0>(probs[i]);
						double v1 = std::get<1>(probs[i]);
						probs[i] = make_tuple(v0, v1*(1-p.wts[l]), 0);
					}
				}
				else if (a==0)
				{
					size_t s = probs.size();
					for (size_t i=0;i<s;++i)
					{
						double v0 = std::get<0>(probs[i]);
						double v1 = std::get<1>(probs[i]);
						probs[i] = make_tuple(v0, v1, 0);
					}
				}
				else // remove a==1 in the above "(a==2 || a==1)"
				{
					// lns<<showw<<"heterozygous genotype on a CNV region that the copy numbre is 1."<<flush_logger;
				}
			}
			else if (cp[l]=='0')
			{
				size_t s = probs.size();
				for (size_t i=0;i<s;++i)
				{
					double v0 = std::get<0>(probs[i]);
					probs[i] = make_tuple(v0, 0, 0);
				}
			}
		}
		p_0D=0; p_1D=0; p_2D=0;
		for (auto &p:probs)
		{
			double v0 = std::get<0>(p);
			double v1 = std::get<1>(p);
			double v2 = std::get<2>(p);
			double both_benign = v1*v2;
			double both_damaging = (1-v1)*(1-v2);
			double otherwise = 1-both_benign-both_damaging;
			p_0D += v0 * both_benign;
			p_1D += v0 * otherwise;
			p_2D += v0 * both_damaging;
		}
		// cerr<<g<<' '<<p_0D<<' '<<p_1D<<' '<<p_2D<<endl;
	}
	else if (ploidy==1)
	{
		double pD00=1;
		for (size_t l=g_beg;l<g_end;++l)
		{
			if 		(cp[l]=='2') pD00 *= ( 1 - GTP::prob_1(g[l]) * p.wts[l] );
			else if (cp[l]=='0') pD00 *= 0;
		}
		double pD01 = 1 - pD00; // prob damaging chromatid a
		p_0D = 1 - pD01; // probability of BenignBenign
		p_2D =	   pD01; // probability of 2 damaging chromosomes
		p_1D = 0;		 // probability of 1 damaging chromosomes
		// cerr<<g<<' '<<p_0D<<' '<<p_1D<<' '<<p_2D<<endl;
	}
	else
	{
		p_0D = 1;
		p_2D = 0;
		p_1D = 0;
		// cerr<<g<<' '<<p_0D<<' '<<p_1D<<' '<<p_2D<<endl;
	}
}//*/

void sum_damaging(const string& g, const string& cp, const JobData& p, double& p_0D, double& p_1D, double& p_2D)
{
	p_0D = 0;
	p_1D = 0;
	p_2D = 0;
	if (g.empty()) return;
	for (size_t g_beg=0,g_end=g_beg+1; g_end<=g.size(); ++g_end)
	{
		if (g_end==g.size() || p.gnm[g_beg]!=p.gnm[g_end])
		{
			double p0,p1,p2;
			prob_damaging(g,cp,g_beg,g_end,p,p0,p1,p2);
			p_0D+=p0;
			p_1D+=p1;
			p_2D+=p2;
			g_beg=g_end;
		}
	}
}

void max_damaging(const string& g, const string& cp, const JobData& p, double& p_0D, double& p_1D, double& p_2D)
{
	p_0D = 1;
	p_1D = 0;
	p_2D = 0;
	if (g.empty()) return;
	for (size_t g_beg=0,g_end=g_beg+1; g_end<=g.size(); ++g_end)
	{
		if (g_end==g.size() || p.gnm[g_beg]!=p.gnm[g_end])
		{
			double p0,p1,p2;
			prob_damaging(g,cp,g_beg,g_end,p,p0,p1,p2);
			p_0D=std::min(p_0D,p0);
			p_1D=std::max(p_1D,p1);
			p_2D=std::max(p_2D,p2);
			g_beg=g_end;
		}
	}
	if 		(	  p_1D+p_2D>=1) { p_1D=1-p_2D; p_0D=0; }
	else if	(p_0D+p_1D+p_2D>=1) { p_0D=1-p_2D-p_1D; }
}

double num_risk_elements(const string& g, const string& cp, const JobData& p)
{
	double p_0D,p_1D,p_2D;
	if (h_test) sum_damaging(g,cp,p,p_0D,p_1D,p_2D); // sum of risk haplotypes in multiple genes, could be >2
	else		max_damaging(g,cp,p,p_0D,p_1D,p_2D); // most severe genotype among multiple genes, 0 to 2
	return p_1D+2*p_2D;
}

int genotype_code(const string& g, const string& cp, const JobData& p)
{
	double p_0D,p_1D,p_2D;
	max_damaging(g,cp,p,p_0D,p_1D,p_2D);
	int genotype = -1;
	double maxp = std::max( p_0D, std::max(p_1D,p_2D));
	if		(p_2D==maxp) genotype=2; // prefer 2 over 1 because it's more likely that 2 RVs are not on the same haplotype
	else if (p_1D==maxp) genotype=1;
	else if (p_0D==maxp) genotype=0;
	else				 genotype=-1;
	return genotype;
}

double penetrance_by_gtp(const string& g, const string& cp, const JobData& p, const vector<double>& pen)
{
	double p_0D,p_1D,p_2D;
	max_damaging(g,cp,p,p_0D,p_1D,p_2D);
	return truncated(pen[0]*p_0D + pen[1]*p_1D + pen[2]*p_2D);
}

double risk_score(const string& g, const string& cp, const JobData& p)
{
	if		(usePRS) return polygenic_score(g, cp, p);
	else if (usePOR) return polygenic_score(g, cp, p);
	else if (h_test) return num_risk_elements(g, cp, p);
	else			 return penetrance_by_gtp(g, cp, p, perch::penetr);
}

double (*penetrance)(const string&, const string&, const JobData&, const vector<double>&) = &penetrance_by_gtp;

vector<string>				AllDataSym; // AllDataSym[GeneID] = Gene Symbol
vector< vector<double> >	AllDataGtp; // AllDataGtp[Sample][GeneID] = probability of damaging

// regression functions
double (*regress)(const Eigen::MatrixXd	& X, const Eigen::VectorXd	& y, bool one_sided) = &pv_1st_linear;
double (*do_SSUw)(const Eigen::MatrixXd	& X, const Eigen::VectorXd	& y, const std::vector<double> &) = &pv_SSUc;

// for within job permutations
Values<double>				wjp_res_mx;			// max result among all genes within a permutation
map< int, Values<double> >	wjp_res_db;			// result database [permutation][results]
std::mutex					wjp_res_db_mutex;	// mutex for

// for between job multithreading
map< string, tuple<string,string,double> >	bjm_res_db;			// result database [permutation][results]
std::mutex									bjm_res_db_mutex;	// mutex for
MtJobs<JobData>								bjmJobs;			// one job for one gene

struct CGinfo {				// collapsed genotype info
	int					cs;	// number of cases
	int					ct;	// number of controls
	double				pf;	// population frequency calculated as min(geno_freq_by_MaxAF) across variants in this genotype
	CGinfo():cs(0),ct(0) {}
};
struct STinfo {				// strata info
	int					cs; // number of cases
	int					ct; // number of controls
	map< tuple<string,string>, CGinfo>	cg;	// collapsed genotype
	STinfo():cs(0),ct(0) {}
};
struct OPNinput
{
	map<string, STinfo> * sg;
	const JobData * p;
};

double _cal_OPN(const std::vector<double> &x, std::vector<double> &grad, void* f_data)
{
	OPNinput* input = static_cast< OPNinput* >(f_data);
	if (x[1]==1 && x[2]==1) return -DBL_MAX;
	vector<double> newPen = { truncated(x[0]), truncated(x[0]*x[1]), truncated(x[0]*x[1]*x[2]) };
	double score=0;	// previously log10LR
	for (auto &s:*input->sg)
	{
		if (!s.second.cs || !s.second.ct) continue;
		if (s.second.cg.size()<2) continue;
		
		vector<double> fp,pn,fs,ft,ns,nt;
		for (auto &g:s.second.cg)
		{
			pn.push_back( penetrance(std::get<0>(g.first), std::get<1>(g.first), *input->p, newPen) );
			ns.push_back( g.second.cs );
			nt.push_back( g.second.ct );
			if (XAF_fp)	fp.push_back( g.second.pf );
			else		fp.push_back( perch::preval * g.second.cs/s.second.cs + (1-perch::preval) * g.second.ct/s.second.ct );
			// cerr<<g.first<<'\t'<<ns.back()<<'\t'<<nt.back()<<'\t'<<fp.back()<<'\t'<<pn.back()<<endl;
		}
		MakeFraction(fp);
		genepi::up_dn_af(fp,pn,fs,ft);
		double llr=0; // log10 likelihood ratio (LOD) for this strata
		for (size_t i=0;i<fp.size();++i)
		{
			llr +=		  ns[i]  * log10(fs[i]) ;
			llr +=		  nt[i]  * log10(ft[i]) ;
			llr -= (ns[i]+nt[i]) * log10(fp[i]) ;
			// cerr<<ns[i]<<'/'<<s.second.cs<<'\t'<< nt[i]<<'/'<<s.second.ct<<'\t'<<fp[i]<<'\t'<<fs[i]<<'\t'<<ft[i]<<'\t'<<llr<<endl;
		}
		score += llr;
	}
	// cerr<<str_of_container(newPen,',',false)<<' '<<score<<endl;
	return score;
}

void _run_job(JobData& p)
{
	string overall_idx=str_of_container(p.idx,',');
	bjm_res_db_mutex.lock();
	map<string,tuple<string,string,double> >::iterator dbit = bjm_res_db.find(overall_idx);
	if (dbit!=bjm_res_db.end())
	{
		p.out=std::get<0>(dbit->second);
		p.xtr=std::get<1>(dbit->second);
		p.res=std::get<2>(dbit->second);
		bjm_res_db_mutex.unlock();
		return;
	}
	bjm_res_db_mutex.unlock();
	
	// need to cal
	p.clear_result();
	string uID=random_string(12);
	vector<int>& ptr=PmtPtr[p.pID];
	
	// skip by de novo count
	if (denovo)
	{
		int observed_denovo=0;
		for (auto &n:p.ndn) observed_denovo+=n;
		if (observed_denovo<denovo)
		{
			p.out="SKIPPED=DeNovoFilter";
			return;
		}
	}
	
	// skip if not SNV/InDel observed.
	if (p.idx.empty())
	{
		p.out="SKIPPED=NoVariant";
		return;
	}
	
	// prepare a new gtM but reserve the p.gtM
	vector<string> newGtM; // has to contain all samples even with missings
	vector<string> newGtC; // has to contain all samples even with missings
	for (auto &i:ptr)
	{
		newGtM.push_back(p.gtM[i]);
		newGtC.push_back(p.gtC[i]);
	}
	
	// for (size_t i=0;i<record.size();++i) cerr<<record[i]<<'\t'<<DepVar[i]<<'\t'<<newGtM[i]<<endl;
	
	// analysis
	if (!ExtCT::ExAC_pfx.empty())
	{
		struct var_summary {
			string idx;
			double xiA;
			double niA;
			double xiU;
			double niU;
			double wt;
			int		cohort; // -1=unknown 0=ExAC 1=study 2=both
			var_summary():xiA(0),niA(0),xiU(0),niU(0),wt(1),cohort(-1){}
		};
		vector<var_summary> to_cal;
		for (size_t g_beg=0,g_end=g_beg+1; g_end<=p.idx.size(); ++g_end)
		{
			if (g_end==p.idx.size() || p.gnm[g_beg]!=p.gnm[g_end])
			{
				if (exist_element(ExtCT::ExAC_dta,p.gnm[g_beg]))
				{
					ExtCT::clear_added(p.gnm[g_beg]);
					double ExAC_ftr=0, Case_ftr=0, Ctrl_ftr=0; // not considering PAR yet
					if		(genepi::is_autosomal(p.chr[g_beg])) { ExAC_ftr=2;					Case_ftr=2;					Ctrl_ftr=2;					}
					else if (genepi::is_chrX	 (p.chr[g_beg])) { ExAC_ftr=1+ExtCT::ExAC_fem;	Case_ftr=1+ExtCT::Case_fem; Ctrl_ftr=1+ExtCT::Ctrl_fem; }
					else if (genepi::is_chrY	 (p.chr[g_beg])) { ExAC_ftr=1-ExtCT::ExAC_fem;	Case_ftr=1-ExtCT::Case_fem; Ctrl_ftr=1-ExtCT::Ctrl_fem; }
					else								 		 { ExAC_ftr=1;					Case_ftr=1;					Ctrl_ftr=1;					}
					
					for (size_t i=g_beg;i<g_end;++i)
					{
						var_summary v;
						v.idx=p.idx[i];
						v.xiU=(ExtCT::ExAC_cso ? 0 : p.xiU[i]);
						v.niU=(ExtCT::ExAC_cso ? 0 : p.niU[i]);
						v.xiA=p.xiA[i];
						v.niA=p.niA[i];
						double max_vqs=0;
						double exac_xi=0;
						double exac_ni=0;
						if (exist_element(ExtCT::ExAC_dta[p.gnm[g_beg]],p.idx[i])) // variant in ExAC and study
						{
							p.qcl[i]+=";ExAC_AC="+ftos(ExtCT::ExAC_dta[p.gnm[g_beg]][p.idx[i]].xiU);
							p.qcl[i]+=";ExAC_AN="+ftos(ExtCT::ExAC_dta[p.gnm[g_beg]][p.idx[i]].niU);
							if (ExtCT::ExAC_dta[p.gnm[g_beg]][p.idx[i]].DoNotUse) { ExtCT::wr_log(p.gnm[g_beg],p.idx[i],"SKIPPED=DoNotUse"); p.qcl[i]+=";SKIPPED=DoNotUse"; continue; } // should not happen
							v.cohort=2;
							exac_xi=ExtCT::ExAC_dta[p.gnm[g_beg]][p.idx[i]].xiU;
							exac_ni=ExtCT::ExAC_dta[p.gnm[g_beg]][p.idx[i]].niU;
							max_vqs=std::max(ExtCT::ExAC_dta[p.gnm[g_beg]][p.idx[i]].vqs,p.vqs[i]);
							ExtCT::ExAC_dta[p.gnm[g_beg]][p.idx[i]].added=true;
						}
						else // variant in study only
						{
							double covered_by_ExAC=ExtCT::covered_samples(ExtCT::ExAC_spl,p.chr[i],p.pos[i],p.ref[i],p.alt[i]);
							if (!covered_by_ExAC) { ExtCT::wr_log(p.gnm[g_beg],p.idx[i],"SKIPPED=study_var_not_covered_by_ExAC"); p.qcl[i]+=";SKIPPED=study_var_not_covered_by_ExAC"; continue; }
							v.cohort=1;
							exac_xi=0;
							exac_ni=ExAC_ftr*ExtCT::ExAC_sub*covered_by_ExAC;
							max_vqs=p.vqs[i];
						}
						v.xiU+=exac_xi;
						v.niU+=exac_ni;
						if (v.niA==0 || v.niU==0) { ExtCT::wr_log(p.gnm[g_beg],p.idx[i],"SKIPPED=study_var_no_observation"); p.qcl[i]+=";SKIPPED=study_var_no_observation"; continue; }
						if (v.xiA==0 && v.xiU==0) { ExtCT::wr_log(p.gnm[g_beg],p.idx[i],"SKIPPED=study_var_no_variation"); p.qcl[i]+=";SKIPPED=study_var_no_variation"; continue; }
						if (v.xiA==v.niA && v.xiU==v.niU) { ExtCT::wr_log(p.gnm[g_beg],p.idx[i],"SKIPPED=study_var_no_variation"); p.qcl[i]+=";SKIPPED=study_var_no_variation"; continue; }
						if (ExtCT::ExAC_fqc)
						{
							int table[4];
							table[0] = exac_ni - exac_xi + 0.5;
							table[1] = exac_xi + 0.5;
							table[2] = p.niA[i] + p.niU[i] - p.xiA[i] - p.xiU[i];
							table[3] = p.xiA[i] + p.xiU[i];
							FisherExactTest fet(true);
							fet.input_RxC(table,2,2);
							double pv = fet.result('2');
							if (pv<ExtCT::ExAC_fqc) { ExtCT::wr_log(p.gnm[g_beg],p.idx[i],"SKIPPED=study_var_allele_freq_qc"); p.qcl[i]+=";SKIPPED=study_var_allele_freq_qc"; continue; }
						}
						v.wt = posterior_given_log10BF( p.del[i] + (max_vqs<0?max_vqs:0) );
						to_cal.push_back(v);
						ExtCT::wr_log(p.gnm[g_beg],p.idx[i],p.qcl[i]);
						// cerr<<v.idx<<'\t'<<p.del[i]<<'\t'<<max_vqs<<'\t'<<v.wt<<'\t'<<v.xiA<<'\t'<<v.niA<<'\t'<<v.xiU<<'\t'<<v.niU<<endl;
					}
					for (auto &d:ExtCT::ExAC_dta[p.gnm[g_beg]]) // variant in ExAC only
					{
						if (d.second.added) { continue; }
						if (d.second.DoNotUse) { ExtCT::wr_log(p.gnm[g_beg],d.second.idx,"SKIPPED=ExAC_var_filtered_out_by_study"); continue; }
						if (exist_element(ExtCT::studyQC,d.second.idx)) { ExtCT::wr_log(p.gnm[g_beg],d.second.idx,"SKIPPED=ExAC_var_removed_by_study_QC"); continue; } // should not happen
						var_summary v;
						v.cohort=0;
						v.idx=d.first;
						v.xiA=0;
						v.niA=0;
						v.xiU=d.second.xiU;
						v.niU=d.second.niU;
						double covered_by_study=ExtCT::covered_samples(ExtCT::StudySpl,d.second.chr,d.second.bp,d.second.ref,d.second.alt);
						if (!covered_by_study) { ExtCT::wr_log(p.gnm[g_beg],d.second.idx,"SKIPPED=ExAC_var_not_covered_by_study"); continue; }
						v.niA += 						Case_ftr*NumCSs*covered_by_study;
						v.niU += (ExtCT::ExAC_cso ? 0 : Ctrl_ftr*NumCTs*covered_by_study);
						if (v.niA==0 || v.niU==0) { ExtCT::wr_log(p.gnm[g_beg],d.second.idx,"SKIPPED=ExAC_var_no_observation"); continue; }
						if (v.xiA==0 && v.xiU==0) { ExtCT::wr_log(p.gnm[g_beg],d.second.idx,"SKIPPED=ExAC_var_no_variation"); continue; }
						if (v.xiA==v.niA && v.xiU==v.niU) { ExtCT::wr_log(p.gnm[g_beg],d.second.idx,"SKIPPED=ExAC_var_no_variation"); continue; }
						if (ExtCT::ExAC_fqc)
						{
							int table[4];
							table[0] = d.second.niU - d.second.xiU;
							table[1] = d.second.xiU;
							table[2] = Case_ftr*NumCSs*covered_by_study + Ctrl_ftr*NumCTs*covered_by_study + 0.5;
							table[3] = 0;
							FisherExactTest fet(true);
							fet.input_RxC(table,2,2);
							double pv = fet.result('2');
							if (pv<ExtCT::ExAC_fqc) { ExtCT::wr_log(p.gnm[g_beg],d.second.idx,"SKIPPED=ExAC_var_allele_freq_qc"); continue; }
						}
						v.wt = posterior_given_log10BF( d.second.del + (d.second.vqs<0?d.second.vqs:0) );
						to_cal.push_back(v);
						ExtCT::wr_log(p.gnm[g_beg],d.second.idx,"ExAC_AC="+ftos(d.second.xiU)+";ExAC_AN="+ftos(d.second.niU));
						// cerr<<v.idx<<'\t'<<d.second.del<<'\t'<<d.second.vqs<<'\t'<<v.wt<<'\t'<<v.xiA<<'\t'<<v.niA<<'\t'<<v.xiU<<'\t'<<v.niU<<endl;
					}
				}
				g_beg=g_end;
			}
		}
		if (AnalysisMethod==GLR)
		{
			double Lambda=0;	// same as that in the VAAST paper
			double Sum_Wt=0;	// sum or variant weights, will be used to normalize weights
			for (auto &v:to_cal)
			{
				double& xiA=v.xiA;
				double& niA=v.niA;
				double& xiU=v.xiU;
				double& niU=v.niU;
				if (xiU/niU>0.5) { xiA=niA-xiA; xiU=niU-xiU; }					// convert to minor allele
				double piP = xiA/niA*perch::preval + xiU/niU*(1-perch::preval);	// allele frequency in population, won't be 0 because "if (!xiA && !xiU) continue;"
				double piA = genepi::up_af(piP,perch::penetr);					// allele frequency in affected
				double piU = genepi::dn_af(piP,perch::penetr);					// allele frequency in unaffected
				double niT = niA + niU;
				double xiT = xiA + xiU;
				double this_l = 0;
				{	this_l = this_l + xiT*log10(piP) + (niT-xiT)*log10(1-piP) ; }
				{	this_l = this_l - xiU*log10(piU) - (niU-xiU)*log10(1-piU) ; }
				{	this_l = this_l - xiA*log10(piA) - (niA-xiA)*log10(1-piA) ; }
				Lambda += this_l * v.wt;
				Sum_Wt += v.wt;
				// cerr<<v.idx<<'\t'<<xiA<<' '<<niA<<' '<<xiU<<' '<<niU<<' '<<piA<<' '<<piU<<' '<<piP<<' '<<this_l<<' '<<Lambda<<' '<<Sum_Wt<<endl;
			}
			if (Sum_Wt)
			{
				double score = -Lambda/Sum_Wt;
				p.res = score;
				if		(out_lp) { double pv = genepi::LOD_to_P(score); p.out=ftos(-log10(pv),MisStr); }
				else if (out_pv) { double pv = genepi::LOD_to_P(score); p.out=ftos(pv		 ,MisStr); }
				else			 {										p.out=ftos(score	 ,MisStr); }
			}
			else
				p.out = MisStr;
		}
		else
		{
			int 	cc[3]={0,0,0};
			double 	agg_xiA=0;
			double	agg_xiU=0;
			double 	agg_riA=0;
			double 	agg_riU=0;
			for (auto &v:to_cal)
			{
				agg_xiA +=  v.xiA;
				agg_xiU +=  v.xiU;
				agg_riA += (v.niA-v.xiA);
				agg_riU += (v.niU-v.xiU);
				++cc[v.cohort];
			}
			int table[4];
			table[0] = agg_riU+0.5;	double HN=table[0]; // Healthy not-exposed
			table[1] = agg_xiU+0.5;	double HE=table[1]; // Healthy Exposed
			table[2] = agg_riA+0.5;	double DN=table[2]; // Diseased not-exposed
			table[3] = agg_xiA+0.5;	double DE=table[3]; // Diseased Exposed
			if ( (DE==0&&HE==0) || (DN==0&&HN==0) || (DE==0&&DN==0) || (HE==0&&HN==0) )
			{
				if (AnalysisMethod==FET) p.out = MisStr;
				else					 p.out = p.chr_str()+"\t"+itos(to_cal.size())+"\tNA\t"+itos(table[2])+"\t"+itos(table[3])+"\t"+itos(table[0])+"\t"+itos(table[1])+"\t"+itos(cc[0])+"\t"+itos(cc[1])+"\t"+itos(cc[2]);
			}
			else
			{
				FisherExactTest fet(true);
				fet.input_RxC(table,2,2);
				double pv = one_sided ? fet.result('R') : fet.result('2') ;
				double mlp = -log10(pv);
				p.res = mlp;
				
				if (AnalysisMethod==FET)
				{
					if		(out_lp) p.out = ftos(mlp					,MisStr);
					else if (out_pv) p.out = ftos(pv					,MisStr);
					else			 p.out = ftos(genepi::P_to_LOD(pv)	,MisStr);
				}
				else // CL1 CL2
				{
					p.out = p.chr_str()+"\t"+itos(to_cal.size())+"\t"+
							ftos(mlp)+"\t"+
							itos(table[2])+"\t"+itos(table[3])+"\t"+itos(table[0])+"\t"+itos(table[1])+"\t"+
							itos(cc[0])+"\t"+itos(cc[1])+"\t"+itos(cc[2]);
				}
			}
		}
	}
	else if (AnalysisMethod==HLR || AnalysisMethod==POL || AnalysisMethod==OPN)
	{
		// convert data, exclude records with missing genotypes
		map<string, STinfo>		sg;	// stratified genotypes sg[strata] (strata defined by relative risk)
		double t_cs=0, t_ct=0;		// total number of samples
		for (size_t i=0;i<newGtM.size();++i)
		{
			if (!GTP::usable(newGtM[i])) continue;
			STinfo& s = sg[strata[i]];
			if (DepVar[i])	{ ++s.cg[make_tuple(newGtM[i],newGtC[i])].cs; ++s.cs; ++t_cs; }
			else			{ ++s.cg[make_tuple(newGtM[i],newGtC[i])].ct; ++s.ct; ++t_ct; }
			double geno_frq=1;
			for (size_t v=0;v<newGtM[i].size();++v)
			{
				int ploidy=GTP::ploidy(newGtM[i][v]);
				int g=GTP::num_alt(newGtM[i][v]);
				if (ploidy==2)
				{
					if (g==1) { double f=(p.maf[v]?p.maf[v]:defMAF); geno_frq*=(f*(1-f)*2); }
					if (g==2) { double f=(p.maf[v]?p.maf[v]:defMAF); geno_frq*=(f*f); }
				}
				if (ploidy==1)
				{
					if (g==1) { double f=(p.maf[v]?p.maf[v]:defMAF); geno_frq*=f; }
				}
			}
			for (size_t g_beg=0,g_end=g_beg+1; g_end<=p.gnm.size(); ++g_end)
			{
				if (g_end==p.gnm.size() || p.gnm[g_beg]!=p.gnm[g_end])
				{
					for (size_t v=g_beg;v<g_end;++v)
						if (newGtC[i][v]==1||newGtC[i][v]==0) { geno_frq*=defMAF; break; }
					g_beg=g_end;
				}
			}
			s.cg[make_tuple(newGtM[i],newGtC[i])].pf=geno_frq;
		}
		
		if (!t_cs || !t_ct)
		{
			p.res = std::numeric_limits<double>::signaling_NaN();
			p.out = MisStr;
		}
		else
		{
			if (AnalysisMethod==OPN)
			{
				static const std::vector<double> lb = { 0.0000001,   1,   1 };
				static const std::vector<double> ub = { 0.9999999, 100, 100 };
				try
				{
					OPNinput f_data;
					f_data.sg = &sg;
					f_data.p= &p;
					double score;
					nlopt::opt oo(nlopt::LN_NELDERMEAD, 3); // LN_NELDERMEAD > LN_SBPLX > LN_BOBYQA >> LN_COBYLA LN_AUGLAG LN_PRAXIS GN_ISRES
					std::vector<double> x = { perch::penetr[0], perch::penetr[1]/perch::penetr[0], perch::penetr[2]/perch::penetr[1] };
					std::vector<double> g;
					oo.set_max_objective(_cal_OPN, &f_data);
					oo.set_lower_bounds(lb);
					oo.set_upper_bounds(ub);
					nlopt::result nlopt_result=oo.optimize(x, score);
					if (nlopt_result<0) { stringstream ss; ss<<"Optimization failed with code "<<nlopt_result; exit_error(ss.str()); }
					p.res = score;
					if		(out_lp) { double pv = genepi::LOD_to_P(score,3); p.out=ftos(-log10(pv),MisStr); } // not sure whether df=3
					else if (out_pv) { double pv = genepi::LOD_to_P(score,3); p.out=ftos(pv        ,MisStr); } // not sure whether df=3
					else			 {										  p.out=ftos(score     ,MisStr); }
				}
				catch (int e) // There is an error/warning code e.
				{
					p.res = std::numeric_limits<double>::signaling_NaN();
					p.out = MisStr;
				}
			}
			else
			{
				double ps_cs=0, ps_ct=0;	// sth to add to denominator : prior strength
				double st_cs=0, st_ct=0;	// sth to add to nominator : ps*t, t is non-informative prior for multinomial (t=1/k, k=#categories)
				double score=0;				// previously log10LR
				for (auto &s:sg)
				{
					if (!s.second.cs || !s.second.ct) continue;
					if (s.second.cg.size()<2) continue;
					if (AdjPrp)
					{
						{ int N=s.second.cs; int k=s.second.cg.size(); ps_cs=prior_strength(N,k); st_cs=ps_cs/k; }
						{ int N=s.second.ct; int k=s.second.cg.size(); ps_ct=prior_strength(N,k); st_ct=ps_ct/k; }
					}
					
					vector<double> fp,pn,fs,ft,ns,nt;
					for (auto &g:s.second.cg)
					{
						pn.push_back( penetrance(std::get<0>(g.first), std::get<1>(g.first), p, perch::penetr) );
						ns.push_back( g.second.cs );
						nt.push_back( g.second.ct );
						if (XAF_fp)	fp.push_back( g.second.pf );
						else		fp.push_back( (g.second.cs+st_cs)/(s.second.cs+ps_cs) * perch::preval + (g.second.ct+st_ct)/(s.second.ct+ps_ct) * (1-perch::preval) );
						// cerr<<g.first<<'\t'<<ns.back()<<'\t'<<nt.back()<<'\t'<<fp.back()<<'\t'<<pn.back()<<endl;
					}
					MakeFraction(fp);
					genepi::up_dn_af(fp,pn,fs,ft);
					double llr=0; // log10 likelihood ratio (LOD) for this strata
					for (size_t i=0;i<fp.size();++i)
					{
						llr +=		  ns[i]  * log10(fs[i]) ;
						llr +=		  nt[i]  * log10(ft[i]) ;
						llr -= (ns[i]+nt[i]) * log10(fp[i]) ;
						// cerr<<ns[i]<<'/'<<s.second.cs<<'\t'<< nt[i]<<'/'<<s.second.ct<<'\t'<<fp[i]<<'\t'<<fs[i]<<'\t'<<ft[i]<<'\t'<<llr<<endl;
					}
					score += llr;
				}
				p.res = score;
				if		(out_lp) { double pv = genepi::LOD_to_P(score); p.out=ftos(-log10(pv),MisStr); }
				else if (out_pv) { double pv = genepi::LOD_to_P(score); p.out=ftos(pv        ,MisStr); }
				else			 {										p.out=ftos(score     ,MisStr); }
			}
		}
	}
	else if (AnalysisMethod==CL1) // count genotype
	{
		vector<int> cgA = { 0,0,0 };	// count 0/0 0/1 1/1 for affected
		vector<int> cgU = { 0,0,0 };	// count 0/0 0/1 1/1 for unaffected
		int  num_variants = 0;
		for (size_t i=0;i<newGtM.size();++i)
		{
			if (!GTP::usable(newGtM[i])) continue;
			int nv=0;
			int het=0, hom=0;
			for (size_t v=0;v<newGtM[i].size();++v)
			{
				++nv;
				int a=GTP::num_alt_recessive(newGtM[i][v]);
				if (a==1) ++het;
				if (a==2) ++hom;
			}
			int g = (hom ? 2 : (het ? 1 : 0));
			if (DepVar[i])	++cgA[g];
			else			++cgU[g];
			num_variants=nv;
		}
		int table[4];
		table[0] = cgU[0];			double HN=table[0]; // Healthy not-exposed
		table[1] = cgU[1] + cgU[2];	double HE=table[1]; // Healthy Exposed
		table[2] = cgA[0];			double DN=table[2]; // Diseased not-exposed
		table[3] = cgA[1] + cgA[2];	double DE=table[3]; // Diseased Exposed
		if ( (DE==0&&HE==0) || (DN==0&&HN==0) || (DE==0&&DN==0) || (HE==0&&HN==0) )
		{
			p.out = p.chr_str()+"\t"+itos(num_variants)+"\tNA\tNA\t[NA,NA]\t"+str_of_container(cgA,'\t')+"\t"+str_of_container(cgU,'\t');
		}
		else
		{
			FisherExactTest fet(true);
			fet.input_RxC(table,2,2);
			double mlp = one_sided ? -log10(fet.result('R')) : -log10(fet.result('2')) ;
			p.res = mlp;
			
			DE+=0.5;
			HE+=0.5;
			DN+=0.5;
			HN+=0.5;
			double OR = (DE*HN)/(DN*HE);
			double SE_logOR=sqrt(1/DE+1/HE+1/DN+1/HN);
			
			p.out = p.chr_str()+"\t"+itos(num_variants)+"\t"+
					ftos(mlp)+"\t"+ftos(OR)+"\t["+ftos(exp(log(OR)-CI*SE_logOR))+","+ftos(exp(log(OR)+CI*SE_logOR))+"]"+"\t"+
					str_of_container(cgA,'\t')+"\t"+str_of_container(cgU,'\t');
		}
	}
	else if ((AnalysisMethod==CL2 || AnalysisMethod==FET) && !FETiwt) // count damaged haplotypes || Fisher's exact test with mid-p method by Irwin's rule
	{
		struct SubGroupData {
			vector<int> cgA;			// count 0/0 0/1 1/1 for affected
			vector<int> cgU;			// count 0/0 0/1 1/1 for unaffected
			SubGroupData() { cgA={0,0,0}; cgU={0,0,0}; }
		};
		map<string,SubGroupData> SubGrp;
		for (size_t i=0;i<newGtM.size();++i)
		{
			if (!GTP::usable(newGtM[i])) continue;
			int g = genotype_code(newGtM[i], newGtC[i], p);
			if (g>=0)
			{
				if (DepVar[i])	{ ++SubGrp[strata[i]].cgA[g]; }
				else			{ ++SubGrp[strata[i]].cgU[g]; }
			}
		}
		
		vector<int> cgA = { 0,0,0 };	// count 0/0 0/1 1/1 for affected
		vector<int> cgU = { 0,0,0 };	// count 0/0 0/1 1/1 for unaffected
		set<string> to_remove;
		for (auto &x:SubGrp)
		{
			int controls=x.second.cgU[0]+x.second.cgU[1]+x.second.cgU[2];
			int cases=x.second.cgA[0]+x.second.cgA[1]+x.second.cgA[2];
			if (cases<FETmin || controls<FETmin) {
				if (perch::_Debug) lns<<showl<<"removed "<<x.first<<' '<<cases<<' '<<controls<<flush_logger;
				to_remove.insert(x.first);
				continue; }
			for (int g=0;g<3;++g) { cgA[g]+=x.second.cgA[g]; cgU[g]+=x.second.cgU[g]; }
		}
		for (auto &x:to_remove) SubGrp.erase(x);
		
		double HN = ( FETrec ? cgU[0]+cgU[1] : cgU[0]		 ); // Healthy not-exposed
		double HE = ( FETrec ? cgU[2]        : cgU[1]+cgU[2] ); // Healthy Exposed
		double DN = ( FETrec ? cgA[0]+cgA[1] : cgA[0]		 ); // Diseased not-exposed
		double DE = ( FETrec ? cgA[2]        : cgA[1]+cgA[2] );	// Diseased Exposed
		bool cc = immedC || (delayC && (HN==0||HE==0||DN==0||DE==0));
		double overall_n=0; // nominator for overall p-value calculation
		double overall_d=0; // denominator for overall p-value calculation
		double overall_z=std::numeric_limits<double>::signaling_NaN(); // z score for overall p-value calculation
		double overall_p=std::numeric_limits<double>::signaling_NaN();
		double mh_n=0; // nominator for Mantel-Haenszel odds ratio
		double mh_d=0; // denominator for Mantel-Haenszel odds ratio
		double mh_or=std::numeric_limits<double>::signaling_NaN();
		double mh_ci_n=0;
		double mh_ci_d=0;
		double mh_lci=std::numeric_limits<double>::signaling_NaN();
		double mh_rci=std::numeric_limits<double>::signaling_NaN();
		for (auto &x:SubGrp)
		{
			int table[4];
			table[0] = ( FETrec ? x.second.cgU[0]+x.second.cgU[1] : x.second.cgU[0]					);	// Healthy not-exposed
			table[1] = ( FETrec ? x.second.cgU[2]				  : x.second.cgU[1]+x.second.cgU[2] );	// Healthy Exposed
			table[2] = ( FETrec ? x.second.cgA[0]+x.second.cgA[1] : x.second.cgA[0]					);	// Diseased not-exposed
			table[3] = ( FETrec ? x.second.cgA[2]				  : x.second.cgA[1]+x.second.cgA[2] );	// Diseased Exposed
			if (table[0]+table[1] && table[2]+table[3] && table[0]+table[2] && table[1]+table[3])
			{
				FisherExactTest fet(true);
				fet.input_RxC(table,2,2);
				double pv = one_sided ? fet.result('R') : fet.result('2');
				double ct = table[0]+table[1];
				double cs = table[2]+table[3];
				double n = int(4/(1/cs+1/ct)+0.5);
				// previously double n = ( ct/cs>4 ? cs*5 : cs+ct );
				// previously double n = cs+ct;
				double w = sqrt(n);
				double z = qnorms(1-pv);
				overall_n += w*z;
				overall_d += w*w;
				double a,b,c,d;
				if ((cc&&(table[0]==0||table[1]==0||table[2]==0||table[3]==0)) || constC) {
					a=table[0]+0.5;
					b=table[1]+0.5;
					c=table[2]+0.5;
					d=table[3]+0.5;
					n=n+2; }
				else {
					a=table[0];
					b=table[1];
					c=table[2];
					d=table[3]; }
				double y = b*c/n;
				mh_n += a*d/n;
				mh_d += y;
				double v = 1/a + 1/b + 1/c + 1/d;
				mh_ci_n += y*y*v;
				mh_ci_d += y;
				if (perch::_Debug) lns<<showl<<str_of_uniq(p.gnm,',')<<' '<<x.first<<' '<<pv<<' '<<w<<' '<<z<<' '<<table[0]<<' '<<table[1]<<' '<<table[2]<<' '<<table[3]<<flush_logger;
			}
			else if (perch::_Debug) lns<<showl<<"not_analyze "<<x.first<<' '<<table[0]<<' '<<table[1]<<' '<<table[2]<<' '<<table[3]<<flush_logger;
		}
		if (overall_d)
		{
			overall_z = overall_n / sqrt(overall_d);
			overall_p = 1-cdf_norms(overall_z);
			if (perch::_Debug) lns<<showl<<str_of_uniq(p.gnm,',')<<' '<<overall_n<<' '<<overall_d<<' '<<overall_z<<' '<<1-cdf_norms(overall_z)<<flush_logger;
		}
		if (mh_d)
		{
			mh_or = mh_n / mh_d;
			double SE = sqrt(mh_ci_n/(mh_ci_d*mh_ci_d));
			mh_lci = exp(log(mh_or)-CI*SE);
			mh_rci = exp(log(mh_or)+CI*SE);
		}
		if (!std::isnan(overall_p))
		{
			double mlp = -log10(overall_p);
			p.res = mlp;
			if (AnalysisMethod==FET)
			{
				if		(out_lp) p.out = ftos(mlp							,MisStr);
				else if (out_pv) p.out = ftos(overall_p						,MisStr);
				else if (out_or) p.out = ftos(mh_or							,MisStr);
				else			 p.out = ftos(genepi::P_to_LOD(overall_p)	,MisStr);
			}
			else
			{
				p.out = p.chr_str()+"\t"+itos(p.idx.size())+"\t"+ftos(mlp)+"\t"+ftos(mh_or)+"\t["+ftos(mh_lci)+","+ftos(mh_rci)+"]"+"\t"+str_of_container(cgA,'\t')+"\t"+str_of_container(cgU,'\t');
				// the 4th column of this output is also read by vGrp --detail
			}
		}
		else
		{
			if (AnalysisMethod==FET) p.out=MisStr;
			else p.out = p.chr_str()+"\t"+itos(p.idx.size())+"\tNA\tNA\t[NA,NA]"+"\t"+str_of_container(cgA,'\t')+"\t"+str_of_container(cgU,'\t');
		}
	}
	else if ((AnalysisMethod==CL2 || AnalysisMethod==FET) && FETiwt) // count damaged haplotypes || Fisher's exact test with mid-p method by Irwin's rule
	{
		struct SubGroupData {
			vector<int> cgA;			// count 0/0 0/1 1/1 for affected
			vector<int> cgU;			// count 0/0 0/1 1/1 for unaffected
			int		ss[2];
			double 	table[4];
			SubGroupData() { ss[0]=ss[1]=0; table[0]=table[1]=table[2]=table[3]=0; cgA={0,0,0}; cgU={0,0,0}; }
			void add(double aff, int g, double wt)
			{
				if (aff)
				{
					table[2] += (2-g)* wt;	// Diseased not-exposed
					table[3] +=    g * wt;	// Diseased Exposed
					++cgA[g];
					++ss[1];
				}
				else
				{
					table[0] += (2-g)* wt;	// Healthy not-exposed
					table[1] +=    g * wt;	// Healthy Exposed
					++cgU[g];
					++ss[0];
				}
			}
		};
		map<string,SubGroupData> SubGrp;
		for (size_t i=0;i<newGtM.size();++i)
		{
			if (!GTP::usable(newGtM[i])) continue;
			int g = genotype_code(newGtM[i], newGtC[i], p);
			if (g>=0) SubGrp[strata[i]].add(DepVar[i],g,origWt[i]);
		}
		
		vector<int> cgA = { 0,0,0 };	// count 0/0 0/1 1/1 for affected
		vector<int> cgU = { 0,0,0 };	// count 0/0 0/1 1/1 for unaffected
		double HN = 0; // Healthy not-exposed
		double HE = 0; // Healthy Exposed
		double DN = 0; // Diseased not-exposed
		double DE = 0; // Diseased Exposed
		set<string> to_remove;
		for (auto &x:SubGrp)
		{
			if (x.second.ss[1]<FETmin || x.second.ss[0]<FETmin) {
				if (perch::_Debug) lns<<showl<<"removed "<<x.first<<' '<<x.second.ss[1]<<' '<<x.second.ss[0]<<flush_logger;
				to_remove.insert(x.first);
				continue; }
			HN+=x.second.table[0];
			HE+=x.second.table[1];
			DN+=x.second.table[2];
			DE+=x.second.table[3];
			for (int g=0;g<3;++g) { cgA[g]+=x.second.cgA[g]; cgU[g]+=x.second.cgU[g]; }
		}
		for (auto &x:to_remove) SubGrp.erase(x);
		bool cc = immedC || (delayC && (HN==0||HE==0||DN==0||DE==0));
		
		double overall_n=0; // nominator for overall p-value calculation
		double overall_d=0; // denominator for overall p-value calculation
		double overall_z=std::numeric_limits<double>::signaling_NaN(); // z score for overall p-value calculation
		double overall_p=std::numeric_limits<double>::signaling_NaN();
		double mh_n=0; // nominator for Mantel-Haenszel odds ratio
		double mh_d=0; // denominator for Mantel-Haenszel odds ratio
		double mh_or=std::numeric_limits<double>::signaling_NaN();
		double mh_ci_n=0;
		double mh_ci_d=0;
		double mh_lci=std::numeric_limits<double>::signaling_NaN();
		double mh_rci=std::numeric_limits<double>::signaling_NaN();
		for (auto &x:SubGrp)
		{
			double table[4];
			table[0]=x.second.table[0];
			table[1]=x.second.table[1];
			table[2]=x.second.table[2];
			table[3]=x.second.table[3];
			if (table[0]+table[1] && table[2]+table[3] && table[0]+table[2] && table[1]+table[3])
			{
				double pv = Fishers_exact_test_2x2(table,true,one_sided?'R':'2');
				double ct = x.second.ss[0];
				double cs = x.second.ss[1];
				double n = int(4/(1/cs+1/ct)+0.5);
				// previously double n = ( ct/cs>4 ? cs*5 : cs+ct );
				// previously double n = cs+ct;
				double w = sqrt(n);
				double z = qnorms(1-pv);
				overall_n += w*z;
				overall_d += w*w;
				double a,b,c,d;
				if ((cc&&(table[0]==0||table[1]==0||table[2]==0||table[3]==0)) || constC) {
					a=x.second.table[0]+0.5;
					b=x.second.table[1]+0.5;
					c=x.second.table[2]+0.5;
					d=x.second.table[3]+0.5;
					n=n+2; }
				else {
					a=x.second.table[0];
					b=x.second.table[1];
					c=x.second.table[2];
					d=x.second.table[3]; }
				double y = b*c/n;
				mh_n += a*d/n;
				mh_d += y;
				double v = 1/a + 1/b + 1/c + 1/d;
				mh_ci_n += y*y*v;
				mh_ci_d += y;
				if (perch::_Debug) lns<<showl<<str_of_uniq(p.gnm,',')<<' '<<x.first<<' '<<pv<<' '<<w<<' '<<z<<' '<<table[0]<<' '<<table[1]<<' '<<table[2]<<' '<<table[3]<<flush_logger;
			}
		}
		if (overall_d)
		{
			overall_z = overall_n / sqrt(overall_d);
			overall_p = 1-cdf_norms(overall_z);
			if (perch::_Debug) lns<<showl<<str_of_uniq(p.gnm,',')<<' '<<overall_n<<' '<<overall_d<<' '<<overall_z<<' '<<1-cdf_norms(overall_z)<<flush_logger;
		}
		if (mh_d)
		{
			mh_or = mh_n / mh_d;
			double SE = sqrt(mh_ci_n/(mh_ci_d*mh_ci_d));
			mh_lci = exp(log(mh_or)-CI*SE);
			mh_rci = exp(log(mh_or)+CI*SE);
		}
		if (!std::isnan(overall_p))
		{
			double mlp = -log10(overall_p);
			p.res = mlp;
			if (AnalysisMethod==FET)
			{
				if		(out_lp) p.out = ftos(mlp							,MisStr);
				else if (out_pv) p.out = ftos(overall_p						,MisStr);
				else if (out_or) p.out = ftos(mh_or							,MisStr);
				else			 p.out = ftos(genepi::P_to_LOD(overall_p)	,MisStr);
			}
			else
			{
				p.out = p.chr_str()+"\t"+itos(p.idx.size())+"\t"+ftos(mlp)+"\t"+ftos(mh_or)+"\t["+ftos(mh_lci)+","+ftos(mh_rci)+"]"+"\t"+str_of_container(cgA,'\t')+"\t"+str_of_container(cgU,'\t');
				// prv: ftos(DN)+"\t"+ftos(DE)+"\t.\t"+ftos(HN)+"\t"+ftos(HE)+"\t.";
				// the 4th column of this output is also read by vGrp --detail
			}
		}
		else
		{
			if (AnalysisMethod==FET) p.out=MisStr;
			else p.out = p.chr_str()+"\t"+itos(p.idx.size())+"\tNA\tNA\t[NA,NA]"+"\t"+str_of_container(cgA,'\t')+"\t"+str_of_container(cgU,'\t');
			// prv: ftos(DN)+"\t"+ftos(DE)+"\t.\t"+ftos(HN)+"\t"+ftos(HE)+"\t.";
		}
	}
	else if (AnalysisMethod==RNK)
	{
		struct SubGroupData {
			set<double> values;
			set<int>	rnkdep;
			deque< pair<double,int> > kwdata;
			SubGroupData() {}
		};
		map<string,SubGroupData> SubGrp;
		for (size_t i=0;i<newGtM.size();++i)
		{
			if (!GTP::usable(newGtM[i])) continue;
			double v=risk_score(newGtM[i], newGtC[i], p);
			int group = DepVar[i]?2:1;
			SubGrp[strata[i]].kwdata.push_back(make_pair(v,group));
			SubGrp[strata[i]].values.insert(v);
			SubGrp[strata[i]].rnkdep.insert(group);
		}
		
		set<string> to_remove;
		for (auto &x:SubGrp)
		{
			if ((int)x.second.kwdata.size()<RNKmin || x.second.values.size()<2 || x.second.rnkdep.size()<2) {
				if (perch::_Debug) lns<<showl<<"removed "<<x.first<<' '<<x.second.kwdata.size()<<' '<<x.second.values.size()<<' '<<x.second.rnkdep.size()<<flush_logger;
				to_remove.insert(x.first);
				continue; }
		}
		for (auto &x:to_remove) SubGrp.erase(x);
		
		double overall_n=0; // nominator for overall p-value calculation
		double overall_d=0; // denominator for overall p-value calculation
		double overall_z=std::numeric_limits<double>::signaling_NaN(); // z score for overall p-value calculation
		double overall_p=std::numeric_limits<double>::signaling_NaN();
		for (auto &x:SubGrp)
		{
			double pv = ranksum(x.second.kwdata,one_sided?'<':'=');
			if (!std::isnan(pv))
			{
				double w = sqrt(x.second.kwdata.size());
				double z = qnorms(1-pv);
				overall_n += w*z;
				overall_d += w*w;
			}
		}
		if (overall_d)
		{
			overall_z = overall_n / sqrt(overall_d);
			overall_p = 1-cdf_norms(overall_z);
			if (perch::_Debug) lns<<showl<<str_of_uniq(p.gnm,',')<<' '<<overall_n<<' '<<overall_d<<' '<<overall_z<<' '<<1-cdf_norms(overall_z)<<flush_logger;
		}
		if (!std::isnan(overall_p))
		{
			double mlp = -log10(overall_p);
			p.res = mlp;
			if		(out_lp) p.out = ftos(mlp							,MisStr);
			else if (out_pv) p.out = ftos(overall_p						,MisStr);
			else			 p.out = ftos(genepi::P_to_LOD(overall_p)	,MisStr);
		}
		else
		{
			p.out=MisStr;
		}
	}
	else if (AnalysisMethod==REG)
	{
		double t_cs=0, t_ct=0;	// total number of samples
		for (size_t i=0;i<newGtM.size();++i)
		{
			if (!GTP::usable(newGtM[i])) continue;
			if (DepVar[i]) ++t_cs; else ++t_ct;
		}
		
		if (!t_cs) // rm !t_ct for linear regression
		{
			p.res = std::numeric_limits<double>::signaling_NaN();
			p.out = MisStr;
		}
		else
		{
			int				np=2+CovVar[0].size();	// number of parameters (1 constance, 1 CollapsedGenotype, # covariates)
			int				ni=t_cs+t_ct;			// number of individuals
			Eigen::VectorXd y(ni);					// y
			Eigen::MatrixXd X(ni,np);				// X
			
			for (size_t i=0,j=0;i<newGtM.size();++i)
			{
				if (!GTP::usable(newGtM[i])) continue;
				double v=risk_score(newGtM[i], newGtC[i], p);
				y(j)  =sqrtWt[i]*DepVar[i];
				X(j,0)=sqrtWt[i];
				X(j,1)=sqrtWt[i]*(v);
				//old:
				//y(j)  =DepVar[i];
				//X(j,0)=1;
				//X(j,1)=v;
				for (size_t c=0;c<CovVar[i].size();++c) X(j,c+2)=CovVar[i][c];
				++j;
			}
			if (!prData.empty()) // for debug only
			{
				openOutFile_or_exit(out,prData+"."+str_of_uniq(p.gnm,','));
				for (int k=0;k<ni;++k)
				{
					out << y(k) ;
					for (int l=0;l<np;++l)
						out << '\t' << X(k,l) ;
					out << endl;
				}
				closefile (out);
			}
			double res = regress(X,y,one_sided); // could be PV or OR depending on out_or
			double mlp = 0;
			if (out_or) { p.res = res; }					// res is OR
			else		{ mlp = -log10(res); p.res = mlp; } // res is PV
			if		(out_lp) p.out = ftos(mlp					,MisStr);
			else if (out_pv) p.out = ftos(res					,MisStr);
			else if (out_or) p.out = ftos(res					,MisStr);
			else			 p.out = ftos(genepi::P_to_LOD(res)	,MisStr);
		}
	}
	else if (AnalysisMethod==FLR)
	{
		double t_cs=0, t_ct=0;	// total number of samples
		for (size_t i=0;i<newGtM.size();++i)
		{
			if (!GTP::usable(newGtM[i])) continue;
			if (DepVar[i]) ++t_cs; else ++t_ct;
		}
		
		if (!t_cs || !t_ct)
		{
			p.res = std::numeric_limits<double>::signaling_NaN();
			p.out = MisStr;
		}
		else
		{
			string group=uID; // str_of_uniq(p.gnm,',') may be too long for GSA
			
			// write data file
			set<double> scores;
			openOutFile_or_exit(data,perch::TMPDIR+"firth_"+group+".data");
			data<<"SeqID\tAff";
			for (size_t c=0;c<CovVar[0].size();++c) data<<"\t"<<CovNID[c];
			data<<"\tgene\tweight"<<endl;
			for (size_t i=0;i<newGtM.size();++i)
			{
				if (!GTP::usable(newGtM[i])) continue;
				double v=risk_score(newGtM[i], newGtC[i], p);
				if (std::isnan(v)) continue;
				scores.insert(v);
				data<<record[i]<<'\t'<<DepVar[i];
				for (size_t c=0;c<CovVar[i].size();++c) data<<'\t'<<CovVar[i][c];
				data<<'\t'<<v<<'\t'<<sqrtWt[i]<<endl;
			}
			closefile(data);
			
			if (scores.size()<2)
			{
				exec("rm \""+perch::TMPDIR+"firth_"+group+".data\"",false);
			}
			else
			{
				// write R script
				openOutFile_or_exit(R,perch::TMPDIR+"firth_"+group+".R");
				R<<"suppressWarnings(library(logistf))"<<endl;
				if (p.MHC) R<<"involve_MHC <- TRUE"<<endl; else R<<"involve_MHC <- FALSE"<<endl;
				R<<"x<-read.table(\""+perch::TMPDIR+"firth_"<<group<<".data\", header=TRUE, colClasses=c(\"character\",\"numeric\"";
				for (size_t c=0;c<CovVar[0].size();++c) R<<",\"numeric\"";
				R<<",\"numeric\",\"numeric\"))"<<endl;
				if (MyCode.empty())
				{
					R<<"fit_firth<-logistf(Aff ~ gene";
					for (size_t c=0;c<CovVar[0].size();++c) R<<"+"<<CovNID[c];
					R<<", data=x, weights=weight"<<(flr_lr?"":", pl=FALSE")<<")"<<endl;
					R<<"summary(fit_firth)"<<endl;
					R<<"pv=fit_firth$prob[2]"<<endl;
					R<<"cat(\"MyR_pval \", pv, \"\\n\", sep = \"\")"<<endl;
				}
				else
				{
					R<<MyCode;
				}
				closefile(R);
				if (perch::_Debug) cerr<<"intermediate files are \""+perch::TMPDIR+"firth_"+group+".*\""<<endl;
				
				// execute
				string result=exec("Rscript \""+perch::TMPDIR+"firth_"+group+".R\"",false);
				
				// read results
				double coef = std::numeric_limits<double>::signaling_NaN();
				double secoef = std::numeric_limits<double>::signaling_NaN();
				double oddsr = std::numeric_limits<double>::signaling_NaN();
				double lower = std::numeric_limits<double>::signaling_NaN();
				double upper = std::numeric_limits<double>::signaling_NaN();
				{
					vector<string> lines;
					boost::split(lines,result,boost::is_any_of("\n"));
					vector<string> words;
					for (auto &l:lines)
						if (str_startsw(l,"gene")) { boost::split(words,l,boost::is_any_of(" "),boost::token_compress_on); break; }
					if (words.size()<5) exit_error("cannot read logistf results for "+group+"\nResult: "+result);
					if (!read_val_noNaN(words[1],coef))   exit_error("cannot read logistf coef "+words[1]+" for "+group);
					if (!read_val_noNaN(words[2],secoef)) exit_error("cannot read logistf se(coef) "+words[2]+" for "+group);
					if (!read_val_noNaN(words[3],lower))  exit_error("cannot read logistf lower95 "+words[3]+" for "+group);
					if (!read_val_noNaN(words[4],upper))  exit_error("cannot read logistf upper95 "+words[4]+" for "+group);
					oddsr=exp(coef);
					lower=exp(lower);
					upper=exp(upper);
					if (out_se)
					{
						if (!std::isnan(oddsr)&&!std::isnan(lower)&&!std::isnan(upper)&&!std::isnan(coef)&&!std::isnan(secoef))
							p.xtr="OR="+ftos(oddsr)+"["+ftos(lower)+","+ftos(upper)+"]|coef="+ftos(coef)+"|se_coef="+ftos(secoef);
						else p.xtr=MisStr;
					}
					else
					{
						if (!std::isnan(oddsr)&&!std::isnan(lower)&&!std::isnan(upper))
							p.xtr="OR="+ftos(oddsr)+"["+ftos(lower)+","+ftos(upper)+"]";
						else p.xtr=MisStr;
					}
				}
				double pv = std::numeric_limits<double>::signaling_NaN();
				double mlp = std::numeric_limits<double>::signaling_NaN();
				{
					bool read_pval=false;
					vector<string> lines;
					boost::split(lines,result,boost::is_any_of("\n"));
					for (auto &l:lines)
					{
						if (str_startsw(l,"MyR_pval"))
						{
							vector<string> words;
							boost::split(words,l,boost::is_any_of(" "),boost::token_compress_on);
							if (words.size()<2) exit_error("cannot read MyR_pval for "+str_of_uniq(p.gnm,',')+"\nResult: "+result);
							if (!read_val(words[1],pv)) exit_error("cannot read MyR_pval "+words[1]+" for "+str_of_uniq(p.gnm,','));
							if (!std::isnan(pv))
							{
								if (one_sided)
								{
									if (oddsr>1.0)	pv = pv/2;
									else			pv = 1-pv/2;
								}
								mlp = -log10(pv);
							}
							read_pval=true;
						}
					}
					if (!read_pval) exit_error("there is no MyR_pval output from the R code");
				}
				if (out_or) { p.res = oddsr; }
				else		{ p.res = mlp; }
				if		(out_lp) p.out = ftos(mlp					,MisStr);
				else if (out_pv) p.out = ftos(pv					,MisStr);
				else if (out_or) p.out = ftos(oddsr					,MisStr);
				else			 p.out = ftos(genepi::P_to_LOD(pv)	,MisStr);
				if (!perch::_Debug) exec("rm \""+perch::TMPDIR+"firth_"+group+".R\"",false);
				if (!perch::_Debug) exec("rm \""+perch::TMPDIR+"firth_"+group+".data\"",false);
			}
		}
	}
	else if (AnalysisMethod==MyR)
	{
		double t_cs=0, t_ct=0;	// total number of samples
		for (size_t i=0;i<newGtM.size();++i)
		{
			if (!GTP::usable(newGtM[i])) continue;
			if (DepVar[i]) ++t_cs; else ++t_ct;
		}
		
		if (!t_cs || !t_ct)
		{
			p.res = std::numeric_limits<double>::signaling_NaN();
			p.out = MisStr;
		}
		else
		{
			string group=uID; // str_of_uniq(p.gnm,',') may be too long for GSA
			
			// write data file
			set<double> scores;
			openOutFile_or_exit(data,perch::TMPDIR+"perch_MyR_"+group+".data");
			data<<"SeqID\tAff";
			for (size_t c=0;c<CovVar[0].size();++c) data<<"\t"<<CovNID[c];
			data<<"\tgene\tweight"<<endl;
			for (size_t i=0;i<newGtM.size();++i)
			{
				if (!GTP::usable(newGtM[i])) continue;
				double v=risk_score(newGtM[i], newGtC[i], p);
				if (std::isnan(v)) continue;
				scores.insert(v);
				data<<record[i]<<'\t'<<DepVar[i];
				for (size_t c=0;c<CovVar[i].size();++c) data<<'\t'<<CovVar[i][c];
				data<<'\t'<<v<<'\t'<<sqrtWt[i]<<endl;
			}
			closefile(data);
			
			if (scores.size()<2)
			{
				exec("rm \""+perch::TMPDIR+"perch_MyR_"+group+".data\"",false);
			}
			else
			{
				// write R script and execute
				openOutFile_or_exit(R,perch::TMPDIR+"perch_MyR_"+group+".R");
				if (p.MHC) R<<"involve_MHC <- TRUE"<<endl; else R<<"involve_MHC <- FALSE"<<endl;
				R<<"x<-read.table(\""+perch::TMPDIR+"perch_MyR_"<<group<<".data\", header=TRUE, colClasses=c(\"character\",\"numeric\"";
				for (size_t c=0;c<CovVar[0].size();++c) R<<",\"numeric\"";
				R<<",\"numeric\",\"numeric\"))"<<endl;
				R<<MyCode;
				closefile(R);
				string result=exec("Rscript \""+perch::TMPDIR+"perch_MyR_"+group+".R\"",false);
				
				// read results
				double pv = std::numeric_limits<double>::signaling_NaN();
				double mlp = std::numeric_limits<double>::signaling_NaN();
				{
					bool read_pval=false;
					vector<string> lines;
					boost::split(lines,result,boost::is_any_of("\n"));
					for (auto &l:lines)
					{
						if (str_startsw(l,"MyR_pval"))
						{
							vector<string> words;
							boost::split(words,l,boost::is_any_of(" "),boost::token_compress_on);
							if (words.size()<2) exit_error("cannot read MyR_pval for "+str_of_uniq(p.gnm,',')+"\nResult: "+result);
							if (!read_val(words[1],pv)) exit_error("cannot read MyR_pval "+words[1]+" for "+str_of_uniq(p.gnm,','));
							if (!std::isnan(pv))
							{
								mlp = -log10(pv);
								p.res = mlp;
							}
							read_pval=true;
						}
						if (str_startsw(l,"MyR_info"))
						{
							vector<string> words;
							boost::split(words,l,boost::is_any_of(" "),boost::token_compress_on);
							if (words.size()<2) exit_error("cannot read MyR_info for "+str_of_uniq(p.gnm,',')+"\nResult: "+result);
							p.xtr=words[1];
						}
					}
					if (!read_pval) exit_error("there is no MyR_pval output from the R code");
				}
				if (out_or) exit_error("--my-R doesn't support --out-or");
				if		(out_lp) p.out = ftos(mlp					,MisStr);
				else if (out_pv) p.out = ftos(pv					,MisStr);
				else			 p.out = ftos(genepi::P_to_LOD(pv)	,MisStr);
				if (!perch::_Debug) exec("rm \""+perch::TMPDIR+"perch_MyR_"+group+".R\"",false);
				if (!perch::_Debug) exec("rm \""+perch::TMPDIR+"perch_MyR_"+group+".data\"",false);
				if (perch::_Debug) cerr<<"intermediate files are \""+perch::TMPDIR+"perch_MyR_"+group+".*\""<<endl;
			}
		}
	}
	else if (AnalysisMethod==OUT)
	{
		double t_cs=0, t_ct=0;	// total number of samples
		for (size_t i=0;i<newGtM.size();++i)
		{
			if (!GTP::usable(newGtM[i])) continue;
			if (DepVar[i]) ++t_cs; else ++t_ct;
		}
		
		/*/ if (!t_cs)
		{
			;
		}
		else //*/
		{
			if (AllDataGtp.empty()) AllDataGtp.assign(record.size(),vector<double>());
			AllDataSym.push_back(p.grp); // AllDataSym.push_back(str_of_uniq(p.gnm,','));
			for (size_t i=0;i<newGtM.size();++i)
			{
				if (!GTP::usable(newGtM[i])) { AllDataGtp[i].push_back(std::numeric_limits<double>::signaling_NaN()); continue; }
				double v=risk_score(newGtM[i], newGtC[i], p);
				AllDataGtp[i].push_back(v);
			}
			p.out="SKIPPED=NoCalc";
		}
	}
	else if (AnalysisMethod==SSU)
	{
		double t_cs=0, t_ct=0;	// total number of samples
		for (size_t i=0;i<newGtM.size();++i)
		{
			if (!GTP::usable(newGtM[i])) continue;
			if (DepVar[i]) ++t_cs; else ++t_ct;
		}
		
		if (!t_cs || !t_ct)
		{
			p.res = std::numeric_limits<double>::signaling_NaN();
			p.out = MisStr;
		}
		else
		{
			int				np = newGtM[0].size()+CovVar[0].size();	// number of parameters
			int				ni = t_cs+t_ct;		// number of individuals
			Eigen::VectorXd y(ni);				// y
			Eigen::MatrixXd X(ni,np);			// X
			vector<double>  W = p.wts;
			W.resize(np,*std::max_element(p.wts.begin(),p.wts.end()));
			
			for (size_t i=0,j=0;i<newGtM.size();++i)
			{
				if (!GTP::usable(newGtM[i])) continue;
				size_t c=0,l=0;
				y(j)=DepVar[i];
				for (l=0;l<newGtM[i].size();++l) X(j,l)  =GTP::num_alt_recessive(newGtM[i][l]);
				for (c=0;c<CovVar[i].size();++c) X(j,l+c)=CovVar[i][c];
				++j;
			}
			double pv = do_SSUw(X,y,W);
			double mlp = -log10(pv);
			p.res = mlp;
			if		(out_lp) p.out = ftos(mlp					,MisStr);
			else if (out_pv) p.out = ftos(pv					,MisStr);
			else			 p.out = ftos(genepi::P_to_LOD(pv)	,MisStr);
		}
	}
	else if (AnalysisMethod==DET) ;
	else exit_error("Wrong AnalysisMethod.");
	
	if (!p.pID && program.nt>1)
	{
		bjm_res_db_mutex.lock();
		bjm_res_db[overall_idx]=make_tuple(p.out,p.xtr,p.res);
		bjm_res_db_mutex.unlock();
	}
	else if (p.pID && !std::isnan(p.res))
	{
		wjp_res_db_mutex.lock();
		wjp_res_db[p.pID].push_back(p.res);
		wjp_res_db_mutex.unlock();
	}
}

void write_results(bool wrINFO, field_numbers& FldInf, field_numbers& FldRes, const string& wrWhat, bool print_last)
{
	if (!permte && program.nt>1)
	{
		bjmJobs.push_back(DefaultJob);
	}
	else if (permte)
	{
		// calculate
		JobData& p = DefaultJob;
		MtJobs<JobData>	wjpJobs;
		for (size_t i=1;i<PmtPtr.size();++i)
		{
			wjpJobs.push_back(p);
			wjpJobs.back().pID=i;
		}
		wjpJobs.run(_run_job,program.nt);
	}
	else
	{
		// calculate
		JobData& p = DefaultJob;
		_run_job(p);
		if (!wjp_res_mx.empty() && !std::isnan(p.res))
		{
			double pv = (wjp_res_mx.num_ge(p.res)+2.0) / (wjp_res_mx.size()+4.0); // make sure pv!=0 and pv!=1
			if		(out_lp) p.res = -log10(pv);
			else if (out_pv) p.res = pv;
			else			 p.res = genepi::P_to_LOD(pv);
		}

		// print results
		// if (!str_startsw(DefaultJob.out,"SKIPPED"))
		{
			if (AnalysisMethod==CL1 || AnalysisMethod==CL2) // print 1 line per gene
			{
				if (!str_startsw(DefaultJob.out,"SKIPPED"))
					program.outf << str_of_uniq(DefaultJob.gnm,',') << DLMTR << DefaultJob.out << endl;
			}
			else if (AnalysisMethod==DET && (showGT||showIF)) // print 1 line per genotype
			{
				int i=0;
				if (print_last)
				{
					for (each_element(program.main_data(),it1),++i)
						if (!str_startsw(DefaultJob.qcl[i],"SKIPPED")) program.outf<<DefaultJob.qcl[i];
				}
				else
				{
					for (each_interval(program.main_data(),it1,it2),++i)
						if (!str_startsw(DefaultJob.qcl[i],"SKIPPED")) program.outf<<DefaultJob.qcl[i];
				}
			}
			else if (AnalysisMethod!=OUT) // print 1 line per variant
			{
				int i=0;
				if (print_last)
				{
					for (each_element(program.main_data(),it1),++i)
					{
						string output;
						if (AnalysisMethod==DET)	output = DefaultJob.qcl[i];
						else						output = DefaultJob.out;
						if (wrINFO) (*it1)[FldRes[0]] += ";" + wrWhat + "=" + output;
						else		(*it1)[FldRes[0]]  = output;
						if (!DefaultJob.xtr.empty()) (*it1)[FldInf[0]]+=";"+perch::h_Axtr+"="+DefaultJob.xtr;
						print_container(*it1,program.outf,DLMTR,true);
					}
				}
				else
				{
					for (each_interval(program.main_data(),it1,it2),++i)
					{
						string output;
						if (AnalysisMethod==DET)	output = DefaultJob.qcl[i];
						else						output = DefaultJob.out;
						if (wrINFO) (*it1)[FldRes[0]] += ";" + wrWhat + "=" + output;
						else		(*it1)[FldRes[0]]  = output;
						if (!DefaultJob.xtr.empty()) (*it1)[FldInf[0]]+=";"+perch::h_Axtr+"="+DefaultJob.xtr;
						print_container(*it1,program.outf,DLMTR,true);
					}
				}
			}
		}
	}
	if (print_last) program.main_data().clear();
	else			program.main_data().clear_ExceptTheLastRow();
}

char _read_trio(const perch::trio& t, const vector<string>& in, int chr_num, int bp, bool DNonly, int& num_trio_has_de_novo, const genepi::gtp_par& gpar)
{
	if (t.col[1] && t.col[2])
	{
		if (t.col[0]>=(int)in.size() || t.col[1]>=(int)in.size() || t.col[2]>=(int)in.size())
		{
			lns<<showe<<"input file doesn't havce sufficient fields: "<<t.col[0]<<' '<<t.col[1]<<' '<<t.col[2]<<' '<<in.size()<<' '<<t.id[0]<<' '<<t.id[1]<<' '<<t.id[2]<<' '<<chr_num<<' '<<bp<<fatal;
		}
		char g_ch = GTP::read(in[t.col[0]],chr_num,bp,t.sex[0],gpar);
		char g_pa = GTP::read(in[t.col[1]],chr_num,bp,t.sex[1],gpar);
		char g_ma = GTP::read(in[t.col[2]],chr_num,bp,t.sex[2],gpar);
		bool is_dn = GTP::denovo(g_ch,g_pa,g_ma);
		if (is_dn) ++num_trio_has_de_novo;
		if (DNonly && !is_dn && GTP::usable(g_ch)) GTP::set_all_ref(g_ch);
		return g_ch;
	}
	else
	{
		if (t.col[0]>=(int)in.size())
		{
			lns<<showe<<"input file doesn't havce sufficient fields: "<<t.col[0]<<' '<<in.size()<<' '<<t.id[0]<<' '<<chr_num<<' '<<bp<<fatal;
		}
		char g_ch = GTP::read(in[t.col[0]],chr_num,bp,t.sex[0],gpar);
		return g_ch;
	}
}

int main (int argc, char * const argv[])
{
	struct ChrReg {
		vector<string>		in; // input filenames
		int					bp; // padding basepair
		genepi::ChrRegions	cr; // chromosome regions
	};

	// basic parameters
	perch::prelude(argc,argv);
	
	// other parameters
	set<string>		aLBF_h={"gene","trait","gene\\trait","trait\\gene","symbol\\trait","trait\\symbol","symbol","gene_symbol","genesymbol"}; // header of column 1 of aLBF_f
	string			spl_in;						// input sample file with header: SeqID Sex Aff/QT Cov1 Cov2 ..
	string			gtp_in;						// input genotype file: VCF 1-8 columns + genotypes
	string			var_in;						// input variant file: VCF 1-8 columns + annotations
	string			ex_chr = "1";				// example chromosome string (1)
	vector<string>	sv_inp;						// input structual variation data
	string			ped_in;						// input pedigree file in LINKAGE format
	string			wrWhat = perch::h_HLRb;		// header of output
	string			h_dels;						// header of BayesDel defined by user
	set<string>		h_grID = {perch::h_symb};	// header of group IDs
	double			MisPVl = 0;					// skip variants if missing % between cases and controls has p-value < MisPVl, 0 = noQC. prev 0.000001
	bool			AddInf = false;				// write result to column 8 (INFO)
	bool			Wt_DEL = true;				// weight variants by deleteriousness
	bool			Wt_VQS = true;				// weight variants by variant call quality
	bool			Wt_GBA = false;				// weight variants by gene relevance
	bool			NegGBA = true;				// add negative GBA
	bool			ToWait = false;				// wait for vSIM to print the first comment line of the genotype file
	bool			Adaptv = false;				// up-lift filSAF. def=false 1) not necessary as MAf is from combined samples, 2) may include too many var, 3) is arbitrary.
	bool			Var_Wt = false;				// apply variant weights, otherwise all filtered variants have a weight 1, i.e., probability of damaging is 100%
	bool			single = false;				// single-variant-based analysis
	bool			skipnm = false;				// skip the variant in the group file that has no match in VCF
	bool			PRSPAF = false;				// PRS calibrate by population allele frequency calcualted from data
	bool			FltLinkageErr = true;
	double			cap_pg = 0;					// put a cap on the total number of variants per gene.   Not helpful by simulation.
	bool			cap_ps = false;				// put a cap on the total number of variants per sample. Not helpful by simulation.
	bool			DNonly = false;				// count only de novo mutations for trios, but read genotypes as usual for independent cs & ct. Test is valid if you have both cs & ct trios and no non-trios.
	double			PGuard = 0.5;				// phenocopy guard. in case penetr[0] is too low that does not allow for phenocopies
	int				DN_max = 1;					// the maximum number of the same de novo mutation, for QC
	bool			HasAll = false;				// for permutation: the input genotype file has all variants (no splitting by chr)
	vector<string>	MxVlFs;						// for permutation: the files storing the max values from each repeat
	int				bigspl = 500;				// definition of big sample size. Beyond this, population frequency will be calculated from data (robust to LD), otherwise from MaxAF (not robust to LD)
	int				mac_ge = 0;					// ignore variants with MAC >= mac_ge
	double			fil_sv = 0.001;				// filter gene-wise deletions by frequency in samples, ie, fil_sv * (cs+ct)
	string			GS2GBA;						// GeneSymbol->GBAlbf file
	double			skip_if_maf_l=0;			// common variant only, useful for QC with --single
	string			iwt_fn;
	bool			do_nothing = false;
	string			MyR_in;

	// internal parameters
	double			internal_MisCut=0.5;
	
#ifdef PERCH_DISTRIBUTION
#else
	program.enable_option("--nt");
#endif

	// handle program arguments
	perch::MisCut=1;
	perch::VQSsnv=-INFINITY;
	perch::VQSidl=-INFINITY;
	perch::filflt.clear();
	program.read_arguments(argc,argv);
	perch::read_arguments();
	for (size_t argi=1; argi<program.arg().size(); ++argi)
	{
		if		(str_startsw(program.arg()[argi],"--group"))			ReadSet(program.arg(),argi,h_grID);
		else if (str_startsw(program.arg()[argi],"--spl"))				ReadArg(program.arg(),argi,spl_in);
		else if (str_startsw(program.arg()[argi],"--ped"))				ReadArg(program.arg(),argi,ped_in);
		else if (str_startsw(program.arg()[argi],"--gtp"))				ReadArg(program.arg(),argi,gtp_in);
		else if (str_startsw(program.arg()[argi],"--example-chr"))		ReadArg(program.arg(),argi,ex_chr);
		else if (str_startsw(program.arg()[argi],"--filt-miss-pv"))		ReadArg(program.arg(),argi,MisPVl);
		else if (str_startsw(program.arg()[argi],"--biol"))				ReadArg(program.arg(),argi,GS2GBA);
		else if (str_startsw(program.arg()[argi],"--default-maf"))		ReadArg(program.arg(),argi,defMAF);
		else if (str_startsw(program.arg()[argi],"--big"))				ReadArg(program.arg(),argi,bigspl);
		else if (str_startsw(program.arg()[argi],"--permute"))			ReadArg(program.arg(),argi,permte);
		else if (str_startsw(program.arg()[argi],"--one-file"))			ReadArg(program.arg(),argi,HasAll);
		else if (str_startsw(program.arg()[argi],"--perm-data"))		ReadSet(program.arg(),argi,MxVlFs);
		else if (str_startsw(program.arg()[argi],"--out-mlp"))			ReadArg(program.arg(),argi,out_lp);
		else if (str_startsw(program.arg()[argi],"--out-pv"))			ReadArg(program.arg(),argi,out_pv);
		else if (str_startsw(program.arg()[argi],"--out-or"))			ReadArg(program.arg(),argi,out_or);
		else if (str_startsw(program.arg()[argi],"--out-se"))			ReadArg(program.arg(),argi,out_se);
		else if (str_startsw(program.arg()[argi],"--one-sided"))		ReadArg(program.arg(),argi,one_sided);
		else if (str_startsw(program.arg()[argi],"--print"))			ReadArg(program.arg(),argi,prData);
		else if (str_startsw(program.arg()[argi],"--snv-only"))			ReadArg(program.arg(),argi,AAA::SNonly);
		else if (str_startsw(program.arg()[argi],"--neg-biol"))			ReadArg(program.arg(),argi,NegGBA);
		else if (str_startsw(program.arg()[argi],"--add-info"))			ReadArg(program.arg(),argi,AddInf);
		else if (str_startsw(program.arg()[argi],"--adj-prop"))			ReadArg(program.arg(),argi,AdjPrp);
		else if (str_startsw(program.arg()[argi],"--wait"))				ReadArg(program.arg(),argi,ToWait);
		else if (str_startsw(program.arg()[argi],"--wt-del"))			ReadArg(program.arg(),argi,Wt_DEL);
		else if (str_startsw(program.arg()[argi],"--wt-vqs"))			ReadArg(program.arg(),argi,Wt_VQS);
		else if (str_startsw(program.arg()[argi],"--single"))			ReadArg(program.arg(),argi,single);
		else if (str_startsw(program.arg()[argi],"--common"))			ReadArg(program.arg(),argi,skip_if_maf_l);
		else if (str_startsw(program.arg()[argi],"--cap"))				ReadArg(program.arg(),argi,cap_pg);
		else if (str_startsw(program.arg()[argi],"--top"))				ReadArg(program.arg(),argi,cap_ps);
		else if (str_startsw(program.arg()[argi],"--phenocopy"))		ReadArg(program.arg(),argi,PGuard);
		else if (str_startsw(program.arg()[argi],"--var-wt"))			ReadArg(program.arg(),argi,Var_Wt);
		else if (str_startsw(program.arg()[argi],"--show-id"))			ReadArg(program.arg(),argi,showID);
		else if (str_startsw(program.arg()[argi],"--show-gt"))			ReadArg(program.arg(),argi,showGT);
		else if (str_startsw(program.arg()[argi],"--show-if"))			ReadArg(program.arg(),argi,showIF);
		else if (str_startsw(program.arg()[argi],"--de-novo"))			ReadArg(program.arg(),argi,denovo);
		else if (str_startsw(program.arg()[argi],"--dn-only"))			ReadArg(program.arg(),argi,DNonly);
		else if (str_startsw(program.arg()[argi],"--filt-dn"))			ReadArg(program.arg(),argi,DN_max);
		else if (str_startsw(program.arg()[argi],"--sv"))				ReadSet(program.arg(),argi,sv_inp);
		else if (str_startsw(program.arg()[argi],"--filt-sv"))			ReadArg(program.arg(),argi,fil_sv);
		else if (str_startsw(program.arg()[argi],"--rm-mac-ge"))		ReadArg(program.arg(),argi,mac_ge);
		else if (str_startsw(program.arg()[argi],"--xct-pfx"))			ReadSet(program.arg(),argi,ExtCT::ExAC_pfx); // previously --pub
		else if (str_startsw(program.arg()[argi],"--xct-pop"))			ReadArg(program.arg(),argi,ExtCT::ExAC_pop);
		else if (str_startsw(program.arg()[argi],"--xct-padding"))		ReadArg(program.arg(),argi,ExtCT::ExAC_pad);
		else if (str_startsw(program.arg()[argi],"--xct-cov-pc"))		ReadArg(program.arg(),argi,ExtCT::ExAC_pCv);
		else if (str_startsw(program.arg()[argi],"--xct-case-only"))	ReadArg(program.arg(),argi,ExtCT::ExAC_cso);
		else if (str_startsw(program.arg()[argi],"--xct-frq-pv"))		ReadArg(program.arg(),argi,ExtCT::ExAC_fqc);
		else if (str_startsw(program.arg()[argi],"--xct-log"))			ReadArg(program.arg(),argi,ExtCT::ExAC_log);
		else if (str_startsw(program.arg()[argi],"--weight"))			ReadArg(program.arg(),argi,iwt_fn);
		else if (str_startsw(program.arg()[argi],"--filt-linkage-err"))	ReadArg(program.arg(),argi,FltLinkageErr);
		else if (str_startsw(program.arg()[argi],"--internal-MisCut"))	ReadArg(program.arg(),argi,internal_MisCut);
		else if (str_startsw(program.arg()[argi],"--min-sample"))		ReadArg(program.arg(),argi,RNKmin);
		else if (str_startsw(program.arg()[argi],"--min-cohort"))		ReadArg(program.arg(),argi,FETmin);
		else if (str_startsw(program.arg()[argi],"--recessive"))		ReadArg(program.arg(),argi,FETrec);
		else if (str_startsw(program.arg()[argi],"--trend"))			ReadArg(program.arg(),argi,h_test);
		else if (str_startsw(program.arg()[argi],"--prs"))				ReadArg(program.arg(),argi,usePRS);
		else if (str_startsw(program.arg()[argi],"--por"))				ReadArg(program.arg(),argi,usePOR);
		else if (str_startsw(program.arg()[argi],"--ok-no-match"))		ReadArg(program.arg(),argi,skipnm);
		else if (str_startsw(program.arg()[argi],"--cal-paf"))			ReadArg(program.arg(),argi,PRSPAF);
		else if (str_startsw(program.arg()[argi],"--do-nothing"))		ReadArg(program.arg(),argi,do_nothing);
		else if (str_startsw(program.arg()[argi],"--my-R"))			{	ReadArg(program.arg(),argi,MyR_in); AnalysisMethod=MyR; }
		else if (str_startsw(program.arg()[argi],"--flr-R"))		{	ReadArg(program.arg(),argi,MyR_in); AnalysisMethod=FLR; }
		else if (program.arg()[argi]=="--no-del")	{	Wt_DEL=false; perch::clear_fltdel(false); }
		else if (program.arg()[argi]=="--no-vqs")	{	Wt_VQS=false; }
		else if (program.arg()[argi]=="--detail")	{	AnalysisMethod=DET; }
		else if (program.arg()[argi]=="--count-gtp"){	AnalysisMethod=CL1; }
		else if (program.arg()[argi]=="--count-hap"){	AnalysisMethod=CL2; Var_Wt=false; }
		else if (program.arg()[argi]=="--collapse")	{	AnalysisMethod=CL2; Var_Wt=false; }
		else if (program.arg()[argi]=="--ranksum")	{	AnalysisMethod=RNK; }
		else if (program.arg()[argi]=="--FET")		{	AnalysisMethod=FET; Var_Wt=false; }
		else if (program.arg()[argi]=="--opt-pen")	{	AnalysisMethod=OPN; penetrance = &penetrance_by_gtp; }
		else if (program.arg()[argi]=="--HLR")		{	AnalysisMethod=HLR; penetrance = &penetrance_by_gtp; }
		else if (program.arg()[argi]=="--SSUw")		{	AnalysisMethod=SSU; do_SSUw=&pv_SSUw; }
		else if (program.arg()[argi]=="--SSUc")		{	AnalysisMethod=SSU; do_SSUw=&pv_SSUc; }
		else if (program.arg()[argi]=="--logistic")	{	AnalysisMethod=FLR; }
		else if (program.arg()[argi]=="--linear")	{	AnalysisMethod=REG; regress=&pv_1st_linear;   }
		else if (program.arg()[argi]=="--write")	{	AnalysisMethod=OUT; }
		else if (program.arg()[argi]=="--immed-cor"){	immedC=true;  delayC=false; constC=false; }
		else if (program.arg()[argi]=="--delay-cor"){	immedC=false; delayC=true;  constC=false; }
		else if (program.arg()[argi]=="--const-cor"){	immedC=false; delayC=false; constC=true;  }
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else if (var_in.empty()) var_in=program.arg()[argi];
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help
	program.help_text_var("_Default_sample_file",spl_in);
	program.help_text_var("_Default_groups",str_of_container(h_grID,',',false));
	program.help_text_var("_Default_default_maf",ftos(defMAF));
	program.help_text_var("_Default_filt_miss_d",ftos(MisPVl));
	program.help_text_var("_Default_filt_sv",ftos(fil_sv));
	program.help_text_var("_Default_snv_only",str_YesOrNo(AAA::SNonly));
	program.help_text_var("_Default_wt_del",str_YesOrNo(Wt_DEL));
	program.help_text_var("_Default_wt_vqs",str_YesOrNo(Wt_VQS));
	program.help_text_var("_Default_single",str_YesOrNo(single));
	program.help_text_var("_Default_out_lp",str_YesOrNo(out_lp));
	program.help_text_var("_Default_out_pv",str_YesOrNo(out_pv));
	program.help_text_var("_Default_out_or",str_YesOrNo(out_or));
	program.help_text_var("_Default_out_se",str_YesOrNo(out_se));
	program.help_text_var("_Default_1sided",str_YesOrNo(one_sided));
	program.help_text_var("_Default_biol",GS2GBA);
	program.help_text_var("_Default_neg_biol",str_YesOrNo(NegGBA));
	program.help_text_var("_Default_var_wt",str_YesOrNo(Var_Wt));
	program.help_text_var("_Default_phenocopy",ftos(PGuard));
	program.help_text_var("_Default_cap",ftos(cap_pg));
	program.help_text_var("_Default_top",str_YesOrNo(cap_ps));
	program.help_text_var("_Default_de_novo",itos(denovo));
	program.help_text_var("_Default_DNonly",str_YesOrNo(DNonly));
	program.help_text_var("_Default_DN_max",itos(DN_max));
	program.help_text_var("_Default_big_spl",itos(bigspl));
	program.help_text_var("_Default_mac_ge",itos(mac_ge));
	program.help_text_var("_Default_ped",ped_in);
	program.help_text_var("_Default_xct_prefix",str_of_container(ExtCT::ExAC_pfx,',',false));
	program.help_text_var("_Default_xct_pop",ExtCT::ExAC_pop);
	program.help_text_var("_Default_xct_padding",itos(ExtCT::ExAC_pad));
	program.help_text_var("_Default_xct_cov_pc",ftos(ExtCT::ExAC_pCv));
	program.help_text_var("_Default_xct_case_only",str_YesOrNo(ExtCT::ExAC_cso));
	program.help_text_var("_Default_xct_frq_pv",ftos(ExtCT::ExAC_fqc));
	program.help_text_var("_Default_xct_log",ExtCT::ExAC_log);
	program.help_text_var("_Default_sv",str_of_container(sv_inp,',',false));
	program.help_text_var("_Default_iwt_fn",iwt_fn);
	program.help_text_var("_Default_my_R",MyR_in);
	program.help_text_var("_Default_permute",itos(permte));
	program.help_text_var("_Default_one_file",str_YesOrNo(HasAll));
	program.help_text_var("_Default_perm_data",str_of_container(MxVlFs,',',false));
	program.help_text_var("_Default_immed_cor",str_YesOrNo(immedC));
	program.help_text_var("_Default_delay_cor",str_YesOrNo(delayC));
	program.help_text_var("_Default_const_cor",str_YesOrNo(constC));
	program.help_text_var("_Default_min_cohort",itos(FETmin));
	program.help_text_var("_Default_min_sample",itos(RNKmin));
	program.help_text_var("_Default_gtp",gtp_in);
	program.help_text_var("_Default_example_chr",ex_chr);
	program.help_text_var("_Default_recessive",str_YesOrNo(FETrec));
	program.help_text_var("_Default_trend",str_YesOrNo(h_test));
	program.help_text_var("_Default_prs",str_YesOrNo(usePRS));
	program.help_text_var("_Default_por",str_YesOrNo(usePOR));
	program.help_text_var("_Default_cal_paf",str_YesOrNo(PRSPAF));
	program.help_text_var("_Default_ok_no_match",str_YesOrNo(skipnm));
	program.help_text_var("_Default_do_nothing",str_YesOrNo(do_nothing));
	program.push_back_help(GTP::help_text());
	perch::check_arguments();
	
	if (do_nothing)
	{
		tfile_format format;
		format.set_delimiters("\t");
		format.set_option(SKIP_NOTES,false);
		for (Lines_in_File(in, var_in, &format)) program.outf << in[0] << endl;
		return 0;
	}
	
	// check parameters
	if ( gtp_in.empty())				exit_error("--gtp not set.");
	if ( h_grID.empty())				exit_error("Groups not set.");
	if (perch::penetr.empty())			exit_error("penetrance is not set.");
	if (perch::preval==0)				exit_error("prevalence is not set");
	if ( perch::penetr[0]==perch::penetr[1] && perch::penetr[0]<perch::penetr[2] ) FETrec = true;
	if (perch::VarCla)
	{
		Wt_DEL=Wt_VQS=false;
		//if (AnalysisMethod!=HLR)	exit_error("For variant classification, analysis method must be HLR.");
		if (permte)					exit_error("For variant classification, permutation is not allowed");
		if (Wt_DEL || Wt_VQS) 		exit_error("For variant classification, analysis should not be weighted");
		if (!GS2GBA.empty())		exit_error("For variant classification, analysis should not be weighted");
		perch::clear_fltdel(false);
		perch::clear_flt_af(false);
		perch::filflt.clear();
	}
	if (usePRS || usePOR)
	{
		Wt_DEL=Wt_VQS=false;
		if (permte)					exit_error("For polygenic score, permutation is not allowed");
		if (Wt_DEL || Wt_VQS) 		exit_error("For polygenic score, analysis should not be weighted");
		if (!GS2GBA.empty())		exit_error("For polygenic score, analysis should not be weighted");
		perch::clear_fltdel(false);
		perch::clear_flt_af(false);
		perch::filflt.clear();
	}
	if (ped_in.empty() && spl_in.empty()) exit_error("please use --spl or --ped to specify samples to be analyzed");
	if (denovo && ped_in.empty()) exit_error("--de-novo must be used together with --ped");
	if (DNonly && ped_in.empty()) exit_error("--dn-only must be used together with --ped");
	//if (DNonly && !denovo) exit_error("--dn-only must be used togethe with --de-novo");
	if (DNonly && !spl_in.empty()) { lns<<showw<<"--dn-only makes a valid test only when you don't have non-trio samples, i.e., do not use --spl to set a sample file."<<flush_logger; }
	if (showID && AnalysisMethod!=DET) exit_error("--show-id must be used together with --detail");
	if (showGT && AnalysisMethod!=DET) exit_error("--show-gt must be used together with --detail");
	if (showIF && AnalysisMethod!=DET) exit_error("--show-if must be used together with --detail");
	if (perch::DomNeg && !sv_inp.empty()) { lns<<showw<<"--dom-neg cannot be used together with --sv. The latter argument is ignored."<<flush_logger; sv_inp.clear(); }
	if (perch::MisCut!=1 && perch::MisCut>internal_MisCut && AnalysisMethod!=DET)
		exit_error("--filt-miss-rate cannot be too high"); // always remove variants that is missing for all samples, otherwise cannot calculate the gene even other variants are not missing
	if (!ExtCT::ExAC_pfx.empty())
	{
		if		(AnalysisMethod==GLR || AnalysisMethod==FET || AnalysisMethod==CL2 || AnalysisMethod==CL1 || AnalysisMethod==DET) ;
		else if (AnalysisMethod==HLR) { AnalysisMethod=GLR; } // previously wrWhat=perch::h_GLRb
		else { lns<<showw<<"unsupported analysis type by --xct-pfx, changed to --FET."<<flush_logger; AnalysisMethod=FET; }
		if (!sv_inp.empty()) { lns<<showw<<"--xct-pfx cannot be used together with --sv. The latter argument is ignored."<<flush_logger; sv_inp.clear(); }
		if (out_or) exit_error("ourput of odds ratio is not supported for --xct-pfx.");
	}
	else
	{
		if (out_or)
		{
			if 	 	(AnalysisMethod==REG) regress=&or_1st_logistic;
			else if (AnalysisMethod==FLR) ;
			else if (AnalysisMethod==FET) ;
			else exit_error("ourput of odds ratio is not supported for this analysis type. You need either --FET or --logistic.");
		}
		if (!iwt_fn.empty())
		{
			if (AnalysisMethod==REG && regress==&pv_1st_linear) ;
			else if (AnalysisMethod==MyR) ;
			else if (AnalysisMethod==FLR) ;
			else if (AnalysisMethod==FET) FETiwt=true;
			else if (AnalysisMethod==CL2) FETiwt=true;
			else if (AnalysisMethod==OUT) ;
			else { lns<<showw<<"--weight only works with --linear / --logistic / --write / --FET / --collapse. This option is ignored."<<flush_logger; iwt_fn.clear(); }
		}
		if (AnalysisMethod==SSU && !sv_inp.empty()) { lns<<showw<<"--SSUs / --SSUc cannot be used together with --sv. The latter argument is ignored."<<flush_logger; sv_inp.clear(); }
		if (AnalysisMethod==CL1 && !sv_inp.empty()) { lns<<showw<<"--count-gtp cannot be used together with --sv. The latter argument is ignored."<<flush_logger; sv_inp.clear(); }
		if (AnalysisMethod==FET || AnalysisMethod==CL2 || AnalysisMethod==RNK || AnalysisMethod==HLR || AnalysisMethod==POL || AnalysisMethod==OPN) perch::a_iPop=true;
	}

	if (!MyR_in.empty())
	{
		tfile_format MyR_format;
		MyR_format.set_option(SKIP_NOTES,true);
		MyR_format.set_option(SKIP_BLANKS,true);
		for (Lines_in_File(in,MyR_in,&MyR_format))
		{
			MyCode += in[0];
			MyCode += '\n';
		}
		if (MyCode.empty()) exit_error("nothing read from "+MyR_in);
	}
	
	// read permutation data
	if (!MxVlFs.empty())
	{
		for (Rows_in_File(in,MxVlFs,1))
		{
			double maxval = -std::numeric_limits<double>::max();
			for (int i=0;i<in.NumFields();++i)
			{
				double gotval;
				if (!read_val(in[i],gotval)) exit_error("cannot read "+in[i]+" in "+in.FileName());
				if (gotval>maxval) maxval=gotval;
			}
			wjp_res_mx.push_back(maxval);
		}
	}
	
	// read sample file
	set< string >					h_csID;			// SeqID of cases
	set< string >					h_ctID;			// SeqID of controls
	map< string, int >				SexMap;			// SeqID => gender (1 for male, 2 for female)
	map< string, double >			DepMap;			// SeqID => dependent variable
	map< string, string >			StrMap;			// SeqID => strata
	map< string, vector<double> >	CovMap;			// SeqID => covariate for analysis, [c] is the same as CovVar CovNID
	if (!spl_in.empty())
	{
		if (ToWait)
		{
			string first_line;
			safeGetline(cin,first_line);
			if (first_line!="## vSIM created genotype file")
				exit_error("Don't use --wait if genotype file is not in standard input or not created by vSIM.");
		}
		
		set< string >				h_ukID;	// SeqID of unknowns
		map< string, string> 		PopMap;	// SeqID => population
		bool						read_cov=true;
		if (!ExtCT::ExAC_pfx.empty() || (AnalysisMethod==SSU && do_SSUw==&pv_SSUw)) read_cov=false;
		perch::read_spl(spl_in,true,true,!read_cov,0.0,h_csID,h_ctID,h_ukID,SexMap,DepMap,PopMap,StrMap,CovMap,CovNID);
	}
	
	// read individual weight file
	map< string, double > iwt_db; // individual weight [SeqID]
	if (!iwt_fn.empty())
	{
		for (Rows_in_File(in,iwt_fn,2)) // required columns: SeqID weight
		{
			if (in.RowNumber()==0)
			{
				boost::to_lower(in[0]); if (in[0]=="sample" || in[0]=="indid" || in[0]=="iid" || in[0]=="seqid")	continue;
				boost::to_lower(in[1]); if (in[1]=="weight" || in[1]=="indwt" || in[1]=="iwt" || in[0]=="iweight")	continue;
			}
			double wt=std::numeric_limits<double>::signaling_NaN();
			if (read_val(in[1],wt))
			{
				if (wt<0) exit_error("individual weight should not be negative");
				if (wt>1) exit_error("individual weight should not be greater than 1");
				iwt_db[in[0]]=wt;
			}
			else exit_error("failed to read "+in[1]+" in "+iwt_fn+" as a weight value");
		}
	}
	
	// check parameters
	bool print_header=true;
	if (permte)
	{
		if (program.nt==1) lns<<showl<<"With --permute, it is better to use multi-threading. Run with 1 thread anyway."<<flush_logger;
		if (var_in.empty()) exit_error("With --permute, input must be a file instead of the standard input.");
		wrWhat=perch::h_Pcrt;
	}
	else
	{
		if (program.nt>1) lns<<showl<<"Multi-threading is only for permutations. Run with 1 thread only."<<flush_logger;
		if (AnalysisMethod==DET)
		{
			wrWhat=perch::h_Adet;
			if (showGT||showIF)
			{
				print_header=false;
				program.outf<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t"<<perch::h_symb<<"\t"<<perch::h_func<<"\t"<<perch::h_fdet<<"\t"<<wrWhat<<endl;
			}
		}
		if (AnalysisMethod==CL1 || AnalysisMethod==CL2)
		{
			print_header=false;
			wrWhat=perch::h_Acol;
			if (ExtCT::ExAC_pfx.empty())	program.outf<<perch::h_symb<<"\t#CHROM\tNo.Var\t"<<perch::h_mLgP<<"\t"<<perch::h_o_r_<<"\t95%_Confidence_Interval\tcs_neg\tcs_het\tcs_hom\tct_neg\tct_het\tct_hom"<<endl;
			else							program.outf<<perch::h_symb<<"\t#CHROM\tNo.Var\t"<<perch::h_mLgP<<"\tcs_ref\tcs_alt\tct_ref\tct_alt\text\tstudy\tboth"<<endl;
		}
		if (AnalysisMethod==OUT)
		{
			print_header=false;
		}
		if (out_lp) { if (AnalysisMethod==FET||AnalysisMethod==MyR||AnalysisMethod==FLR||AnalysisMethod==REG||AnalysisMethod==SSU||AnalysisMethod==GLR||AnalysisMethod==RNK) wrWhat=perch::h_mLgP; else exit_error("--out-mlp only works with --linear or --SSUw or --FET"); }
		if (out_pv) { if (AnalysisMethod==FET||AnalysisMethod==MyR||AnalysisMethod==FLR||AnalysisMethod==REG||AnalysisMethod==SSU||AnalysisMethod==GLR||AnalysisMethod==RNK) wrWhat=perch::h_pVal; else exit_error("--out-pv only works with --linear or --SSUw or --FET"); }
		if (out_or) { if (AnalysisMethod==FET||AnalysisMethod==FLR||(AnalysisMethod==REG&&regress==&or_1st_logistic)) wrWhat=perch::h_o_r_; else exit_error("--out-or only works with --logistic or --FET"); }
		if (AnalysisMethod==SSU && !CovMap.empty()) lns<<showw<<"For SSU with covariates, please use permutations to obtain p-values."<<flush_logger;
	}
	
	// GBA
	map<string,double> S2G;
	if (!GS2GBA.empty())
	{
		int num_columns=0; for (Rows_in_File(in,GS2GBA,1)) { num_columns=in.NumFields(); break; }
		if (num_columns>2)
		{
			vector<string> GeneSymbols;
			for (Rows_in_File(in,GS2GBA,1))
			{
				if (in.RowNumber()==0) {	GeneSymbols=in.contents(); continue; }
				if (in.RowNumber()==1)
				{
					for (int i=1;i<in.NumFields();++i)
					{
						double v = 0;
						try { v=boost::lexical_cast<double>(in[i]); }
						catch (...) { exit_error("Error reading "+GS2GBA+": failed to read "+in[i]+" as a number."); }
						if (!NegGBA && v<0) v=0;
						S2G[GeneSymbols[i]] = v;
					}
					break;
				}
			}
		}
		else if (num_columns==2)
		{
			for (Rows_in_File(in,GS2GBA,2))
			{
				if (exist_element(aLBF_h,to_lower_copy(in[0]))) continue;
				double v = 0;
				try { v=boost::lexical_cast<double>(in[1]); }
				catch (...) { exit_error("Error reading "+GS2GBA+": failed to read "+in[1]+" as a number."); }
				if (!NegGBA && v<0) v=0;
				S2G[in[0]] = v;
			}
		}
		else exit_error(GS2GBA+" has the wrong number of columns.");
		Wt_GBA = true;
	}
	
	if (!sv_inp.empty()) AAA_SV::sv_read_wiSymb(sv_inp);
	if (showIF) AAA_IF::read_panels();
	ExtCT::read_files();

	// read VCF header
	vector<int>		SexCSs;				// same order as FldCSs
	vector<int>		SexCTs;				// same order as FldCTs
	field_numbers	FldCSs(false,true);	// field numb for cases
	field_numbers	FldCTs(false,true);	// field numb for controls
	field_numbers	FldChr(false,true);	// field numb for #CHROM
	field_numbers	FldPos(false,true);	// field numb for POS
	field_numbers	FldRef(false,true);	// field numb for REF
	field_numbers	FldAlt(false,true);	// field numb for ALT
	field_numbers	FldFmt(false,true);	// field numb for FORMAT
	tfile_format	vcf_format;
	vcf_format.set_delimiters("\t");
	vcf_format.set_option(SKIP_NOTES,false);
	bool vcf_header_not_read = true;
	size_t	max_pg=0;	// max number of variants per gene, 0=no_limit
	int     max_ps=0;	// max number of variants per sample, 0=no_limit
	vector<int> vc_cs, vc_ct;
	vector<perch::trio> tr_cs, tr_ct;// col_pa[sampleID] is 0-based column of dad, 0=not_sequenced
	map<string,perch::trio> trios;// trios[SeqID]=trio (affected child unaffected parents)
	bool trios_read=false;
	string chr_spe_gtp = boost::algorithm::replace_all_copy(gtp_in,"@",ex_chr);
	for (Rows_in_File(in, chr_spe_gtp, &vcf_format))
	{
		// read header
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
		{
			in.clear_nf();
			continue;
		}
		if (exist_any(perch::h_col1, in.contents()))
		{
			if (!vcf_header_not_read) { continue; }
			vcf_header_not_read=false;
			vcf_format.clear_field_nums();
			FldCSs.clear();
			FldCTs.clear();
			FldChr.clear();
			FldPos.clear();
			FldRef.clear();
			FldAlt.clear();
			FldFmt.clear();

			// read trios
			if (!ped_in.empty() && !trios_read)
			{
				perch::get_trios(ped_in,in,trios);
				if ((denovo||DNonly) && trios.empty()) exit_error("No trios found.");
				// treat them as independent cases & controls. Two sibs with the same phenotype does not fit with the hypothesis of de novo mutation.
				for (auto &i:trios)
				{
					perch::trio& t = i.second;
					if		(t.aff==2) { SexMap[t.id[0]]=t.sex[0]; DepMap[t.id[0]]=1; h_csID.insert(t.id[0]); }
					else if (t.aff==1) { SexMap[t.id[0]]=t.sex[0]; DepMap[t.id[0]]=0; h_ctID.insert(t.id[0]); }
				}
				trios_read=true;
			}

			// read case-control and other fields
			for (int i=0;i<in.NumFields();++i)
			{
				if (exist_element(h_csID,in[i])) {	FldCSs.push_back(i+1); SexCSs.push_back(SexMap[in[i]]); }
				if (exist_element(h_ctID,in[i])) {	FldCTs.push_back(i+1); SexCTs.push_back(SexMap[in[i]]); }
				if (in[i]=="Chr"   && FldChr.no_input())	FldChr.push_back(i+1); if (in[i]=="#CHROM")		{	FldChr.clear(); FldChr.push_back(i+1); }
				if (in[i]=="Start" && FldPos.no_input())	FldPos.push_back(i+1); if (in[i]=="POS")		{	FldPos.clear(); FldPos.push_back(i+1); }
				if (in[i]=="Ref"   && FldRef.no_input())	FldRef.push_back(i+1); if (in[i]=="REF")		{	FldRef.clear(); FldRef.push_back(i+1); }
				if (in[i]=="Alt"   && FldAlt.no_input())	FldAlt.push_back(i+1); if (in[i]=="ALT")		{	FldAlt.clear(); FldAlt.push_back(i+1); }
				if (in[i]=="FORMAT") FldFmt.push_back(i+1);
			}

			// check errors or warnings
			if (FldChr.no_input()) exit_error("The #CHROM/Chr column is missing.");
			if (FldPos.no_input()) exit_error("The POS/Start column is missing.");
			if (FldRef.no_input()) exit_error("The REF/Ref column is missing.");
			if (FldAlt.no_input()) exit_error("The ALT/Alt column is missing.");
			if (FldCSs.no_input() &&  (AnalysisMethod==SSU||AnalysisMethod==FLR||AnalysisMethod==REG||AnalysisMethod==HLR||AnalysisMethod==FET||AnalysisMethod==OPN||AnalysisMethod==POL||AnalysisMethod==RNK)) exit_error("Find no cases.");
			if (FldCTs.no_input() &&  (AnalysisMethod==SSU||AnalysisMethod==FLR||AnalysisMethod==HLR||AnalysisMethod==FET||AnalysisMethod==OPN||AnalysisMethod==POL||AnalysisMethod==RNK) && ExtCT::ExAC_pfx.empty()) exit_error("Find no controls.");
			//if (perch::is_qtl() && !(AnalysisMethod==DET||AnalysisMethod==SSU||(AnalysisMethod==REG&&regress==&pv_1st_linear))) exit_error("QTL not supported");
			if (FldCSs.size()<h_csID.size()) lns<<showl<<h_csID.size()-FldCSs.size()<<" cases in the Sample File cannot be found in the Genotype File."<<flush_logger;
			if (FldCTs.size()<h_ctID.size()) lns<<showl<<h_ctID.size()-FldCTs.size()<<" controls in the Sample File cannot be found in the Genotype File."<<flush_logger;
			lns<<showl<<"There are "<<FldCSs.size()<<" cases and "<<FldCTs.size()<<" controls remaining for analysis."<<flush_logger;
			if (AnalysisMethod==FET)
			{
				if (!ExtCT::ExAC_pfx.empty())
				{
					if ((int)FldCSs.size()<FETmin)
						lns<<showe<<"The sample size is small. Please reduce the --min-cohort argument to run analysis, or increase the sample size."<<fatal;
				}
				else
				{
					if ((int)FldCSs.size()<FETmin||(int)FldCTs.size()<FETmin)
						lns<<showe<<"The sample size is small. Please reduce the --min-cohort argument to run analysis, or increase the sample size."<<fatal;
				}
			}
			if (AnalysisMethod==HLR || AnalysisMethod==POL || AnalysisMethod==OPN || AnalysisMethod==CL2 || AnalysisMethod==FET || (AnalysisMethod==SSU && do_SSUw==&pv_SSUw))
			{
				set<string> SubGrp;
				for (size_t i=0;i<FldCSs.size();++i) SubGrp.insert(StrMap[in[FldCSs[i]]]);
				if (((double)FldCSs.size())/SubGrp.size()<FETmin)
					exit_error("Too many strata; "+itos(FldCSs.size())+" cases are divided into "+itos(SubGrp.size())+" groups. Make sure covariates are NOT continuous variables.");
				
				SubGrp.clear();
				for (size_t i=0;i<FldCTs.size();++i) SubGrp.insert(StrMap[in[FldCTs[i]]]);
				if (((double)FldCTs.size())/SubGrp.size()<FETmin)
					exit_error("Too many strata; "+itos(FldCTs.size())+" controls are divided into "+itos(SubGrp.size())+" groups. Make sure covariates are NOT continuous variables.");
			}

			// prepare data
			if (!SexCSs.empty()) { ExtCT::Case_fem=0; for (size_t i=0;i<SexCSs.size();++i) if (SexCSs[i]==2) ++ExtCT::Case_fem; ExtCT::Case_fem/=SexCSs.size(); }
			if (!SexCTs.empty()) { ExtCT::Ctrl_fem=0; for (size_t i=0;i<SexCTs.size();++i) if (SexCTs[i]==2) ++ExtCT::Ctrl_fem; ExtCT::Ctrl_fem/=SexCTs.size(); }
			for (size_t i=0;i<FldCSs.size();++i) {
				CovVar.push_back(CovMap[in[FldCSs[i]]]);
				DepVar.push_back(DepMap[in[FldCSs[i]]]);
				strata.push_back(StrMap[in[FldCSs[i]]]);
				seq2co[in[FldCSs[i]]]=record.size();
				record.push_back(in[FldCSs[i]]);
				if (exist_element(iwt_db,in[FldCSs[i]])) { sqrtWt.push_back(sqrt(iwt_db[in[FldCSs[i]]])); origWt.push_back(iwt_db[in[FldCSs[i]]]); } else { sqrtWt.push_back(1); origWt.push_back(1); }
			}
			for (size_t i=0;i<FldCTs.size();++i) {
				CovVar.push_back(CovMap[in[FldCTs[i]]]);
				DepVar.push_back(DepMap[in[FldCTs[i]]]);
				strata.push_back(StrMap[in[FldCTs[i]]]);
				seq2co[in[FldCTs[i]]]=record.size();
				record.push_back(in[FldCTs[i]]);
				if (exist_element(iwt_db,in[FldCTs[i]])) { sqrtWt.push_back(sqrt(iwt_db[in[FldCTs[i]]])); origWt.push_back(iwt_db[in[FldCTs[i]]]); } else { sqrtWt.push_back(1); origWt.push_back(1); }
			}
			for (size_t i=0;i<FldCSs.size();++i) { perch::trio t=trios[in[FldCSs[i]]]; t.id[0]=in[FldCSs[i]]; t.col[0]=FldCSs[i]; t.sex[0]=SexMap[in[FldCSs[i]]]; tr_cs.push_back(t); }
			for (size_t i=0;i<FldCTs.size();++i) { perch::trio t=trios[in[FldCTs[i]]]; t.id[0]=in[FldCTs[i]]; t.col[0]=FldCTs[i]; t.sex[0]=SexMap[in[FldCTs[i]]]; tr_ct.push_back(t); }
			NumCSs = FldCSs.size();
			NumCTs = FldCTs.size();
			if ((AnalysisMethod==HLR || AnalysisMethod==POL || AnalysisMethod==OPN) && PGuard && FldCSs.size()>1 && !FldCTs.empty())
			{
				double phenocopy = std::ceil(FldCSs.size()*PGuard) / (FldCSs.size()+FldCTs.size());
				if (perch::penetr[0]<phenocopy)
				{
					lns<<showl<<"The penetrance ("<<str_of_container(perch::penetr,',')<<") is too low, not allowing a sufficient phenocopy rate in the sequenced samples.";
					perch::penetr[1] = truncated(perch::penetr[1]+(phenocopy-perch::penetr[0]));
					perch::penetr[2] = truncated(perch::penetr[2]+(phenocopy-perch::penetr[0]));
					perch::penetr[0] = truncated(phenocopy);
					lns<<" Now changed to "<<str_of_container(perch::penetr,',')<<flush_logger;
				}
			}
			if ( FldCTs.size()/(1-perch::preval) >=bigspl ) // smallest_maf = min(P/cs,(1-P)/ct), because P<0.5, it's (1-P)/ct that matters. So sample size = 1/smallest_maf = ct/(1-P).
			{
				XAF_fp = false;
				// cerr << "# Sample size is big enough (>="<<bigspl<<"), population frequency will be calculated from data.\n";
			}
			else
			{
				XAF_fp = true;
				// cerr << "# Sample size not big enough (<"<<bigspl<<"), population frequency will be calculated from MaxAF.\n";
			}
			for (size_t i=0;i<(FldCSs.size()+FldCTs.size());++i) OriPtr.push_back(i);
			DefaultJob.setup(FldCSs.size()+FldCTs.size(),"");	// guard for the situation that no variant has passed to the CalHLR in the whole file, write_results() will be called anyway after "for (Rows_in_File)" has ended.
			PmtPtr.push_back(OriPtr);
			for (int i=0;i<permte;++i)
			{
				PmtPtr.push_back(OriPtr);
				std::random_shuffle(PmtPtr.back().begin(),PmtPtr.back().end(),MyRandom);
			}
			if (cap_pg) max_pg = ceil(FldCSs.size()*cap_pg);
			if (cap_ps) max_ps = (perch::penetr[1]==perch::penetr[2] ? 1 : 2);

			// prepare adaptive MAF cutoff
			if (Adaptv && perch::filSAF)
			{
				// increase the frequency cutoff to accommodate the fact that the risk allele has a higher frequency in cases than in the population,
				vector<double> newpen = perch::penetr;
				newpen[1] = truncated(newpen[1]==newpen[0] ? newpen[0] : 10 * newpen[1]);
				newpen[2] = truncated(newpen[2]==newpen[0] ? newpen[0] : 10 * newpen[2]);
				double CsFcut = genepi::up_af(perch::filSAF,newpen);
				double CtFcut = perch::filSAF;

				// no up-lift MAF to at least 1 observation among samples.
				// perch::filSAF = (FldCSs.size()*2*CsFcut + FldCTs.size()*2*CtFcut) / (FldCSs.size()*2 + FldCTs.size()*2);

				//    up-lift MAF to at least 1 observation among samples.
				perch::filSAF = std::max( (FldCSs.size()*2*CsFcut + FldCTs.size()*2*CtFcut), 1.0 ) / (FldCSs.size()*2 + FldCTs.size()*2);
			}
			
			continue;
		}
		if (vcf_header_not_read) exit_error("Header lines missing in the input Genotype File.");
		break;
	}
	if (vcf_header_not_read) exit_error("Header lines missing in the input Genotype File.");

	field_numbers	FldFlt(false,true);	// field numb for FILTER
	field_numbers	FldInf(false,true);	// field numb for INFO
	field_numbers	FldGrp(false,true);	// field numb for group ID
	field_numbers	FldSym(false,true);	// field numb for GeneSymbol
	field_numbers	FldFun(false,true);	// field numb for FuncConseq
	field_numbers	FldDet(false,true);	// field numb for FuncDetail
	field_numbers	FldXAF(false,true);	// field numb for MaxAF
	field_numbers	FldRes(false,true);	// field numb for result
	tfile_format	var_format;
	var_format.set_delimiters("\t");
	var_format.set_option(SKIP_NOTES,false);
	var_format.set_storage_to(program.main_data());
	program.main_data().clear();
	program.main_data().keep_all_input();
	int	ColDel = -2; // -2 = not annotated; -1 = in INFO; 0+ = column.
	int ColSEG = -2; // column number for vSEG analysis results as log10 Bayes factors
	int warning1 = elog.get_token("duplicated variants within a group were skipped.");
	int warning2 = elog.get_token("duplicated groups.");
	int warning3 = elog.get_token("variants have no match in VCF");
	bool var_header_not_read = true;
	string		prev_group;	// the group right before this, to see whether it's a new group
	set<string> old_groups; // all groups before this, to check whether rows are sorted by groups
	set<string> prev_index;	// chr-pos-id-ref-alt
	string ExAC_VAG;

read_var_group:
	old_groups.clear();
	prev_group.clear();
	prev_index.clear();
	for (Rows_in_File(in, var_in, &var_format))
	{
		// read header
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
		{
			in.clear_nf();
			perch::read_meta(in[0]);
			if (str_startsw(in[0],"##INFO=<ID="+perch::h_SEGb+",")) ColSEG=-1;
			if (str_has(in[0],"##INFO=<ID="+wrWhat+",")) { program.main_data().clear(); continue; }
			if (str_has(in[0],"##INFO=<ID="+perch::h_Axtr+",")) { program.main_data().clear(); continue; }
			if (str_has(in[0],"##INFO=<ID=vAAA_SWAP,")) { program.main_data().clear(); continue; }
			if (!var_header_not_read) { program.main_data().clear(); continue; }
			if ( h_dels.empty() && str_startsw(in[0],"##INFO=<ID=BayesDel")) ColDel=-1;
			if (!h_dels.empty() && str_startsw(in[0],"##INFO=<ID="+h_dels+",")) ColDel=-1;
			if (str_has(in[0],"##VictorCommandLine=<ID=vAnnGene,")) { ExAC_VAG=in[0].substr(52); ExAC_VAG.pop_back(); }
			if (!ExtCT::ExAC_VAG.empty() && !ExAC_VAG.empty() && ExAC_VAG!=ExtCT::ExAC_VAG) exit_error("The study and ExAC VCF were not annotated with the same parameters.");
			if (print_header) print_container(in.contents(),program.outf,' ',true);
			program.main_data().clear();
			continue;
		}

		if (exist_any(perch::h_col1, in.contents()))
		{
			if (AddInf) program.outf<<"##INFO=<ID="<<wrWhat<<",Number=1,Type=Float,Description=\"vSEG co-segregation analysis result.\">"<<endl;
			program.outf<<"##INFO=<ID="+perch::h_Axtr+",Number=1,Type=String,Description=\"vAAA addition results.\">"<<endl;
			program.outf<<"##INFO=<ID=vAAA_SWAP,Number=0,Type=Flag,Description=\"vAAA decided to swap REF and ALT due to allele frequency.\">"<<endl;
			if (!var_header_not_read) { program.main_data().clear(); continue; }
			lns<<showl<<"BayesDel filter cutoff value is "<<perch::FltDel<<flush_logger;
			var_header_not_read=false;
			var_format.clear_field_nums();
			FldGrp.clear();
			FldFlt.clear();
			FldInf.clear();
			FldSym.clear();
			FldFun.clear();
			FldDet.clear();
			FldXAF.clear();
			FldRes.clear();
			
			// read case-control and other fields
			for (int i=0;i<in.NumFields();++i)
			{
				if (in[i]==perch::h_SEGb)			ColSEG=i;
				if (in[i]=="FILTER")				FldFlt.push_back(i+1);
				if (in[i]=="INFO")					FldInf.push_back(i+1);
				if (in[i]==perch::h_MxAF)			FldXAF.push_back(i+1);
				if (in[i]==perch::h_symb)			FldSym.push_back(i+1);
				if (in[i]==perch::h_func)			FldFun.push_back(i+1);
				if (in[i]==perch::h_fdet)			FldDet.push_back(i+1);
				if (exist_element(h_grID,in[i]))	FldGrp.push_back(i+1);
				if ( h_dels.empty() && str_startsw(in[i],"BayesDel"))	{ if (ColDel==-2) ColDel=i; else exit_error("multiple columns for BayesDel"); }
				if (!h_dels.empty() && in[i]==h_dels)					{ if (ColDel==-2) ColDel=i; else exit_error("multiple columns for "+h_dels);	}
				if (in[i]=="Chr"   || in[i]=="#CHROM")	{	if (FldChr[0]!=i) exit_error("The #CHROM column in Group File not match with the Genotype File."); 	}
				if (in[i]=="Start" || in[i]=="POS")		{	if (FldPos[0]!=i) exit_error("The POS column in Group File not match with the Genotype File."); 	}
				if (in[i]=="Ref"   || in[i]=="REF")		{	if (FldRef[0]!=i) exit_error("The REF column in Group File not match with the Genotype File.");		}
				if (in[i]=="Alt"   || in[i]=="ALT")		{	if (FldAlt[0]!=i) exit_error("The ALT column in Group File not match with the Genotype File.");		}
			}
			
			// check errors or warnings
			if ( AddInf && FldInf.no_input()) exit_error("The INFO column is missing.");
			if ( Wt_VQS && FldInf.no_input()) exit_error("The INFO  column is missing.");
			if ( Wt_GBA && FldSym.no_input()) exit_error("The Gene Symbol column is missing.");
			if ( Wt_DEL && ColDel==-2)   exit_error("The BayesDel annotation is missing.");
			if ( perch::LFonly && FldFun.no_input()) exit_error("The Function column is missing.");
			if ((perch::VarCla || single) && FldSym.no_input()) exit_error("The Gene Symbol column is missing.");
			if ( perch::VarCla || single) FldGrp = FldChr + FldPos + FldRef + FldAlt + FldSym;
			if (FldGrp.no_input()) exit_error("The group ID column is missing.");
			
			// prepare output
			if (AddInf)	FldRes=FldInf;
			else	  { FldRes.push_back(in.NumFields()+1); in.contents().push_back(wrWhat); var_format.set_field_nums(FldRes,"",tfile_format::Expand); }
			if (print_header) print_container(in.contents(),program.outf,DLMTR,true);
			program.main_data().clear();
			continue;
		}
		if (var_header_not_read) exit_error("Header lines missing.");

		// check file
		string VarIndex = in[FldChr[0]]+'_'+in[FldPos[0]]+'_'+in[FldRef[0]]+'_'+in[FldAlt[0]];
		if (!FldSym.no_input() && in[FldSym[0]].find_first_of(",;|")!=std::string::npos) exit_error("Multiple gene symbols for a variant is not allowed.");
		if (in[FldAlt[0]].find(',')!=std::string::npos) exit_error("The input file has not been split by alternative alleles.");

		string GeneSymbol;
		if (!FldSym.no_input()) GeneSymbol=in[FldSym[0]];
		if (str_has(GeneSymbol,"_CHR")) GeneSymbol=substr_before_find(GeneSymbol,"_CHR");
		if (str_startsw(GeneSymbol,"ENSG00") && str_has(GeneSymbol,"(")) { GeneSymbol=substr_after_find(GeneSymbol,"("); GeneSymbol.pop_back(); }

		// skip chrUnknown
		int chr_num = genepi::read_chr_num(in[FldChr[0]]);

		// group & index
		string this_group;
		FldGrp.contents_to_a_string(in.contents(),this_group,','); // previously DLMTR
		if (this_group!=prev_group)
		{
			if (exist_element(old_groups,this_group)) elog.add(warning2, this_group); // exit_error("The Genotype File is not sorted by groups.");
			else old_groups.insert(this_group);
			if (!prev_group.empty()) write_results(AddInf,FldInf,FldRes,wrWhat,false);
			prev_group=this_group;
			prev_index.clear();
			DefaultJob.setup(FldCSs.size()+FldCTs.size(),this_group);
			vc_cs.assign(FldCSs.size(),0);
			vc_ct.assign(FldCTs.size(),0);
		}
		else
		{
			if (exist_element(prev_index,VarIndex)) { elog.add(warning1); lns<<showw<<"Duplicated variant "<<VarIndex<<flush_logger; DefaultJob.qcl.push_back("SKIPPED=DupVar"); continue; }
		}
		prev_index.insert(VarIndex);
		
		// skip if not an SNV || is an SV
		boost::to_upper(in[FldRef[0]]);
		boost::to_upper(in[FldAlt[0]]);
		string& ref = in[FldRef[0]];
		string& alt = in[FldAlt[0]];
		bool is_snv = ( (ref=="A" || ref=="T" || ref=="C" || ref=="G") && (alt=="A" || alt=="T" || alt=="C" || alt=="G") );
		bool is_lof = ( FldFun.no_input() ? false : str_has(in[FldFun[0]],"LoF")||str_has(in[FldFun[0]],"NMD") );
		bool is_vks = ( FldFun.no_input() ? false : str_has(in[FldFun[0]],"knClinSig=1") );
		bool is_DNg = ( FldFun.no_input() ? false : perch::is_DomNeg(in[FldFun[0]]) );
		bool is_cds = ( FldFun.no_input() ? false : perch::is_coding(in[FldFun[0]]) );
		if (chr_num<1)												{ DefaultJob.qcl.push_back("SKIPPED=Chromosome"); continue; }
		if (perch::LFonly && is_lof && perch::is_LoFtol(GeneSymbol)){ DefaultJob.qcl.push_back("SKIPPED=LoF-tolerant"); continue; }
		if (perch::LFonly && !is_lof)								{ DefaultJob.qcl.push_back("SKIPPED=NotLoF"); continue; }
		if (perch::CDonly && !is_cds)								{ DefaultJob.qcl.push_back("SKIPPED=NotCDS"); continue; }
		if (perch::DomNeg && !is_DNg)								{ DefaultJob.qcl.push_back("SKIPPED=NotDomNeg"); continue; }
		if (AAA::SNonly && !is_snv)									{ DefaultJob.qcl.push_back("SKIPPED=NotSNV"); continue; }
		
		// skip if not in a targeted region (cs) or well covered region (ct), or in an excluded region (eg, genomicSuperDups )
		int bp=-1; try { bp = boost::lexical_cast<int>(in[FldPos[0]]); } catch (...) { exit_error("Failed to read "+in[FldPos[0]]+" as a position in basepairs."); }
		if (!perch::within_covered_region(chr_num,bp)) { DefaultJob.qcl.push_back("SKIPPED=ChrRegion"); continue; }
		bool is_MHC = genepi::is_MHC(chr_num,bp);

		// skip by ExAC_qc, used for GLR
		if (!ExtCT::ExAC_pfx.empty())
		{
			// skip by exclusion
			if (exist_element(ExtCT::ExAC_QC,VarIndex)) { DefaultJob.qcl.push_back("SKIPPED=GLR_Ref_QC"); continue; }
			
			// skip by coverage
			bool covered=true;
			int beg=bp, end=bp;
			if (ref.size()>1) end=end+ref.size()-1;
			if (alt.size()>1) end=end+1;
			for (int j=beg;j<=end;++j)
			{
				if (!exist_element(ExtCT::StudySpl,chr_num) || !exist_element(ExtCT::StudySpl[chr_num],j)) { covered=false; break; }
				if (!exist_element(ExtCT::ExAC_spl,chr_num) || !exist_element(ExtCT::ExAC_spl[chr_num],j)) { covered=false; break; }
			}
			if (!covered) { DefaultJob.qcl.push_back("SKIPPED=GLR_Ref_QC"); continue; }
		}

		// get INFO
		vector<string> INFO;
		if (!FldInf.no_input())
		{
			if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";"));
			bool info_modified=false;
			for (vector<string>::iterator it = INFO.begin(); it != INFO.end(); it++)
				if (str_startsw(*it,wrWhat+"=")) { INFO.erase(it); info_modified=true; break; }
			for (vector<string>::iterator it = INFO.begin(); it != INFO.end(); it++)
				if (*it=="vAAA_SWAP") { INFO.erase(it); info_modified=true; break; }
			if (info_modified)
			{
				in[FldInf[0]] = str_of_container(INFO,';');
				if (in[FldInf[0]].empty()) in[FldInf[0]]=".";
			}
		}
		
		// skip by MaxAF
		double MaxAF = std::numeric_limits<double>::signaling_NaN();
		if (!FldXAF.no_input())	read_val(in[FldXAF[0]],MaxAF);
		else					MaxAF=get_value(INFO,perch::h_MxAF);
		if (std::isnan(MaxAF)) MaxAF=0;
		if (perch::filXAF)
		{
			double f = (MaxAF>0.5 ? 1-MaxAF : MaxAF);
			if (!perch::rf_XAF && f> perch::filXAF) { DefaultJob.qcl.push_back("SKIPPED=MaxAF"); continue; }
			if ( perch::rf_XAF && f<=perch::filXAF) { DefaultJob.qcl.push_back("SKIPPED=MaxAF_rev"); continue; }
		}
		
		// skip by BayesDel
		double BayesDel = std::numeric_limits<double>::signaling_NaN();
		if		(ColDel==-1)	BayesDel = get_value_sw(INFO,"BayesDel");
		else if (ColDel>=0)	read_val(in[ColDel],BayesDel);
		if (perch::filter_AnnAF(BayesDel,GeneSymbol,is_lof,is_vks)) { DefaultJob.qcl.push_back("SKIPPED=BayesDel"); continue; }
		
		// Above filter also applies to ExAC, so no need to do ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear();
		
		// skip if has Mendelian error.
		if (ColSEG!=-2 && FltLinkageErr)
		{
			double v=std::numeric_limits<double>::signaling_NaN();
			perch::read_variable(in,ColSEG,INFO,perch::h_SEGb,v);
			if (std::isnan(v)) { DefaultJob.qcl.push_back("SKIPPED=LinkageErr"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
		}

		// skip by VICTOR_QC
		if (!FldInf.no_input())
		{
			string vQC = get_string(INFO,"VICTOR_QC");
			if (!vQC.empty() && vQC!="PASS") { DefaultJob.qcl.push_back("SKIPPED=vQC="+vQC); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
		}
		
		// skip by VQSLOD or hard filter. Use double so that it is nan if not exist. nan automatically return false in comparison, so no need to check isnan().
		double VQSLOD = get_value(INFO,"VQSLOD");
		/*if (std::isnan(VQSLOD))
		{
			if (perch::VQSnan)	{ DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
			if (perch::HFnoVQ || perch::hardft)
			{
				if (perch::FiltQD) { if (get_value(INFO,"QD")				<perch::FiltQD) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
				if (perch::FiltMQ) { if (get_value(INFO,"MQ")				<perch::FiltMQ) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
				if (perch::FiltFS) { if (get_value(INFO,"FS")				>perch::FiltFS) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
				if (perch::FiltHS) { if (get_value(INFO,"HaplotypeScore")	>perch::FiltHS) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
				if (perch::FiltMR) { if (get_value(INFO,"MQRankSum")		<perch::FiltMR) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
				if (perch::FiltRP) { if (get_value(INFO,"ReadPosRankSum")	<perch::FiltRP) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
			}
		}
		else
		{
			if (is_snv) { if (VQSLOD<perch::VQSsnv) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
			else		{ if (VQSLOD<perch::VQSidl) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
			if (perch::hardft)
			{
				if (perch::FiltQD) { if (get_value(INFO,"QD")				<perch::FiltQD) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
				if (perch::FiltMQ) { if (get_value(INFO,"MQ")				<perch::FiltMQ) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
				if (perch::FiltFS) { if (get_value(INFO,"FS")				>perch::FiltFS) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
				if (perch::FiltHS) { if (get_value(INFO,"HaplotypeScore")	>perch::FiltHS) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
				if (perch::FiltMR) { if (get_value(INFO,"MQRankSum")		<perch::FiltMR) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
				if (perch::FiltRP) { if (get_value(INFO,"ReadPosRankSum")	<perch::FiltRP) { DefaultJob.qcl.push_back("SKIPPED=VQSLOD|HardFilter"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } }
			}
		} //*/
		if (std::isnan(VQSLOD)) VQSLOD=0; // it's important for GLR, although it should not happen

		/*/ skip by FILTER
		if (!FldFlt.no_input())
			if (!perch::filflt.empty() && !exist_element(perch::filflt,in[FldFlt[0]])) { DefaultJob.qcl.push_back("SKIPPED=FILTER"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } //*/
		
		// find matched row in Genotype File
		bool matched=false;
		string vcf_tbx;	// return by tabix
		string chr = in[FldChr[0]];
		string chr_spe_gtp = boost::algorithm::replace_all_copy(gtp_in,"@",chr);
		try { vcf_tbx = exec("tabix "+chr_spe_gtp+" "+chr+":"+itos(bp)+"-"+itos(bp),false); }
		catch (const std::exception& error) { exit_error("tabix "+chr_spe_gtp+" failed"); }
		vector<string> gt;
		if (!vcf_tbx.empty())
		{
			vcf_tbx.pop_back(); // vcf_tbx ends with \n
			vector<string> vcf_row;
			boost::split(vcf_row,vcf_tbx,boost::is_any_of("\n"));
			for (auto &row:vcf_row)
			{
				boost::split(gt,row,boost::is_any_of("\t"));
				if (gt[3]==ref && gt[4]==alt && gt[1]==in[FldPos[0]]) { matched=true; break; }
			}
		}
		if (!matched)
		{
			if (skipnm) { DefaultJob.qcl.push_back("SKIPPED=NoMatch"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); elog.add(warning3); continue; }
			exit_error("cannot find a match when tabix "+chr_spe_gtp+" "+chr+":"+itos(bp)+"-"+itos(bp)+" and look for "+ref+" > "+alt);
		}
		
		// read FORMAT
		genepi::gtp_par gpar;
		if (!FldFmt.no_input()) gpar.read(in[FldFmt[0]]);

		// read genotypes and skip if no variant observed.
		double	xiA=0, xiU=0;	// count alternative alleles
		double  niA=0, niU=0;	// count chromosomes w/o missing
		int		iiA=0, iiU=0;	// count samples valid
		int		miA=0, miU=0;	// count samples missing
		int		aiA=0, aiU=0;	// count samples with 1+ alt
		string	g_cs, g_ct;
		int		ndn=0; // number of trios has de novo
		double PAF = std::numeric_limits<double>::signaling_NaN(); // Population allele frequency
		if (trios.empty()) // old version, will never happen
		{
			size_t j=FldCSs.size();
			GTP::rewind(); for (size_t i=0;i<FldCSs.size();++i) { char g=GTP::read(gt[FldCSs[i]],chr_num,bp,SexCSs[i],gpar); g_cs+=g; GTP::count(g,iiA,miA,niA,xiA,origWt[i]); }
			GTP::rewind(); for (size_t i=0;i<FldCTs.size();++i) { char g=GTP::read(gt[FldCTs[i]],chr_num,bp,SexCTs[i],gpar); g_ct+=g; GTP::count(g,iiU,miU,niU,xiU,origWt[i+j]); }
		}
		else
		{
			size_t j=FldCSs.size();
			GTP::rewind(); for (size_t i=0;i<FldCSs.size();++i) { char g=_read_trio(tr_cs[i],gt,chr_num,bp,DNonly,ndn,gpar); g_cs+=g; GTP::count(g,iiA,miA,niA,xiA,origWt[i]); }
			GTP::rewind(); for (size_t i=0;i<FldCTs.size();++i) { char g=_read_trio(tr_ct[i],gt,chr_num,bp,DNonly,ndn,gpar); g_ct+=g; GTP::count(g,iiU,miU,niU,xiU,origWt[i+j]); }
			if (DN_max && ndn>DN_max) { DefaultJob.qcl.push_back("SKIPPED=MultiDeNovo"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } // ndn should not apply to ExAC
		}
		if ( iiA+iiU==0 ) { DefaultJob.qcl.push_back("SKIPPED=NoData"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; } // no data observed
		
		if ( (xiA==0&&xiU==0) || (xiA==niA&&xiU==niU) ) { if (ExtCT::ExAC_pfx.empty()) { DefaultJob.qcl.push_back("SKIPPED=NoVariation"); continue; } } // no variation observed
		
		// skip by MAC
		int xiT=xiA+xiU;
		int niT=niA+niU;
		int spMAC = std::min(xiT, niT-xiT);
		if (mac_ge && spMAC>=mac_ge) { DefaultJob.qcl.push_back("SKIPPED=MAC"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
		
		// skip by observed allele frequency
		if (perch::filSAF)
		{
			double SplAF = get_value(INFO,"SplAF");
			if (std::isnan(SplAF)) exit_error("cannot apply --filt-SplAF because SplAF is not calculated by vQC");
			if (SplAF>0.5) SplAF=1-SplAF;
			if (!perch::rf_SAF && SplAF> perch::filSAF) { DefaultJob.qcl.push_back("SKIPPED=SplAF"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
			if ( perch::rf_SAF && SplAF<=perch::filSAF) { DefaultJob.qcl.push_back("SKIPPED=SplAF_rev"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
		}
		
		// skip by founder allele frequency
		if (perch::filFAF)
		{
			double FdrAF = get_value(INFO,"AF_Founder");
			double FdrAN = get_value(INFO,"AN_Founder");
			if (!std::isnan(FdrAF) && !std::isnan(FdrAN))
			{
				if (FdrAF>0.5) FdrAF=1-FdrAF;
				if (FdrAN>=perch::filFAFminAN)
				{
					if (!perch::rf_FAF && FdrAF> perch::filFAF) { DefaultJob.qcl.push_back("SKIPPED=FdrAF"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
					if ( perch::rf_FAF && FdrAF<=perch::filFAF) { DefaultJob.qcl.push_back("SKIPPED=FdrAF_rev"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
				}
			}
		}
		
		// skip by founder allele frequency but keep common variant only (useful for QQ-plot)
		if (skip_if_maf_l)
		{
			double FdrAF = get_value(INFO,"AF_Founder");
			if (FdrAF>0.5) FdrAF=1-FdrAF;
			if (FdrAF<skip_if_maf_l) { DefaultJob.qcl.push_back("SKIPPED=UncommonVar"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
		}
		
		/*/ skip by missingness
		double MsgAll = (double)(miA+miU) / (iiA+iiU);
		double MsgCSs = iiA ? (double)miA / iiA : std::numeric_limits<double>::signaling_NaN();
		double MsgCTs = iiU ? (double)miU / iiU : std::numeric_limits<double>::signaling_NaN();
		if (perch::MisCut!=1)
		{
			// if (MsgCSs > perch::MisCut)	{ DefaultJob.qcl.push_back("SKIPPED=MissingRateCs"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
			// if (MsgCTs > perch::MisCut)	{ DefaultJob.qcl.push_back("SKIPPED=MissingRateCt"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
			if (MsgAll > perch::MisCut)	{ DefaultJob.qcl.push_back("SKIPPED=MissingRate"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
		}
		if (MisPVl)
		{
			double zScore = (MsgCSs-MsgCTs) / sqrt(MsgAll*(1-MsgAll)*(1.0/iiA+1.0/iiU));
			double pValue = cdf_norms_2sided_pv(zScore);
			if (pValue<=MisPVl) { DefaultJob.qcl.push_back("SKIPPED=MissingRateTest"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
		} //*/
		
		// skip by --cap. Both --cap and --top must be after all other "continue"s.
		if (max_pg)
		{
			if (DefaultJob.del.size()>=max_pg) { DefaultJob.qcl.push_back("SKIPPED=NumVarPerGene"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
			if (!DefaultJob.del.empty() && !std::isnan(BayesDel))
			{
				if (std::isnan(DefaultJob.del.back())) exit_error("in sorting Genotype File by BayesDel, NA should be treated as the minimum value.");
				else if (DefaultJob.del.back()<BayesDel) exit_error("the Genotype File is not sorted by BayesDel. See "+this_group);
			}
		}
		
		// swap genotype if ref allele is rare
		if (niA && niU)
		{
			PAF = xiA/niA * perch::preval + xiU/niU * (1-perch::preval);
		}
		else if (niU)
		{
			PAF = xiU/niU;
		}
		if (PAF>0.5 && !usePRS && !usePOR)
		{
			PAF=1-PAF;
			for (auto &ind:g_cs) GTP::swap_allele(ind);
			for (auto &ind:g_ct) GTP::swap_allele(ind);
			if (!FldInf.no_input())
			{
				INFO.push_back("vAAA_SWAP");
				in[FldInf[0]] = str_of_container(INFO,';');
			}
		}

		// skip by --top. Must be after swapping common variants.
		// count only the 1 or 2 most deleterious variants per sample, rely on genotype file being sorted by Del (descending) and MAF (ascending)
		// Do not erase genotypes, otherwise a variant may be counted for some samples but not the others, not fair. Instead, erase the variant.
		if (max_ps)
		{
			for (size_t i=0;i<FldCSs.size();++i) if (GTP::usable(g_cs[i])&&GTP::num_alt(g_cs[i])) { if (vc_cs[i]<max_ps) { ++vc_cs[i]; ++aiA; } /*else GTP::set_all_ref(g_cs[i]);*/ }
			for (size_t i=0;i<FldCTs.size();++i) if (GTP::usable(g_ct[i])&&GTP::num_alt(g_ct[i])) { if (vc_ct[i]<max_ps) { ++vc_ct[i]; ++aiU; } /*else GTP::set_all_ref(g_ct[i]);*/ }
			if (!aiA && !aiU) { DefaultJob.qcl.push_back("SKIPPED=NoVarAfterPruning"); ExtCT::ExAC_dta[GeneSymbol][VarIndex].clear(); continue; }
		}

		// Add variant to DefaultJob. After this, there should not be any "continue".
		// add structual variant
		if (!sv_inp.empty())
		{
			string prev_gsymb;
			if (!DefaultJob.gnm.empty()) prev_gsymb=DefaultJob.gnm.back();
			if (prev_gsymb!=GeneSymbol)
			{
				double sv_frq = AAA_SV::sv_exist(chr_num,GeneSymbol,record)/record.size();
				if (sv_frq && sv_frq<=fil_sv)
				{
					size_t j=FldCSs.size();
					for (size_t i=0;i<FldCSs.size();++i) { int ploidy=genepi::expected_ploidy(chr_num,bp,SexCSs[i]); char g=(ploidy==2 ? genepi::genotype::HoR : genepi::genotype::HpR); DefaultJob.gtM[i]  +=g; }
					for (size_t i=0;i<FldCTs.size();++i) { int ploidy=genepi::expected_ploidy(chr_num,bp,SexCTs[i]); char g=(ploidy==2 ? genepi::genotype::HoR : genepi::genotype::HpR); DefaultJob.gtM[i+j]+=g; }
					for (size_t i=0;i<FldCSs.size();++i) { int cp=AAA_SV::sv_copy(chr_num,GeneSymbol,record[i]  ); char g=itos(cp)[0]; DefaultJob.gtC[i]  +=g; }
					for (size_t i=0;i<FldCTs.size();++i) { int cp=AAA_SV::sv_copy(chr_num,GeneSymbol,record[i+j]); char g=itos(cp)[0]; DefaultJob.gtC[i+j]+=g; }
					DefaultJob.idx.push_back("SV_"+GeneSymbol);
					DefaultJob.wts.push_back(1);
					DefaultJob.maf.push_back(defMAF);
					DefaultJob.del.push_back(5);
					DefaultJob.vqs.push_back(5);
					DefaultJob.xiU.push_back(0);
					DefaultJob.niU.push_back(0);
					DefaultJob.xiA.push_back(0);
					DefaultJob.niA.push_back(0);
					DefaultJob.ndn.push_back(0);
					DefaultJob.gnm.push_back(GeneSymbol);
					DefaultJob.chr.push_back(chr_num);
					DefaultJob.pos.push_back(bp-1);
					DefaultJob.ref.push_back("SV_REF");
					DefaultJob.alt.push_back("SV_ALT");
					DefaultJob.rsk.push_back('A');
				}
			}
		}
		// add this variant
		vector<int> cgA = { 0,0,0 };	// count genotype 0/0 0/1 1/1 for affected
		vector<int> cgU = { 0,0,0 };	// count genotype 0/0 0/1 1/1 for unaffected
		{
			size_t j=FldCSs.size();
			for (size_t i=0;i<FldCSs.size();++i) { char& g=g_cs[i]; DefaultJob.gtM[i]  +=g; if (GTP::usable(g)) ++cgA[GTP::num_alt_recessive(g)]; }
			for (size_t i=0;i<FldCTs.size();++i) { char& g=g_ct[i]; DefaultJob.gtM[i+j]+=g; if (GTP::usable(g)) ++cgU[GTP::num_alt_recessive(g)]; }
			for (size_t i=0;i<FldCSs.size();++i) { int cp=AAA_SV::sv_copy(chr_num,GeneSymbol,record[i]  ,bp); char g=itos(cp)[0]; DefaultJob.gtC[i]  +=g; }
			for (size_t i=0;i<FldCTs.size();++i) { int cp=AAA_SV::sv_copy(chr_num,GeneSymbol,record[i+j],bp); char g=itos(cp)[0]; DefaultJob.gtC[i+j]+=g; }
		}
		double L10BF = 0; // BF in log10 scale
		if (Wt_VQS) { double v = (VQSLOD<0?VQSLOD:0);							 if (!std::isnan(v)) L10BF+=v; }
		if (Wt_DEL) { double v = BayesDel;										 if (!std::isnan(v)) L10BF+=v; }
		if (Wt_GBA) { double v = S2G[GeneSymbol];								 if (!std::isnan(v)) L10BF+=v; }
		double vw;
		if (!Var_Wt)					vw=1;
		else if (AnalysisMethod==SSU)	vw=pow(10,L10BF);
		else							vw=posterior_given_log10BF(L10BF);
		double maf=std::min(MaxAF,1-MaxAF); // prv PAF but not good for small studies (e.g, 15 cs 100 ct, prevalence 0.00002, a novel variant in ct would have PAF~0.01, not right)
		char rsk='U'; // effect allele. U=unknown A=alternative R=reference
		{
			string PRS_allele = get_string(INFO,"PRS_allele");
			if (!PRS_allele.empty())
			{
				if 		(PRS_allele=="ALT") rsk='A'; // maf=0;
				else if	(PRS_allele=="REF") rsk='R'; // maf=1;
				else exit_error("PRS_allele should be either ALT or REF.");
			}
		}

		if (usePRS || usePOR)
		{
			double PRS_beta = get_value(INFO,"PRS_beta"); // beta from publication
			double PRS_freq = get_value(INFO,"PRS_freq"); // freq from publication
			double loc_freq = std::numeric_limits<double>::signaling_NaN(); // freq from local population for PRS calibration and PRS_allele checking
			for (auto &id : perch::h_afID)
			{
				double x = get_value(INFO,id);
				if (!std::isnan(x))
				{
					if (rsk=='A')	loc_freq=x;
					else			loc_freq=1-x;
					break;
				}
			}
			if (!std::isnan(PRS_freq) && !std::isnan(loc_freq))
			{
				if ((PRS_freq>0.5&&loc_freq<0.5)||(PRS_freq<0.5&&loc_freq>0.5))
					if ((ref=="C"&&alt=="G")||(ref=="G"&&alt=="C")||(ref=="A"&&alt=="T")||(ref=="T"&&alt=="A"))
						lns<<showw<<"PRS_freq doesn't match with allele frequency in --af, please check "<<VarIndex<<flush_logger;
			}
			if (PRSPAF) { if (rsk=='A') maf=PAF; else maf=1-PAF; }
			else		maf=(!std::isnan(loc_freq)?loc_freq:PRS_freq);
			if (rsk=='U') 				exit_error("PRS_allele not exist or empty. This INFO field is required for --prs or --por.");
			if (std::isnan(PRS_beta))	exit_error("PRS_beta not exist or not a number. This INFO field is required for --prs or --por.");
			if (std::isnan(maf) && usePOR) exit_error("Allele frequency not known for "+VarIndex);
			vw=PRS_beta;
		}
		
		DefaultJob.idx.push_back(VarIndex);
		DefaultJob.wts.push_back(vw);
		DefaultJob.maf.push_back(maf);
		DefaultJob.del.push_back(BayesDel);
		DefaultJob.vqs.push_back(VQSLOD);
		DefaultJob.xiU.push_back(xiU);
		DefaultJob.niU.push_back(niU);
		DefaultJob.xiA.push_back(xiA);
		DefaultJob.niA.push_back(niA);
		DefaultJob.ndn.push_back(ndn);
		DefaultJob.gnm.push_back(GeneSymbol);
		DefaultJob.chr.push_back(chr_num);
		DefaultJob.pos.push_back(bp);
		DefaultJob.ref.push_back(ref);
		DefaultJob.alt.push_back(alt);
		DefaultJob.rsk.push_back(rsk);
		if (is_MHC) DefaultJob.MHC=true;
		
		// Add qcl
		if (showID)
		{
			vector<string> cs_het, cs_hom, ct_het, ct_hom;
			size_t j=FldCSs.size();
			for (size_t i=0;i<FldCSs.size();++i) if (GTP::usable(g_cs[i])) { int g=GTP::num_alt_recessive(g_cs[i]); if (g==1) cs_het.push_back(record[i]);   if (g==2) cs_hom.push_back(record[i]); }
			for (size_t i=0;i<FldCTs.size();++i) if (GTP::usable(g_ct[i])) { int g=GTP::num_alt_recessive(g_ct[i]); if (g==1) ct_het.push_back(record[i+j]); if (g==2) ct_hom.push_back(record[i+j]); }
			DefaultJob.qcl.push_back("het_cs:"+str_of_container(cs_het,',')+";hom_cs:"+str_of_container(cs_hom,',')+";het_ct:"+str_of_container(ct_het,',')+";hom_ct:"+str_of_container(ct_hom,','));
		}
		else if (showGT)
		{
			string this_var=in[0]+DLMTR+in[1]+DLMTR+in[2]+DLMTR+in[3]+DLMTR+in[4]+DLMTR+in[5]+DLMTR+in[6]+DLMTR+in[7]+DLMTR+GeneSymbol+DLMTR+in[FldFun[0]]+DLMTR+in[FldDet[0]];
			string this_qcl;
			size_t j=FldCSs.size();
			for (size_t i=0;i<FldCSs.size();++i) if (GTP::usable(g_cs[i]) && GTP::num_alt(g_cs[i])) { this_qcl+=(this_var+DLMTR+record[i]  +"="+gt[FldCSs[i]]+"\n"); }
			for (size_t i=0;i<FldCTs.size();++i) if (GTP::usable(g_ct[i]) && GTP::num_alt(g_ct[i])) { this_qcl+=(this_var+DLMTR+record[i+j]+"="+gt[FldCTs[i]]+"\n"); }
			if (this_qcl.empty()) this_qcl="SKIPPED=NoVariants";
			DefaultJob.qcl.push_back(this_qcl);
		}
		else if (showIF)
		{
			string this_qcl;
			if (is_vks)
			{
				string this_var=in[0]+DLMTR+in[1]+DLMTR+in[2]+DLMTR+in[3]+DLMTR+in[4]+DLMTR+in[5]+DLMTR+in[6]+DLMTR+in[7]+DLMTR+GeneSymbol+DLMTR+in[FldFun[0]]+DLMTR+in[FldDet[0]];
				size_t j=FldCSs.size();
				for (size_t i=0;i<FldCSs.size();++i) if (GTP::usable(g_cs[i])) AAA_IF::to_rpt(GeneSymbol,GTP::num_alt_recessive(g_cs[i]),i,  (this_var+DLMTR+record[i]  +"="+gt[FldCSs[i]]),this_qcl);
				for (size_t i=0;i<FldCTs.size();++i) if (GTP::usable(g_ct[i])) AAA_IF::to_rpt(GeneSymbol,GTP::num_alt_recessive(g_ct[i]),i+j,(this_var+DLMTR+record[i+j]+"="+gt[FldCTs[i]]),this_qcl);
				if (this_qcl.empty()) this_qcl="SKIPPED=NoIncidentalFindings";
			}
			else
			{
				this_qcl="SKIPPED=NotPathogenic";
			}
			DefaultJob.qcl.push_back(this_qcl);
		}
		else
		{
			DefaultJob.qcl.push_back(str_of_container(cgA,',')+"/"+str_of_container(cgU,','));
		}
	}
	write_results(AddInf,FldInf,FldRes,wrWhat,true);
	
	if (!permte && program.nt>1)
	{
		bjmJobs.run(_run_job,program.nt);
		program.nt=1;
		var_format.forbid_nf_rpt();
		goto read_var_group;
	}
	else if (permte)
	{
		if (HasAll)
		{
			for (auto &repeat:wjp_res_db) wjp_res_mx.push_back(repeat.second.get(STAT::MAX));
			permte=0;
			program.nt=1;
			var_format.forbid_nf_rpt();
			goto read_var_group;
		}
		else
		{
			for (auto &repeat:wjp_res_db) program.outf<<repeat.second.get(STAT::MAX)<<endl;
		}
	}
	else if (AnalysisMethod==OUT)
	{
		if (AllDataSym.empty()) exit_error("no data. possible reasons include no cases or no variants");
		program.outf<<"SeqID\tAff";
		for (size_t c=0;c<CovVar[0].size();++c) program.outf<<"\tcov"<<c+1;
		for (size_t g=0;g<AllDataSym.size();++g) program.outf<<'\t'<<AllDataSym[g];
		program.outf<<"\tweight"<<endl;
		for (size_t i=0;i<DepVar.size();++i)
		{
			bool has_missing=false;
			for (size_t g=0;g<AllDataSym.size();++g) if (std::isnan(AllDataGtp[i][g])) { has_missing=true; break; }
			if (has_missing) continue;
			
			program.outf<<record[i]<<'\t'<<DepVar[i];
			for (size_t c=0;c<CovVar[i].size();++c) program.outf<<'\t'<<CovVar[i][c];
			for (size_t g=0;g<AllDataSym.size();++g) program.outf<<'\t'<<AllDataGtp[i][g];
			program.outf<<'\t'<<sqrtWt[i]<<endl;
		}
	}
	
	return 0;
}
