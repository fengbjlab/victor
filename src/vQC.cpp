/*
 Options not in help text:
 --mem S           Memory quota {_Default_memory}
 --qq F            Draw a QQ plot of allelic association test in common variants to file F (.pdf) {}
 --maf-for-qq D    In drawing the QQ plot, use variants whose MAF >= D {_Default_qq_maf}
 --excl-for-qq F   In drawing the QQ plot, exclude regions defined in F {_Default_qq_xcl}
 --qc-gtp B        Rewrite genotypes after QC {_Default_qc_gtp}
 --qc-mnp-only B   QC MNP variants only {_Default_qc_mnp}
 --qc-prv B        QC variants that have passed previous vQC {_Default_qc_prv}
 --x-distance I    Basepairs between X markers for sex inference {_Default_x_distance}
 --no-filt-known B  Do not exclude known variants (MaxAF>0) unless otherwise filtered by the above 5 options (not recommended) {_Default_noFtKn}
 --no-filt-bad B    Do not exclude bad variants (--filt-xx) unless otherwise filtered by the above 5 options (not recommended) {_Default_noFtBad}
 --qq F             Draw a QQ plot to file F (something.pdf, will also save p-values to something_pv.txt) {_Default_qq_file}
 --maf-for-qq D     In drawing the QQ plot, use the variants with MAF >= D {_Default_qq_maf}
 --excl-for-qq F    In drawing the QQ plot, exclude regions defined in file F {_Default_qq_xcl}
 --cadd-cutoff D      unless CADD raw score is >= D (NaN to turn off) {_Default_cadd_cutoff} or
 --fmnc-cutoff D      unless fathmm-MKL non-coding is >= D (NaN to turn off) {_Default_fmnc_cutoff}
 --qc-ploidy B       Correct genotypes according to the expected ploidy of chrX, chrY, chrM {_Default_qc_pld}

 Todo:
 --filt-info
 */
#include <cinttypes>
#include <tft/libfbj_base.hpp>
#include <tft/libfbj_file.hpp>
#include <tft/libfbj_math.hpp>
#include <tft/libfbj_genepi.hpp>
#include <tft/libfbj_program.hpp>
#include <tft/libfbj_mtjobs.hpp>
#include <tft/libfbj_regress.hpp>
#include "victor_par.hpp"

using namespace std;
typedef genepi::genotype GTP;

// HWE test, returns p-value, modified from http://csg.sph.umich.edu/abecasis/Exact/snp_hwe.c
double my_SNPHWE(const int obs_hets, const int obs_hom1, const int obs_hom2, vector<double>& het_probs)
{
	int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;
	int rare_copies = 2 * obs_homr + obs_hets;
	int genotypes   = obs_hets + obs_homc + obs_homr;
	if (genotypes==0 || rare_copies==0) return 1;
	
	het_probs.assign (rare_copies+1, 0);
	
	/* start at midpoint */
	// int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
	// the above is not robust to large sample size (e.g., Het=25768, Hom=22444, AC=70656, AN=111176). So change to below
	int mid = rare_copies * ( 1.0 - rare_copies / (2.0 * genotypes) );
	
	/* check to ensure that midpoint and rare alleles have same parity */
	if ((rare_copies & 1) ^ (mid & 1)) mid++;
	
	int curr_hets = mid;
	int curr_homr = (rare_copies - mid) / 2;
	int curr_homc = genotypes - curr_hets - curr_homr;
	
	
	het_probs[mid] = 1.0;
	double sum = het_probs[mid];
	for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
	{
		het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
		sum += het_probs[curr_hets - 2];
		
		/* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
		++curr_homr;
		++curr_homc;
	}
	
	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;
	for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
	{
		het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0));
		sum += het_probs[curr_hets + 2];
		
		/* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
		--curr_homr;
		--curr_homc;
	}
	
	for (int i = 0; i <= rare_copies; ++i) het_probs[i] /= sum;
	
	double p_hwe = 0.0;
	for (int i = 0; i <= rare_copies; ++i)
	{
		if (het_probs[i] > het_probs[obs_hets]) continue;
		p_hwe += het_probs[i];
	}
	p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;
	return p_hwe;
}

// parameters
struct parameters {
	string			spl_in;	// sample file for QC
	string			coh_in;	// sample cohort
	string			ped_in;	// pedigree file in LINKAGE format
	string			dup_in;	// duplicated sample file for QC
	vector<string>	rep_in;	// replication by SNP array: sample file, genotype file
	string			h_dels;	// header of dellbf defined by user
	string			spl_qc; // sample-wise QC output filename
	bool			rm_var; // removed QC not passed variants
	bool			qc_prv; // QC variants labeled VICTOR_QC=PASS
	bool			skipSV;	// exclude structural variations
	bool			SNonly;	// exclude non-SNVs
	bool			IDonly;	// exclude non-InDels
	bool			AutoCh; // exclude non-autosomal chrosomes
	bool			qc_gtp;	// rewrite genotypes after QC
	bool			qc_mnp;	// rewrite genotypes after QC, for MNP only
	bool			qc_pld;	// fix ploidy for chrX, chrY, chrM
	bool			noFtKn;	// do not filter known variants (MaxAF!=0 or ID=rs#)
	bool			noFtAl; // do not filter any variants
	bool			ftNoGt; // filter variants if NoGenotypes
	bool			ftNoVr;	// filter variants if NoVariation
	bool			hwesnf;	// skip non-founders when doing HWE test in unknown samples
	bool			hweict; // do HWE test in controls if population is not defined
	double			ftAFNG;	// filter variants by MaxAF (less than ftAFNG) for non-genic (intergenic / intronic) variants. 0.2 is good (informative, not too many).
	int				disDup; // filter variants by the number of discordances among duplicated samples (0 = no filter)
	int				FltMAC; // filter variants by minor allele acount (0 = no filter)
	int				MenErr; // filter variants by # Mendelian errors
	int				denovo; // filter variants by # de novo mutations. 2+ de novo of the same variant is unlikely.
	int				uno_dn; // filter variants by # non-de_novo if a de novo exist.
	int				min_AN; // filter variants by INFO::AN<min_AN (0 = no filter)
	double			hh_prp; // filter variants by proportion of heterozygous haploidy (0 = no filter)
	double			filRos; // filter variants by Rosenberg "Informativeness of Genetic Markers for Inference of Ancestry". 0 = no filter.
	double			i2info;	// filter variants by IMPUTE2 INFO score (>0.7)
	double			MR_pvl;	// skip variants if missing rate cs-ct w p-value < MR_pvl, 0=noQC.
	double			HW_pvl;	// skip variants if HWE test in controls p-value < HW_pvl. 0=noQC.
	double			AF_pvl;	// skip variants if prob(obs) in samples p-value < AF_pvl. 0=noQC.
	double			AF_dif;	// skip variants if AF between Cs-Ct has p-value < AF_dif. 0=noQC.
	double			AF_cut;	// allele frequency cutoff for test of probability of seeing a rare variant
	double			sum_RV; // definition of rare variant. MAF is calculated from MaxAF and data (if sample size is big).
	double			cov_pc; // well covered regions - percent of samples
	int				cov_dp; // well covered regions - DP cutoff
	int				mnAvDP; // min averaged DP
	int				mxAvDP; // max averaged DP
	int				TotSpl; // total number of samples for sum_RV
	int				y_n_mn;	// Y call numb cutoff for man
	double			y_r_mn;	// Y call rate cutoff for man
	int				y_n_wm;	// Y numb var cutoff for woman
	double			y_r_wm;	// Y call rate cutoff for woman
	double			xSxLOD; // X sex LOD score threshold
	double			x_err1; // heterozygous haploid rate
	double			x_err2; // het->hom error rate. To make things simple, I use x_err1 instead.
	double			x_err3; // hom->het error rate. To make things simple, I use x_err1 instead.
	int				x_dist; // basepairs between X markers for sex inference
	string			memcap;	// memory capacity, must ends with G/g/M/m
	vector<string>	InfoHW; // INFO field for Hardy Weinberg test
	string			InfoAC;	// INFO field for AC
	string			InfoAN;	// INFO field for AN
	parameters() {
		rm_var = true;				// removed QC not passed variants
		skipSV = false;				// exclude structural variations
		SNonly = false;				// exclude non-SNVs
		IDonly = false;				// exclude non-InDels
		AutoCh = false;				// exclude non-autosomal chrosomes
		qc_prv = true;				// QC variants labeled VICTOR_QC=PASS
		qc_gtp = true;				// rewrite genotypes after QC
		qc_mnp = false;				// rewrite genotypes after QC, for MNP only
		qc_pld = true;				// fix ploidy for chrX, chrY, chrM
		noFtKn = false;				// do not filter known variants (MaxAF!=0)
		noFtAl = false;				// do not filter any variants
		ftNoGt = true;				// filter variants if NoGenotypes
		ftNoVr = true;				// filter variants if NoVariation
		hwesnf = true;				// skip non-founders when doing HWE test in unknown samples
		hweict = false;				// do HWE test in controls if population is not defined
		ftAFNG = 0;					// filter variants by MaxAF (less than ftAFNG) for non-genic (intergenic / intronic) variants. 0.2 is good (informative, not too many).
		disDup = 1;					// number of discrepancies among duplicated samples to filter variants (0 = no filter)
		FltMAC = 0;					// filter variants by minor allele acount (0 = no filter)
		MenErr = 1;					// filter variants by # Mendelian errors excluding de novo mutations.
		denovo = 2;					// filter variants by # de novo mutations. 2+ de novo of the same variant is unlikely.
		uno_dn = 0;					// filter variants by # non-de_novo if a de novo exist.
		min_AN = 0;					// filter variants by INFO::AN<min_AN (0 = no filter)
		hh_prp = 0.1;				// filter variants by proportion of heterozygous haploidy (0 = no filter). In RGC, 328/341 variants have 0 HH, others in [0.000592417,0.9869] evenly distributed.
		filRos = 0;					// filter variants by Rosenberg "Informativeness of Genetic Markers for Inference of Ancestry". 0 = no filter.
		i2info = 0.7;				// filter variants by IMPUTE2 INFO score (>0.7)
		MR_pvl = 0.000001;			// skip variants if missing rate cs-ct w p-value < MR_pvl, 0=noQC.
		HW_pvl = 0.000001;			// skip variants if HWE test in controls p-value < HW_pvl. 0=noQC.
		AF_pvl = 0;					// skip variants if prob(obs) in samples p-value < AF_pvl. 0=noQC.
		AF_dif = 0;					// skip variants if AF between Cs-Ct has p-value < AF_dif. 0=noQC.
		AF_cut = 0.01;				// allele frequency cutoff for the test of probability of seeing a rare variant
		sum_RV = 0.05;				// def of rare var (use 0.05 to be robust to small sample size, the prob of no obs in 200 chr for MAF=0.05 is 0.000035). 0=no summary. Required by sample-wise QC
		cov_pc = 0;					// well covered regions - percent of samples. Previously 0.8
		cov_dp = 0;					// well covered regions - DP cutoff. Previously 10
		mnAvDP = 0;					// min averaged DP. Previously 10.
		mxAvDP = 0;					// max averaged DP. Previously 500.
		TotSpl = 0;					// total number of samples for sum_RV
		y_n_mn = 5;					// Y call numb cutoff for man.
		y_r_mn = 0.8;				// Y call rate cutoff for man.   In RGC, Ycall% in mn = 0.99714.    Ycall% in wm range [0,0.57].
		y_n_wm = 100;				// Y numb var cutoff for woman
		y_r_wm = 0.02;				// Y call rate cutoff for woman. In RGC, Ycall% in wm = 0.00553145. Ycall% in mn range [0.40,1]. Zero means don't use this. I have seen data (Kathleen Cooney project 1) with Ycall%=1/46=0.02 in men.
		xSxLOD = 3;					// X sex LOD score threshold
		x_err1 = 0.0009;			// heterozygous haploid rate. 3.89465e-05 in the true males (removing contaminated / sex wrong samples) of RGC. Better to use a bigger number.
		x_err2 = 0.0006;			// dizygotic het->hom error rate. This number was estimated by comparing ExomeXX with UGP/NPF-WES. To make things simple, I use x_err1 instead.
		x_err3 = 0.0002;			// dizygotic hom->het error rate. This number was estimated by comparing ExomeXX with UGP/NPF-WES. To make things simple, I use x_err1 instead.
		x_dist = 1000;				// basepairs between X markers for sex inference
		memcap = "200m";			// memory capacity for each squad, must ends with G/g/M/m
		InfoAC = "AC";				// ExAC variable name to find alternative allele count, should ends with =
		InfoAN = "AN";				// ExAC variable name to find chromosome count, should ends with =
	}
} par;

// stable data

genepi::ChrRegions min_bed;

struct DupGrp {
	vector<string>	IDs; // SeqID, not used after reading dup_in/rep_in
	vector<int>		FNs; // 0-based field number, 0 is allowed, -1 is not set
	vector<int>		SXs; // 0/1/2 sex
	vector<int>		Spl; // coordinate in samples
};
vector<DupGrp>	DupSpl;
vector<DupGrp>	RepSpl;

string						dnv_in; // de novo mutation input file name
map<string, vector<int> >	dnv_db; // de novo mutation database

struct Sample {
	// input
	int		sex;
	int		aff;
	string	ID;
	int		fld; // 0-based field number, 0 is allowed, -1 is not set
	bool	fdr; // is founder
	// statistics
	int mss; // missing
	int nms; // non-missing
	int het; // heterozygous
	int aho; // non-reference homozygous
	int rhe; // heterozygous RV
	int rho; // homozygous RV
	int nes; // novel heterozygous singletons
	int nos; // novel homozygous singletons
	int oes; // observed heterozygous singletons
	int oos; // observed homozygous singletons
	int eti; // (protein-coding) exon transition
	int etv; // (protein-coding) exon transversion
	int msy; // missing Y genotypes
	int nmy; // non-missing Y genotypes
	int n1y; // ref/alt Y genotypes
	int n2y; // alt/alt Y genotypes
	int msx; // missing X genotypes
	int nmx; // non-missing X genotypes used for sex inference
	int n1x; // ref/alt X genotypes
	int n2x; // alt/alt X genotypes
	int rnm; // replication not missing
	int rbg; // replication both genotyped
	int rcc; // replication concordance
	int nGQ; // number of GQ
	int nDP; // number of DP
	int nGP; // number of GP
	double tGQ; // mean GQ
	double tDP; // mean DP
	double tGP; // mean GP
	double lrx;
	Sample(int i_sex, int i_aff, string& i_id, int i_fld, bool i_fdr):fld(-1),mss(0),nms(0),het(0),aho(0),rhe(0),rho(0),nes(0),nos(0),oes(0),oos(0),eti(0),etv(0),msy(0),nmy(0),n1y(0),n2y(0),msx(0),nmx(0),n1x(0),n2x(0),rnm(0),rbg(0),rcc(0),nGQ(0),nDP(0),nGP(0),tGQ(0),tDP(0),tGP(0),lrx(0) { sex=i_sex; aff=i_aff; ID=i_id; fld=i_fld; fdr=i_fdr; }
};
vector<Sample>				samples;
map<string, vector<int> >	PopLoc;		// PopLoc[pop][] = coordinate in samples[]
map<string, vector<int> >	CohLoc;		// CohLoc[coh][] = coordinate in samples[]
set<int>		chrX_bp;
std::mutex		spl_mt;

// pedigrees
struct nuclear_ped {
	vector<string>	id;		// DAD MOM CH1 CH2 ..
	vector<int>		sex;	// DAD MOM CH1 CH2 ..
	vector<int>		col;	// DAD MOM CH1 CH2 .. 0-based
};
vector<nuclear_ped> pedigrees;	// at least one parent sequenced

field_numbers	FldChr(false,true);	// field numb for #CHROM
field_numbers	FldPos(false,true);	// field numb for POS
field_numbers	FldSNP(false,true);	// field numb for ID
field_numbers	FldRef(false,true);	// field numb for REF
field_numbers	FldAlt(false,true);	// field numb for ALT
field_numbers	FldFlt(false,true);	// field numb for FILTER
field_numbers	FldInf(false,true);	// field numb for INFO
field_numbers	FldFmt(false,true);	// field numb for FORMAT
field_numbers	FldXAF(false,true);	// field numb for MaxAF

// row data. In reading VCF, write to rsrc & rows. In calculation, create flag and modify rows. In writing output, print rows and clear everything.
const int				kNSq=2;			// number of squads
vector<int_fast64_t>	rsrc(kNSq,0);	// resource (memory) taken
deque< vector<string> > rows[kNSq];		// data for in and out. size() does not need to be N times of num_threads.
deque< int >			flag[kNSq];		// whether to print to program.outf (0/1)
std::mutex				iome[kNSq];		// input/output mutex
field_numbers			ofld(false,true);
int						CADD_column = -2; // -2 = not annotated; -1 = in INFO; 0+ = column.
int						FMNC_column = -2; // -2 = not annotated; -1 = in INFO; 0+ = column.
int						ColDel = -2; // -2 = not annotated; -1 = in INFO; 0+ = column.
int 					ColSEG = -2; // column number for vSEG analysis results as log10 Bayes factors

// outputs
string								log_fn;			// file name, set by program option
boost::iostreams::filtering_ostream logout;			// file, set by main if !log_fn.empty()
std::mutex							log_mt;			// mutex

genepi::ChrRegions					qq_xbed;
double								qq_maf=0.01;	// maf cutoff for QQ
string								qq_xcl;			// exclude region for QQ
string								pqq_fn;			// file name, set by program option
bool								savePV=true;	// save p-values to file
multiset<double>					pqq_pv;			// data
stringstream						pqq_log;		// output
std::mutex							pqq_mt;			// mutex

int									fnd_RV[3]={0,0,0};	// number of rare variants     per founder, control, case
int									fnd_SV[3]={0,0,0};	// number of singleton variant per founder, control, case
int									fnd_PV[3]={0,0,0};	// number of personal variant  per founder, control, case
int									fnd_RS[3]={0,0,0};	// number of rare SNVs         per founder, control, case
int									fnd_SS[3]={0,0,0};	// number of singleton SNVs    per founder, control, case
int									fnd_PS[3]={0,0,0};	// number of personal SNVs     per founder, control, case
std::mutex							nrv_mt;				// mutex

string								dup_fn;
vector<int>							concord_n;
vector<int>							concord_p;
vector<int>							false_pos;
vector<int>							false_neg;
vector<int>							het_2_hom;
vector<int>							hom_2_het;
vector<int>							SplSet1na;
vector<int>							soft_miss;
vector<int>							hard_miss;
vector<int>							both12_na;
boost::iostreams::filtering_ostream	GQ_DP_log;
std::mutex							gtp_mt;			// mutex
void test_dup(const vector<char>& GenoVec, const vector<int>& GQ, const vector<int>& DP)
{
	if (!par.rep_in.empty() && !par.dup_in.empty()) return; // can't let two analyses write to the same data above
	if (GTP::valid(GenoVec[0]))
	{
		for (size_t i=1;i<GenoVec.size();++i)
		{
			if (GTP::valid(GenoVec[i]))
			{
				if (GTP::usable(GenoVec[0]))
				{
					if (GTP::usable(GenoVec[i]))
					{
						int g1=GTP::num_alt(GenoVec[0]);
						int g2=GTP::num_alt(GenoVec[i]);
						if (g1==g2)
						{
							if (g1==0)			{ gtp_mt.lock(); ++concord_n[i]; if (!dup_fn.empty()&&GQ[i]>=0&&DP[i]>=0) GQ_DP_log<<"CNeg\t"<<GQ[i]<<'\t'<<DP[i]<<endl; gtp_mt.unlock(); }
							else				{ gtp_mt.lock(); ++concord_p[i]; if (!dup_fn.empty()&&GQ[i]>=0&&DP[i]>=0) GQ_DP_log<<"CPos\t"<<GQ[i]<<'\t'<<DP[i]<<endl; gtp_mt.unlock(); }
						}
						else
						{
							if (g1==0)			{ gtp_mt.lock(); ++false_pos[i]; if (!dup_fn.empty()&&GQ[i]>=0&&DP[i]>=0) GQ_DP_log<<"DPos\t"<<GQ[i]<<'\t'<<DP[i]<<endl; gtp_mt.unlock(); }
							else
							{
								if		(g2==0) { gtp_mt.lock(); ++false_neg[i]; if (!dup_fn.empty()&&GQ[i]>=0&&DP[i]>=0) GQ_DP_log<<"DNeg\t"<<GQ[i]<<'\t'<<DP[i]<<endl; gtp_mt.unlock(); }
								else if (g2==2) { gtp_mt.lock(); ++het_2_hom[i]; if (!dup_fn.empty()&&GQ[i]>=0&&DP[i]>=0) GQ_DP_log<<"DHom\t"<<GQ[i]<<'\t'<<DP[i]<<endl; gtp_mt.unlock(); }
								else			{ gtp_mt.lock(); ++hom_2_het[i]; if (!dup_fn.empty()&&GQ[i]>=0&&DP[i]>=0) GQ_DP_log<<"DHet\t"<<GQ[i]<<'\t'<<DP[i]<<endl; gtp_mt.unlock(); }
							}
						}
					}
					else
					{
						if (GTP::num_alt(GenoVec[0])>0) { gtp_mt.lock(); ++hard_miss[i]; gtp_mt.unlock(); }
						else							{ gtp_mt.lock(); ++soft_miss[i]; gtp_mt.unlock(); }
					}
				}
				else
				{
					if (GTP::usable(GenoVec[i]))	{ gtp_mt.lock(); ++SplSet1na[i]; gtp_mt.unlock(); }
					else							{ gtp_mt.lock(); ++both12_na[i]; gtp_mt.unlock(); }
				}
			}
		}
	}
}

string								vqs_fn;			// file name
stringstream						vqs_log;		// data
std::mutex							vqs_mt;			// mutex

string								dnv_fn;			// file name
stringstream						dnv_log;		// data
std::mutex							dnv_mt;			// mutex

void do_log(string& reasons, const string& idx, const int code)
{
	static const vector<string> message = {
		"PASS",
		"VariantType",
		"StructuralVar",
		"UnknChr|nonAutoChr",
		"ChrRegion",
		"ANis0",
		"ACis0",
		"MinorAlleleCount",
		"TotalAlleleCount",
		"MissingRateInINFO",
		"HWD_InINFO",
		"vSEG_Error",
		"MaxAF|SplAF|PopAF|FdrAF",
		"PreviousQC",
		"VQSLOD_missing",
		"VQSLOD_low",
		"HardFilter",
		"FILTER",
		"IMPUTE2_INFO",
		"DiscordanceInDup",
		"DiscordanceInRep",
		"min_avg_DP_ct",
		"min_avg_DP_cs",
		"max_avg_DP_ct",
		"max_avg_DP_cs",
		"ExcessiveObs|FreqDiff",
		"MinorAlleleCount",
		"PctCovered_ct",
		"PctCovered_cs",
		"MissingRate",
		"Missing_Cs>Ct",
		"Missing_Ct>Cs",
		"NoGenotype",
		"NoVariation",
		"HetHaploidy",
		"HWD_InControls",
		"HWD_InOnePop",
		"MendelianInconsistency",
		"MultipleDeNovo",
		"UnoriginalDeNovo",
		"NonGenicUncommonLowDel",
		"NonAIMs"
	};
	static const vector<int> qc_filter  = {
		elog.get_token(message[0]),
		elog.get_token(message[1]),
		elog.get_token(message[2]),
		elog.get_token(message[3]),
		elog.get_token(message[4]),
		elog.get_token(message[5]),
		elog.get_token(message[6]),
		elog.get_token(message[7]),
		elog.get_token(message[8]),
		elog.get_token(message[9]),
		elog.get_token(message[10]),
		elog.get_token(message[11]),
		elog.get_token(message[12]),
		elog.get_token(message[13]),
		elog.get_token(message[14]),
		elog.get_token(message[15]),
		elog.get_token(message[16]),
		elog.get_token(message[17]),
		elog.get_token(message[18]),
		elog.get_token(message[19]),
		elog.get_token(message[20]),
		elog.get_token(message[21]),
		elog.get_token(message[22]),
		elog.get_token(message[23]),
		elog.get_token(message[24]),
		elog.get_token(message[25]),
		elog.get_token(message[26]),
		elog.get_token(message[27]),
		elog.get_token(message[28]),
		elog.get_token(message[29]),
		elog.get_token(message[30]),
		elog.get_token(message[31]),
		elog.get_token(message[32]),
		elog.get_token(message[33]),
		elog.get_token(message[34]),
		elog.get_token(message[35]),
		elog.get_token(message[36]),
		elog.get_token(message[37]),
		elog.get_token(message[38]),
		elog.get_token(message[39]),
		elog.get_token(message[40]),
		elog.get_token(message[41])
	};
	elog.add(qc_filter[code]);
	if (!log_fn.empty() && code)
	{
		log_mt.lock();
		logout<<idx<<'\t'<<message[code]<<endl;
		log_mt.unlock();
	}
	if (reasons.empty())	reasons =message[code];
	else					reasons+=("/"+message[code]);
}

// job data
struct JobData {
	size_t sqd_id;	// input: one of the kNSq squads, 0 to kNSq-1
	size_t row_id;	// input: coordinate, 0 to num_threads-1
	size_t skip_n;	// input: jump over num_threads rows
	vector<double> het_probs; // for HWE
	JobData(size_t j, size_t r, size_t s) { sqd_id=j; row_id=r; skip_n=s; }
};

void compute(JobData& p)
{
	for (size_t section=0; section+p.row_id<rows[p.sqd_id].size(); section+=p.skip_n)
	{
		vector<string>& in = rows[p.sqd_id][section+p.row_id];
		
		// check file
		string this_index = in[FldChr[0]]+'\t'+in[FldPos[0]]+'\t'+in[FldRef[0]]+'\t'+in[FldAlt[0]];
		if (in[FldAlt[0]].find(',')!=std::string::npos)
		{
			string message="# The input Genotype File has not been split alternative alleles. Culprit: "+this_index;
			lns<<showe<<message<<fatal;
			throw message;
		}
		
		// insert de novo
		for (map<string, vector<int> >::const_iterator it=dnv_db.find(this_index); it!=dnv_db.end(); )
		{
			for (auto &i:dnv_db[this_index]) GTP::to_dn(in[i]);
			break;
		}
		
		string reasons;
		
		// skip by var type
		string& ref = in[FldRef[0]];
		string& alt = in[FldAlt[0]];
		bool is_snv = ( (ref=="A" || ref=="T" || ref=="C" || ref=="G") && (alt=="A" || alt=="T" || alt=="C" || alt=="G") );
		if (par.SNonly && !is_snv)									{ do_log(reasons,this_index,1); continue; }
		if (par.IDonly &&  is_snv)									{ do_log(reasons,this_index,1); continue; }
		if (par.skipSV && str_has(alt,"<"))							{ do_log(reasons,this_index,2); continue; }
		
		// skip by chr
		int chr_num = genepi::read_chr_num(in[FldChr[0]]);
		if (chr_num<1) { do_log(reasons,this_index,3); continue; }
		bool is_auto = genepi::is_autosomal(chr_num);
		bool is_X = genepi::is_chrX(chr_num);
		bool is_Y = genepi::is_chrY(chr_num);
		bool is_M = genepi::is_chrM(chr_num);
		if (par.AutoCh && !is_auto) { do_log(reasons,this_index,3); continue; }
		
		// skip by chr region
		int bp=-1; try { bp = boost::lexical_cast<int>(in[FldPos[0]]); } catch (...) { string message="Failed to read "+in[FldPos[0]]+" as a position in basepairs."; lns<<showe<<message<<fatal; throw message; }
		if (!perch::within_covered_region(chr_num,bp)) { do_log(reasons,this_index,4); continue; }

		// get INFO
		vector<string> INFO;
		if (!FldInf.no_input()) { if (!in[FldInf[0]].empty() && in[FldInf[0]]!=".") boost::split(INFO,in[FldInf[0]],boost::is_any_of(";")); }
		bool pr_qc = false; for (auto &i:INFO) if (i=="VICTOR_QC=PASS") { pr_qc=true; break; }
		bool isMNP = false; for (auto &i:INFO) if (i=="MNP_by_VANNER")  { isMNP=true; break; }
		if (!par.qc_prv && pr_qc) { flag[p.sqd_id][section+p.row_id]=1; continue; }
		if (par.qc_mnp && !isMNP) { flag[p.sqd_id][section+p.row_id]=1; continue; }
		double VQSLOD = get_value(INFO,"VQSLOD");
		double infoQD = get_value(INFO,"QD"); // use double so that it is nan if not exist
		double infoMQ = get_value(INFO,"MQ"); // nan automatically return false in comparison
		double infoFS = get_value(INFO,"FS"); // so no need to check isnan()
		double infoHS = get_value(INFO,"HaplotypeScore");
		double infoMR = get_value(INFO,"MQRankSum");
		double infoRP = get_value(INFO,"ReadPosRankSum");
		
		// get FORMAT
		genepi::gtp_par gpar;
		if (!FldFmt.no_input()) gpar.read(in[FldFmt[0]]);
		
		// calculate Ti/Tv for coding SNV only. Some coding SNV may be excluded due to "SpliceAltering" or "SpliceSite", but it's ok to exclude only a few.
		// The annotation is best suitable for vAnnGene --no-split. If overlapping genes are split into multiple lines, some variant may be counted twice. But that's rare.
		string func_ann = get_string(INFO,perch::i_func);
		vector<string> func_vec;
		boost::split(func_vec,func_ann,boost::is_any_of(","));
		bool is_coding=false;
		for (size_t i=1; i<func_vec.size(); i+=3)
			if (func_vec[i]=="Synonymous"||func_vec[i]=="Missense"||func_vec[i]=="StopLoss"||func_vec[i]=="StopGain"||func_vec[i]=="Translation") { is_coding=true; break; }
		bool coding_snv = (is_snv && is_coding);
		bool is_ti = (coding_snv && ( (str_has("AG",ref)&&str_has("AG",alt)) || (str_has("CT",ref)&&str_has("CT",alt)) ));
		bool is_tv = (coding_snv && ( (str_has("AC",ref)&&str_has("AC",alt)) || (str_has("AT",ref)&&str_has("AT",alt)) || (str_has("GC",ref)&&str_has("GC",alt)) || (str_has("GT",ref)&&str_has("GT",alt)) ));
		
		// determine genic variant, robust to vAnnGene field not sorted by func_type
		bool is_ii = false;
		if (par.ftAFNG)
		{
			set<string> func_type;
			for (size_t i=1; i<func_vec.size(); i+=3)
				if (func_vec[i]!="Intergenic"&&func_vec[i]!="Intronic") func_type.insert(func_vec[i]);
			is_ii = func_type.empty();
		}
		
		// get MaxAF. Nan if not annotated, 0 if not match.
		double MaxAF = std::numeric_limits<double>::signaling_NaN();
		if (!FldXAF.no_input())
		{
			read_val(in[FldXAF[0]],MaxAF);
			if (std::isnan(MaxAF)) MaxAF=0;
		}
		else
		{
			for (auto& f:INFO)
				if (str_startsw(f,perch::h_MxAF+"="))
				{
					if (!read_val(f.substr(6),MaxAF)) MaxAF=0;
					break;
				}
		}
		if (std::isnan(MaxAF)) // not annotated
		{
			if (!vqs_fn.empty()) { string message="# --out-vqs requires the MaxAF annotation, but it is missing."; lns<<showe<<message<<fatal; throw message; }
			if (par.noFtKn) { string message="# --no-filt-known requires the MaxAF annotation, but it is missing."; lns<<showe<<message<<fatal; throw message; }
			if (perch::filXAF) { string message="# --filt-MaxAF requires the MaxAF annotation, but it is missing."; lns<<showe<<message<<fatal; throw message; }
			if (par.ftAFNG) { string message="# --rm-ii-MaxAF-lt requires the MaxAF annotation, but it is missing."; lns<<showe<<message<<fatal; throw message; }
			if (par.sum_RV) { string message="# summary of RV requires the MaxAF annotation, but it is missing."; lns<<showe<<message<<fatal; throw message; }
		}
		
		// get known status
		bool known = ( MaxAF>0 ); // || str_startsw(in[FldSNP[0]],"rs") );

		// decide whether to filter
		bool toFilt=par.rm_var;
		if (par.noFtAl) toFilt=false;
		if (par.noFtKn && known) toFilt=false;
		
		// filters based on AC and AN, applies to VCF w/o genotypes (ExAC)
		double MsAll = std::numeric_limits<double>::signaling_NaN(); // missing rate in all sequenced samples
		bool het_singleton = false;
		bool hom_singleton = false;
		if (par.TotSpl && !FldInf.no_input())
		{
			int INFO_AC=-1; get_int(INFO,par.InfoAC,INFO_AC);
			int INFO_AN=-1; get_int(INFO,par.InfoAN,INFO_AN);
			if (INFO_AC<0) exit_error(par.InfoAC+" sub-field missing or not a valid number (non-negative integer)");
			if (INFO_AN<0) exit_error(par.InfoAN+" sub-field missing or not a valid number (non-negative integer)");
			int SeqSpl=INFO_AN/2;
			int mac = std::min(INFO_AC, INFO_AN-INFO_AC);
			
			if (SeqSpl>par.TotSpl) exit_error("There are more sequenced samples than you spcified. Check the argument for --tot-spl or --info-an.");
			if (par.ftNoGt &&  INFO_AN==0)						{ do_log(reasons,this_index,5); if (toFilt) continue; } // NoGenotypes
			if (par.ftNoVr && (INFO_AC==0||INFO_AC==INFO_AN))	{ do_log(reasons,this_index,6); if (toFilt) continue; } // NoVariation
			if (mac==1) het_singleton=true;
			
			if (par.FltMAC)
			{
				if (mac<par.FltMAC) { do_log(reasons,this_index,7); if (toFilt) continue; }
			}
			if (par.min_AN)
			{
				if (INFO_AN<par.min_AN) { do_log(reasons,this_index,8); if (toFilt) continue; }
			}
			MsAll=1-(double)SeqSpl/par.TotSpl;
			if (perch::MisCut!=1)
			{
				if (MsAll>perch::MisCut) { do_log(reasons,this_index,9); if (toFilt) continue; }
			}
			if (par.HW_pvl)
			{
				double HWE = get_value(INFO,"HWE");
				if (!std::isnan(HWE) && HWE<=par.HW_pvl) { do_log(reasons,this_index,10); if (toFilt) continue; };
			}
		}

		// Hardy-Weinberg test with data in INFO
		if (is_auto && !par.InfoHW.empty())
		{
			bool HWD_INFO=false;
			for (size_t i=0;i<par.InfoHW.size();i+=3)
			{
				int AN, Het, Hom;
				bool good = (read_val_ge(get_string(INFO,par.InfoHW[i]),AN,2) &&
							 read_val_ge(get_string(INFO,par.InfoHW[i+1]),Het,0) &&
							 read_val_ge(get_string(INFO,par.InfoHW[i+2]),Hom,0));
				if (good)
				{
					double pValue = my_SNPHWE(Het,AN/2-Het-Hom,Hom,p.het_probs);
					if (pValue<=par.HW_pvl) { HWD_INFO=true; break; }
				}
			}
			if (HWD_INFO) { do_log(reasons,this_index,10); if (toFilt) continue; };
		}
		
		// skip if has Mendelian error.
		if (ColSEG!=-2)
		{
			double v=std::numeric_limits<double>::signaling_NaN();
			perch::read_variable(in,ColSEG,INFO,perch::h_SEGb,v);
			if (std::isnan(v)) { do_log(reasons,this_index,11); if (toFilt) continue; }
		}
		
		// skip by QC
		if (!FldInf.no_input())
		{
			if (str_has(in[FldInf[0]],"VICTOR_QC=BAD")) { do_log(reasons,this_index,13); if (toFilt) continue; }
		}
		
		// skip by VQSLOD or hard filtering
		if (!FldInf.no_input())
		{
			if (std::isnan(VQSLOD))
			{
				if (perch::VQSnan)	{ do_log(reasons,this_index,14); if (toFilt) continue; }
				if (perch::HFnoVQ || perch::hardft)
				{
					if (perch::FiltQD) { if (infoQD<perch::FiltQD)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
					if (perch::FiltMQ) { if (infoMQ<perch::FiltMQ)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
					if (perch::FiltFS) { if (infoFS>perch::FiltFS)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
					if (perch::FiltHS) { if (infoHS>perch::FiltHS)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
					if (perch::FiltMR) { if (infoMR<perch::FiltMR)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
					if (perch::FiltRP) { if (infoRP<perch::FiltRP)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
				}
			}
			else
			{
				if (vqs_fn.empty())
				{
					if (is_snv) { if (VQSLOD<perch::VQSsnv)		{ do_log(reasons,this_index,15); if (toFilt) continue; } }
					else		{ if (VQSLOD<perch::VQSidl)		{ do_log(reasons,this_index,15); if (toFilt) continue; } }
					if (perch::hardft)
					{
						if (perch::FiltQD) { if (infoQD<perch::FiltQD)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
						if (perch::FiltMQ) { if (infoMQ<perch::FiltMQ)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
						if (perch::FiltFS) { if (infoFS>perch::FiltFS)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
						if (perch::FiltHS) { if (infoHS>perch::FiltHS)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
						if (perch::FiltMR) { if (infoMR<perch::FiltMR)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
						if (perch::FiltRP) { if (infoRP<perch::FiltRP)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
					}
				}
			}
		}
		
		// skip by FILTER
		if (!FldFlt.no_input())
			if (!perch::filflt.empty() && !exist_element(perch::filflt,in[FldFlt[0]])) { do_log(reasons,this_index,17); if (toFilt) continue; }
		
		// skip by IMPUTE2 INFO
		if (par.i2info)
		{
			double i2info = get_value(INFO,"INFO");
			if (!std::isnan(i2info) && i2info<par.i2info) { do_log(reasons,this_index,18); if (toFilt) continue; }
		}
		
		// skip by BayesDel
		double BayesDel = std::numeric_limits<double>::signaling_NaN();
		if		(ColDel==-1)	BayesDel = get_value_sw(INFO,"BayesDel");
		else if (ColDel>=0)	read_val(in[ColDel],BayesDel);

		// duplications, put it before case-control so that genotypes are not filtered by GQ,DP,GP
		if (!DupSpl.empty())
		{
			int num_discrepancies=0;
			for (auto &DG:DupSpl)
			{
				set<char> GenoSet;
				vector<char> GenoVec;
				vector<int> GQ_Vec;
				vector<int> DP_Vec;
				for (size_t i=0;i<DG.FNs.size();++i)
				{
					char g=GTP::CNA;
					int GQ=-1,DP=-1; // initialization to -1 is not necessary, but do it anyway
					double GP=-1;
					if (DG.FNs[i]>=0)
					{
						if (i)	g=GTP::read_vcf(in[DG.FNs[i]],chr_num,bp,DG.SXs[i],false,GQ,DP,GP,gpar);
						else	g=GTP::read_vcf(in[DG.FNs[i]],chr_num,bp,DG.SXs[i],true,GQ,DP,GP,gpar);
						GTP::dephase(g);
					}
					if ( GTP::usable(g) && (GQ<0||GQ>=GTP::GQ_cut) && (DP<0||DP>=GTP::DP_cut) && (GP<0||GP>=GTP::GP_cut) ) GenoSet.insert(g);
					GenoVec.push_back(g);
					GQ_Vec.push_back(GQ);
					DP_Vec.push_back(DP);
				}
				if (GenoSet.size()>1) num_discrepancies+=(GenoSet.size()-1);
				
				// summarize genotypes
				test_dup(GenoVec,GQ_Vec,DP_Vec);
			}
			if (num_discrepancies>=par.disDup) { do_log(reasons,this_index,19); if (toFilt) continue; }
		}
		
		// replications, put it before case-control so that genotypes are not filtered by GQ,DP. This is useful for finding the cutoff values for GQ and DP.
		if (!RepSpl.empty() && !dup_fn.empty())
		{
			bool matched=false;
			vector<string>	rp;			// matched line
			string			rep_tbx;	// return by tabix
			string chr = genepi::convert_chr_num(chr_num);
			string chr_spe_gtp = boost::algorithm::replace_all_copy(par.rep_in[1],"@",chr);
			if (!FileExists(chr_spe_gtp+".tbi")) exit_error("cannot find "+chr_spe_gtp+".tbi");
			try { rep_tbx = exec("tabix "+chr_spe_gtp+" "+chr+":"+itos(bp)+"-"+itos(bp),false); }
			catch (const std::exception& error) { exit_error("tabix "+chr_spe_gtp+" "+chr+":"+itos(bp)+"-"+itos(bp)+" failed"); }
			if (!rep_tbx.empty())
			{
				rep_tbx.pop_back(); // rep_tbx ends with \n
				vector<string> rep_row;
				boost::split(rep_row,rep_tbx,boost::is_any_of("\n"));
				for (auto &row:rep_row)
				{
					boost::split(rp,row,boost::is_any_of("\t"));
					if (rp[3]==ref && rp[4]==alt && rp[1]==in[FldPos[0]]) { matched=true; break; }
				}
			}
			
			if (matched)
			{
				int num_discrepancies=0;
				vector<char> GenoVec(2);
				vector<int> GQ_Vec(2,-1); // initialization to -1 is not necessary, but do it anyway
				vector<int> DP_Vec(2,-1); // initialization to -1 is not necessary, but do it anyway
				vector<double> GP_Vec(2,-1); // initialization to -1 is not necessary, but do it anyway
				for (auto &DG:RepSpl)
				{
					GenoVec[0] = GTP::read_vcf(rp[DG.FNs[0]],chr_num,bp,DG.SXs[0],true, GQ_Vec[0],DP_Vec[0],GP_Vec[0],gpar);
					GenoVec[1] = GTP::read_vcf(in[DG.FNs[1]],chr_num,bp,DG.SXs[1],false,GQ_Vec[1],DP_Vec[1],GP_Vec[1],gpar);
					test_dup(GenoVec,GQ_Vec,DP_Vec);
					if ( GQ_Vec[1]>=0 && GQ_Vec[1]<GTP::GQ_cut ) GTP::set_missing(GenoVec[1]);
					if ( DP_Vec[1]>=0 && DP_Vec[1]<GTP::DP_cut ) GTP::set_missing(GenoVec[1]);
					if ( GP_Vec[1]>=0 && GP_Vec[1]<GTP::GP_cut ) GTP::set_missing(GenoVec[1]);
					if (GTP::usable(GenoVec[0])&&GTP::usable(GenoVec[1])&&GenoVec[0]!=GenoVec[1]) ++num_discrepancies;
				}
				if (num_discrepancies>=par.disDup) { do_log(reasons,this_index,20); if (toFilt) continue; }
			}
		}

		// QC by genotypes and summarize
		double MsRCS = std::numeric_limits<double>::signaling_NaN(); // missing rate in cs
		double MsRCT = std::numeric_limits<double>::signaling_NaN(); // missing rate in ct
		double SplAF = std::numeric_limits<double>::signaling_NaN();
		double FdrAF = std::numeric_limits<double>::signaling_NaN();
		double PopAF = std::numeric_limits<double>::signaling_NaN();
		double Cs_AF = std::numeric_limits<double>::signaling_NaN();
		double Ct_AF = std::numeric_limits<double>::signaling_NaN();
		double IntAF = MaxAF; // AF integrated MaxAF && ( PopAF || FdrAF )
		double Rosenberg = std::numeric_limits<double>::signaling_NaN(); // Rosenberg et al. "Informativeness of Genetic Markers for Inference of Ancestry"
		int    FdrAN = 0;
		int    FdrAC = 0;
		int    FdrHO = 0;
		int	   spMAC = -1;  // minor allele count among cs-ct    (no duplicated samples)
		int    osMAC = -1;	// minor allele count among overall  (no duplicated samples)
		int    fdMAC = -1;	// minor allele count among founders (no duplicated samples)
		vector<int> geno;	// coordinate is the same as samples. -1 or number of ALT.
		vector<int> r_GQ;	// record of GQ
		vector<int> r_DP;	// record of DP
		vector<int> r_GP;	// record of GP
		int ploidy_err=0;	// number of heterozygous haploidy
		int n_haploidy=0;	// number of haploidy individuals with genotype called
		double	xiA=0, xiU=0, xiX=0, xiF=0, xiT=0, xiO=0;	// count alternative alleles, for affected, unaffected, unknown, founders(in cs/ct/uk), respectively
		double  niA=0, niU=0, niX=0, niF=0, niT=0, niO=0;	// count non-missing alleles
		int		iiA=0, iiU=0, iiX=0, iiF=0, iiT=0, iiO=0;	// count valid samples, include missing
		int		miA=0, miU=0, miX=0, miF=0, miT=0, miO=0;	// count missing genotypes
		vector<int> cgA = { 0,0,0 };		// count genotype 0/0 0/1 1/1 for affected
		vector<int> cgU = { 0,0,0 };		// count genotype 0/0 0/1 1/1 for unaffected
		vector<int> cgX = { 0,0,0 };		// count genotype 0/0 0/1 1/1 for unknown
		vector<int> cgF = { 0,0,0 };		// count genotype 0/0 0/1 1/1 for founders
		int last_i=-1;	// sample coordinate associated with het_singleton / hom_singleton
		double testOA = std::numeric_limits<double>::signaling_NaN();
		double testOU = std::numeric_limits<double>::signaling_NaN();
		double testHW = std::numeric_limits<double>::signaling_NaN();
		double testMR = std::numeric_limits<double>::signaling_NaN();
		double csctGQ = std::numeric_limits<double>::signaling_NaN(); // ratio of mean ranks for GQ in cs/ct
		double csctDP = std::numeric_limits<double>::signaling_NaN(); // ratio of mean ranks for DP in cs/ct
		double psc[2] = {std::numeric_limits<double>::signaling_NaN(),std::numeric_limits<double>::signaling_NaN()}; // percent_of_samples_covered

		if (!samples.empty())
		{
			psc[0]=0; psc[1]=0;
			int totDP[2] = {0,0};
			int numDP[2] = {0,0};
			for (size_t i=0;i<samples.size();++i)
			{
				Sample& spl=samples[i];
				int ploidy = genepi::expected_ploidy(chr_num,bp,spl.sex);
				int GQ=-1;
				int DP=-1;
				double GP=-1;
				char g0 = GTP::read_vcf(in[spl.fld],true,GQ,DP,GP,gpar);// genotype original
				char g1 = GTP::fix_ploidy(g0,ploidy,ploidy_err);		// genotype corrected (ploidy)
				char gc = ( par.qc_pld ? g1 : g0 );						// genotype chosen for output.
				if (ploidy==1 && GTP::usable(g0)) ++n_haploidy;
				int a1 = (GTP::usable(gc) ? GTP::num_alt(gc) : -1);
				if		(spl.fdr)	 { GTP::count(gc,iiF,miF,niF,xiF); if (GTP::usable(gc)) ++cgF[a1]; if (a1>0) last_i=i; }
				if		(spl.aff==2) { GTP::count(gc,iiA,miA,niA,xiA); if (GTP::usable(gc)) ++cgA[a1]; if (par.qc_gtp) GTP::to_rewr(gc,in[spl.fld]); if (DP!=-1) { totDP[1]+=DP; ++numDP[1]; if (DP>=par.cov_dp) ++psc[1]; } }
				else if (spl.aff==1) { GTP::count(gc,iiU,miU,niU,xiU); if (GTP::usable(gc)) ++cgU[a1]; if (par.qc_gtp) GTP::to_rewr(gc,in[spl.fld]); if (DP!=-1) { totDP[0]+=DP; ++numDP[0]; if (DP>=par.cov_dp) ++psc[0]; } }
				else				 { GTP::count(gc,iiX,miX,niX,xiX); if (GTP::usable(gc)) ++cgX[a1]; if (par.qc_gtp) GTP::to_rewr(gc,in[spl.fld]); }
				
				if	(is_Y || is_X)
				{
					// Count geno on X,Y to infer sex. No fix_ploidy because sex unknown.
					int a0 = (GTP::usable(g0) ? GTP::num_alt(g0) : -1);
					geno.push_back(a0);
				}
				else
				{
					geno.push_back(a1);
				}
				r_GQ.push_back(GQ);
				r_DP.push_back(DP);
				r_GP.push_back(GP);
			}
			if (par.mnAvDP && numDP[0] && totDP[0]/numDP[0] < par.mnAvDP) { do_log(reasons,this_index,21); if (toFilt) continue; }
			if (par.mnAvDP && numDP[1] && totDP[1]/numDP[1] < par.mnAvDP) { do_log(reasons,this_index,22); if (toFilt) continue; }
			if (par.mxAvDP && numDP[0] && totDP[0]/numDP[0] >=par.mxAvDP) { do_log(reasons,this_index,23); if (toFilt) continue; }
			if (par.mxAvDP && numDP[1] && totDP[1]/numDP[1] >=par.mxAvDP) { do_log(reasons,this_index,24); if (toFilt) continue; }
			if (iiA && numDP[1]) psc[1]/=iiA; else psc[1]=std::numeric_limits<double>::signaling_NaN();
			if (iiU && numDP[0]) psc[0]/=iiU; else psc[0]=std::numeric_limits<double>::signaling_NaN();
			het_singleton = (cgF[1]==1 && cgF[2]==0);
			hom_singleton = (cgF[1]==0 && cgF[2]==1);
			
			xiT=xiA+xiU, xiO=xiT+xiX;
			niT=niA+niU, niO=niT+niX;
			iiT=iiA+iiU, iiO=iiT+iiX;
			miT=miA+miU, miO=miT+miX;
			if (niF) FdrAF = xiF / niF; // allele frequency in founders
			if (niT) SplAF = xiT / niT; // allele frequency in cases+controls
			if (niA) Cs_AF = xiA / niA; // allele frequency in cases
			if (niU) Ct_AF = xiU / niU; // allele frequency in controls
			FdrAN=niF;
			FdrAC=xiF;
			FdrHO=cgF[2];
			
			// skip by FdrAN
			if (par.min_AN)
			{
				if (FdrAN<par.min_AN) { do_log(reasons,this_index,8); if (toFilt) continue; }
			}
			
			if (niT)
			{
				if (niA && niU) // case-control study
				{
					double one_in_cs = perch::preval/niA;
					double one_in_ct = (1-perch::preval)/niU;
					PopAF = xiA * one_in_cs + xiU * one_in_ct; // Population allele frequency
					// If sample size is large enough (MAC=2 is RV), then calculate MAF from data
					// Use MAC=2 because homozygous singleton would still be counted as a rare variant.
					if (2*one_in_ct<=par.sum_RV && 2*one_in_cs<=par.sum_RV && (one_in_ct+one_in_cs)<=par.sum_RV)
					{
						if (!std::isnan(MaxAF)) { if (std::abs(0.5-PopAF)<std::abs(0.5-MaxAF)) IntAF=PopAF; } else IntAF=PopAF;
					}
				}
				else // case-only or control-only
				{
					if (2.0/niT<=par.sum_RV) // sample size is large enough (MAC=2 is RV)
					{
						if (!std::isnan(MaxAF)) { if (std::abs(0.5-FdrAF)<std::abs(0.5-MaxAF)) IntAF=FdrAF; } else IntAF=FdrAF;
					}
				}
				
				// skip by Prob obs. Actually not very effective. Disabled by default. Use AF test between Cs-Ct instead.
				if (par.AF_pvl)
				{
					// if (std::isnan(MaxAF)) { string message="# --filt-obs-pv requires the MaxAF annotation, but it is missing."; lns<<showe<<message<<fatal; throw message; }
					double af = (std::isnan(MaxAF) ? 0 : MaxAF);
					if		(af<=  par.AF_cut)	af=  par.AF_cut;
					else if (af>=1-par.AF_cut)	af=1-par.AF_cut;
					else						af=std::numeric_limits<double>::signaling_NaN();
					if (af<0.5)
					{
						if (niA && xiA/niA>  par.AF_cut && (testOA=cdf_binomial_ge(xiA,af,niA))<=par.AF_pvl) { do_log(reasons,this_index,25); if (toFilt) continue; }
						if (niU && xiU/niU>  par.AF_cut && (testOU=cdf_binomial_ge(xiU,af,niU))<=par.AF_pvl) { do_log(reasons,this_index,25); if (toFilt) continue; }
					}
					else if (af>0.5)
					{
						if (niA && xiA/niA<1-par.AF_cut && (testOA=cdf_binomial_le(xiA,af,niA))<=par.AF_pvl) { do_log(reasons,this_index,25); if (toFilt) continue; }
						if (niU && xiU/niU<1-par.AF_cut && (testOU=cdf_binomial_le(xiU,af,niU))<=par.AF_pvl) { do_log(reasons,this_index,25); if (toFilt) continue; }
					}
				}
				
				// AF test between Cs-Ct.
				if (par.AF_dif && niA && niU)
				{
					int table[4];
					table[0]=xiA;
					table[1]=niA-xiA;
					table[2]=xiU;
					table[3]=niU-xiU;
					double testAD = Fishers_exact_test_2x2(table,false);
					if (perch::_Debug) cerr<<table[0]<<' '<<table[1]<<' '<<table[2]<<' '<<table[3]<<' '<<testAD<<' '<<par.AF_dif<<endl;
					if (testAD<=par.AF_dif) { do_log(reasons,this_index,25); if (toFilt) continue; }
				}
			}
			
			// skip by MAC (minor allele count)
			if (niT) spMAC = std::min(xiT, niT-xiT);
			if (niO) osMAC = std::min(xiO, niO-xiO);
			if (niF) fdMAC = std::min(xiF, niF-xiF);
			if (par.FltMAC && osMAC!=-1)
			{
				if (osMAC<par.FltMAC) { do_log(reasons,this_index,26); if (toFilt) continue; }
			}
			
			// skip by batch effect (these measures are not good indicators, so by default they are not used)
			if (par.cov_pc && par.cov_dp && psc[0]<par.cov_pc) { do_log(reasons,this_index,27); if (toFilt) continue; }
			if (par.cov_pc && par.cov_dp && psc[1]<par.cov_pc) { do_log(reasons,this_index,28); if (toFilt) continue; }
			
			// skip by missing rate discrepancy. Put it before missingrate so that I can see which sample (ct/cs) have more missing.
			bool Mis_applied=false;
			if (iiT)
			{
				MsRCS = iiA ? (double)miA/iiA : std::numeric_limits<double>::signaling_NaN();
				MsRCT = iiU ? (double)miU/iiU : std::numeric_limits<double>::signaling_NaN();
				if (CohLoc.empty() && perch::Mis_ea && perch::MisCut!=1)
				{
					if (!std::isnan(MsRCS)) { Mis_applied=true; if (MsRCS>perch::MisCut) { do_log(reasons,this_index,29); if (toFilt) continue; } }
					if (!std::isnan(MsRCT)) { Mis_applied=true; if (MsRCT>perch::MisCut) { do_log(reasons,this_index,29); if (toFilt) continue; } }
				}
				if (par.MR_pvl && iiA && iiU)
				{
					if (MsRCS || MsRCT)
					{
						int table[4];
						table[0]=miA;
						table[1]=iiA-miA;
						table[2]=miU;
						table[3]=iiU-miU;
						testMR = Fishers_exact_test_2x2(table,false);
					}
					else
					{
						testMR = 1;
					}
					if (testMR<=par.MR_pvl)
					{
						if (MsRCS>MsRCT)	{ do_log(reasons,this_index,30); if (toFilt) continue; }
						else				{ do_log(reasons,this_index,31); if (toFilt) continue; }
					}
				}
			}
			if (!CohLoc.empty() && perch::MisCut!=1)
			{
				bool not_pass=false;
				for (auto &coh:CohLoc)
				{
					if (coh.first.empty()) continue;
					int missing=0, nonmiss=0;
					for (auto &i:coh.second)
					{
						if	(geno[i]==-1)	++missing;
						else 				++nonmiss;
					}
					if (missing+nonmiss)
					{
						Mis_applied=true;
						double MsRate = (double)missing / (missing+nonmiss);
						if (MsRate>perch::MisCut) { not_pass=true; break; }
					}
				}
				if (not_pass) { do_log(reasons,this_index,29); if (toFilt) continue; }
			}
			
			// skip by missingness
			if (iiO) // some or all valid samples
			{
				MsAll = (double)miO/iiO;
				if (perch::MisCut!=1 && !Mis_applied)
				{
					// Test missing rate in overall samples instead of case/control separately.
					// Pro: 1) retain the variants that have missing in a small cohort (eg, 1 out 4 cases);
					//      2) it is more reasonable to test aff+unaff+unknown rather than aff/unaff.
					// Con: 1) not robust to very imbalanced design (eg, 4 cases 1000 controls).
					//         In this case, it will lead to a high % of variants that are missing in cases, then the cases will be removed due to sample-wise missing rate.
					//         Although there is a missing rate discrepancy test, it will not protect this situation because of the small sample size.
					if (MsAll>perch::MisCut)	{ do_log(reasons,this_index,29); if (toFilt) continue; }
				}
			}
			else // no valid samples, they should be treated as missing
			{
				if (perch::MisCut!=1) { do_log(reasons,this_index,29); if (toFilt) continue; }
			}

			// skip by variability
			if ( par.ftNoGt &&  niO==0 )			{ do_log(reasons,this_index,32); if (toFilt) continue; } // NoGenotypes (this should not happen if missingness filter was applied)
			if ( par.ftNoVr && (xiO==0||xiO==niO) ) { do_log(reasons,this_index,33); if (toFilt) continue; } // NoVariation (not robust to related individuals in pedigrees).
			
			// skip by heterozygous haploidy
			if (par.hh_prp && n_haploidy)
			{
				double proportion = (double)ploidy_err / n_haploidy;
				if (proportion>=par.hh_prp) { do_log(reasons,this_index,34); if (toFilt) continue; }
			}
			
			// skip by HWD
			if (par.HW_pvl && is_auto)
			{
				// HWE test in each population if population is defined
				if (!PopLoc.empty())
				{
					bool HWD=false;
					vector<double> freq;
					for (auto &pop:PopLoc)
					{
						if (pop.first.empty()) continue;
						int obs_hets=0, obs_hom1=0, obs_hom2=0;
						for (auto &i:pop.second)
						{
							if (par.hwesnf && !samples[i].fdr) continue;
							if (samples[i].aff!=1) continue;
							if		(geno[i]==0) ++obs_hom1;
							else if (geno[i]==1) ++obs_hets;
							else if (geno[i]==2) ++obs_hom2;
						}
						double pValue = my_SNPHWE(obs_hets, obs_hom1, obs_hom2, p.het_probs);
						if (pValue<=par.HW_pvl) HWD=true;
						int N = (obs_hets + obs_hom1 + obs_hom2) * 2;
						if (N) freq.push_back( (obs_hets + obs_hom2 * 2.0) / N );
					}
					if (HWD) { do_log(reasons,this_index,36); if (toFilt) continue; }
					if (freq.size()>1)
					{
						double pj = 0;
						for (auto &x:freq) pj+=x;
						pj/=freq.size();
						if (pj!=0 && pj!=1)
						{
							Rosenberg = ( -pj*log(pj) ) + ( -(1-pj)*log(1-pj) );
							double sum=0;
							for (auto &pij:freq)
							{
								if (pij)	sum += pij*log(pij);
								if (1-pij)	sum += (1-pij)*log(1-pij);
							}
							sum /= freq.size();
							Rosenberg += sum;
						}
					}
				}
				// HWE test in controls. Do it only when population is not defined. It's not a good idea because at this time PCA has not been done yet. So hweict defaults to false.
				else if (xiU && xiU!=niU && par.hweict)
				{
					testHW = my_SNPHWE(cgU[1],cgU[0],cgU[2],p.het_probs);
					if (testHW<=par.HW_pvl) { do_log(reasons,this_index,35); if (toFilt) continue; }
				}
			}
		}
		
		// pedigrees.
		int num_MenErr=0;
		int num_denovo=0;
		stringstream dn_ss;
		if (!pedigrees.empty())
		{
			for (auto &p:pedigrees)
			{
				string gtp;
				for (size_t i=0; i<p.id.size(); ++i)
				{
					int ploidy = genepi::expected_ploidy(chr_num,bp,p.sex[i]);
					char g0 = ( p.col[i] ? GTP::read(in[p.col[i]],gpar) : GTP::mss);
					char g1 = GTP::fix_ploidy(g0,ploidy,ploidy_err);
					gtp.push_back(g1);
				}
				int errors=0;
				int denovo=0;
				for (size_t i=2; i<p.id.size(); ++i)
				{
					if (GTP::denovo(gtp[i],gtp[0],gtp[1]))
					{
						++denovo;
						if (!dnv_fn.empty()) dn_ss<<this_index<<'\t'<<p.id[i]<<endl;
					}
					else
					{
						if	(is_auto)	errors += GTP::MenErr_auto(gtp[i],gtp[0],gtp[1]);
						else if (is_X)	errors += GTP::MenErr_chrX(gtp[i],gtp[0],gtp[1]);
						else if (is_Y)	errors += GTP::MenErr_chrY(gtp[i],gtp[0],gtp[1]);
						else if (is_M)	errors += GTP::MenErr_chrM(gtp[i],gtp[0],gtp[1]);
					}
				}
				if (errors) ++num_MenErr; // 2+ Mendelian errors in one nuclear pedigree is likely explained by a genotyping error in the parent
				if		(denovo >1) ++num_MenErr; // 2+ de novo in one nuclear pedigree is likely explained by a genotyping error in the parent
				else if (denovo==1) ++num_denovo; // 1  de novo mutation
			}
			if (par.MenErr && num_MenErr>=par.MenErr) { do_log(reasons,this_index,37); if (toFilt) continue; }
			if (par.denovo && num_denovo>=par.denovo) { do_log(reasons,this_index,38); if (toFilt) continue; }
			if (par.uno_dn && num_denovo && fdMAC>=par.uno_dn) { do_log(reasons,this_index,39); if (toFilt) continue; } // not considering samples that are not founders and whose parent was not sequenced!!!
		}

		// replications, put it after case-control so that all QC are done, this is useful for the estimation of concordance rate
		if (!RepSpl.empty() && dup_fn.empty())
		{
			bool matched=false;
			vector<string>	rp;			// matched line
			string			rep_tbx;	// return by tabix
			string chr = genepi::convert_chr_num(chr_num);
			string chr_spe_gtp = boost::algorithm::replace_all_copy(par.rep_in[1],"@",chr);
			if (!FileExists(chr_spe_gtp+".tbi")) exit_error("cannot find "+chr_spe_gtp+".tbi");
			try { rep_tbx = exec("tabix "+chr_spe_gtp+" "+chr+":"+itos(bp)+"-"+itos(bp),false); }
			catch (const std::exception& error) { exit_error("tabix "+chr_spe_gtp+" "+chr+":"+itos(bp)+"-"+itos(bp)+" failed"); }
			if (!rep_tbx.empty())
			{
				rep_tbx.pop_back(); // rep_tbx ends with \n
				vector<string> rep_row;
				boost::split(rep_row,rep_tbx,boost::is_any_of("\n"));
				for (auto &row:rep_row)
				{
					boost::split(rp,row,boost::is_any_of("\t"));
					if (rp[3]==ref && rp[4]==alt && rp[1]==in[FldPos[0]]) { matched=true; break; }
				}
			}
			
			if (matched)
			{
				int num_discrepancies=0;
				vector<char> GenoVec(2);
				vector<int> GQ_Vec(2,-1); // initialization to -1 is not necessary, but do it anyway
				vector<int> DP_Vec(2,-1); // initialization to -1 is not necessary, but do it anyway
				vector<double> GP_Vec(2,-1); // initialization to -1 is not necessary, but do it anyway
				for (auto &DG:RepSpl)
				{
					GenoVec[0] = GTP::read_vcf(rp[DG.FNs[0]],chr_num,bp,DG.SXs[0],true,GQ_Vec[0],DP_Vec[0],GP_Vec[0],gpar);
					GenoVec[1] = GTP::read_vcf(in[DG.FNs[1]],chr_num,bp,DG.SXs[1],true,GQ_Vec[1],DP_Vec[1],GP_Vec[1],gpar); // this is different from above
					test_dup(GenoVec,GQ_Vec,DP_Vec);
					// if ( GQ_Vec[1]>=0 && GQ_Vec[1]<GTP::GQ_cut ) GTP::set_missing(GenoVec[1]); // this is different from above
					// if ( DP_Vec[1]>=0 && DP_Vec[1]<GTP::DP_cut ) GTP::set_missing(GenoVec[1]); // this is different from above
					// below is different from above
					if (!par.spl_qc.empty() && DG.Spl[1]>=0)
					{
						spl_mt.lock();
						if (GTP::usable(GenoVec[0]))
						{
							++samples[DG.Spl[1]].rnm;
							if (GTP::usable(GenoVec[1]))
							{
								++samples[DG.Spl[1]].rbg;
								if (GenoVec[0]==GenoVec[1]) ++samples[DG.Spl[1]].rcc;
								else						++num_discrepancies;
							}
						}
						spl_mt.unlock();
					}
					else
					{
						if (GTP::usable(GenoVec[0]) && GTP::usable(GenoVec[1]) && GenoVec[0]!=GenoVec[1]) ++num_discrepancies;
					}
				}
				if (num_discrepancies>=par.disDup) { do_log(reasons,this_index,20); if (toFilt) continue; }
			}
		}
		
		// write VQSLOD of the variants that passed all QC except VQSLOD
		if (!vqs_fn.empty() && !std::isnan(VQSLOD))
		{
			if (is_snv)
			{
				string more_info=".";
				if (is_ti) more_info="coding_ti";
				if (is_tv) more_info="coding_tv";
				if (known)	{ vqs_mt.lock(); vqs_log << "known\tSNV\t" << VQSLOD << "\t" << more_info << endl; vqs_mt.unlock(); }
				else		{ vqs_mt.lock(); vqs_log << "new\tSNV\t" << VQSLOD << "\t" << more_info << endl;   vqs_mt.unlock(); }
			}
			else
			{
				string more_info=".";
				if (known)	{ vqs_mt.lock(); vqs_log << "known\tInDel\t" << VQSLOD << "\t" << more_info << endl; vqs_mt.unlock(); }
				else		{ vqs_mt.lock(); vqs_log << "new\tInDel\t" << VQSLOD << "\t" << more_info << endl;   vqs_mt.unlock(); }
			}
			if (is_snv) { if (VQSLOD<perch::VQSsnv)		{ do_log(reasons,this_index,15); if (toFilt) continue; } }
			else		{ if (VQSLOD<perch::VQSidl)		{ do_log(reasons,this_index,15); if (toFilt) continue; } }
			if (perch::hardft)
			{
				if (perch::FiltQD) { if (infoQD<perch::FiltQD)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
				if (perch::FiltMQ) { if (infoMQ<perch::FiltMQ)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
				if (perch::FiltFS) { if (infoFS>perch::FiltFS)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
				if (perch::FiltHS) { if (infoHS>perch::FiltHS)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
				if (perch::FiltMR) { if (infoMR<perch::FiltMR)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
				if (perch::FiltRP) { if (infoRP<perch::FiltRP)	{ do_log(reasons,this_index,16); if (toFilt) continue; } }
			}
		}
		
		// summary of variants that passed all QC
		if (num_denovo && !dnv_fn.empty())
		{
			dnv_mt.lock();
			dnv_log<<dn_ss.str();
			dnv_mt.unlock();
		}

		if (iiT)
		{
			if (!pqq_fn.empty() && niA && niU)
			{
				int table[4];
				table[0]=niU-xiU;
				table[1]=xiU;
				table[2]=niA-xiA;
				table[3]=xiA;
				double pv = std::numeric_limits<double>::signaling_NaN();
				test_2x2_array(table,pv);
				if (SplAF>=qq_maf && (qq_xcl.empty() || !qq_xbed.contain(chr_num,bp)))
				{
					pqq_mt.lock();
					pqq_pv.insert(pv);
					pqq_log	<<pv<<'\t'<<psc[0]*100<<'\t'<<psc[1]*100<<'\t'<<csctGQ<<'\t'<<csctDP<<'\t'<<testMR<<'\t'<<MsAll<<'\t'<<testHW<<'\t'<<testOA<<'\t'<<testOU<<'\t'
							<<in[FldChr[0]]+'\t'+in[FldPos[0]]+'\t'+'.'+'\t'+in[FldRef[0]]+'\t'+in[FldAlt[0]]<<'\t'<<MaxAF<<'\t'<<SplAF<<'\n';
					pqq_mt.unlock();
				}
			}
			if (!par.qc_gtp)
			{
				stringstream ss;
				ss	<< in[FldChr[0]]<<DLMTR<<in[FldPos[0]]<<DLMTR<<in[FldRef[0]]<<DLMTR<<in[FldAlt[0]]<<DLMTR<<VQSLOD
				<< DLMTR << iiA-miA << DLMTR << cgA[0] << DLMTR << cgA[1] << DLMTR << cgA[2]
				<< DLMTR << iiU-miU << DLMTR << cgU[0] << DLMTR << cgU[1] << DLMTR << cgU[2];
				in.clear();
				in.push_back(ss.str());
			}
		}

		double maf = (std::isnan(IntAF) ? 0 : std::min (IntAF, 1-IntAF));
		bool is_exome_rv = (maf<=par.sum_RV && min_bed.contain(chr_num,bp));
		if (is_auto && par.sum_RV && is_exome_rv)
		{
			nrv_mt.lock();
			if ( (xiF!=0 && xiF!=niF) || par.TotSpl )
			{
							  ++fnd_RV[0]; if (hom_singleton||het_singleton) { ++fnd_SV[0]; if (!known) ++fnd_PV[0]; }
				if (is_snv) { ++fnd_RS[0]; if (hom_singleton||het_singleton) { ++fnd_SS[0]; if (!known) ++fnd_PS[0]; } }
			}
			if (xiU!=0 && xiU!=niU)
			{
							  ++fnd_RV[1]; if (hom_singleton||het_singleton) { ++fnd_SV[1]; if (!known) ++fnd_PV[1]; }
				if (is_snv) { ++fnd_RS[1]; if (hom_singleton||het_singleton) { ++fnd_SS[1]; if (!known) ++fnd_PS[1]; } }
			}
			if (xiA!=0 && xiA!=niA)
			{
							  ++fnd_RV[2]; if (hom_singleton||het_singleton) { ++fnd_SV[2]; if (!known) ++fnd_PV[2]; }
				if (is_snv) { ++fnd_RS[2]; if (hom_singleton||het_singleton) { ++fnd_SS[2]; if (!known) ++fnd_PS[2]; } }
			}
			nrv_mt.unlock();
		}
		
		if (!par.spl_qc.empty())
		{
			spl_mt.lock();
			// for mean GQ and DP
			for (size_t i=0;i<samples.size();++i)
			{
				if (r_GQ[i]!=-1) { ++samples[i].nGQ; samples[i].tGQ+=r_GQ[i]; }
				if (r_DP[i]!=-1) { ++samples[i].nDP; samples[i].tDP+=r_DP[i]; }
				if (r_GP[i]!=-1) { ++samples[i].nGP; samples[i].tGP+=r_GP[i]; }
			}
			
			// Ti/Tv for haploid is different from diploid; het/hom assumes HWE; rhe/rho/oes/oos/nes/nos should be comparable between studies.
			if (is_auto)
			{
				for (size_t i=0;i<samples.size();++i)
				{
					if		(geno[i]==1) { ++samples[i].nms; ++samples[i].het; if (is_exome_rv)              ++samples[i].rhe; if (is_ti) ++samples[i].eti; if (is_tv) ++samples[i].etv; }
					else if (geno[i]==2) { ++samples[i].nms; ++samples[i].aho; if (is_exome_rv && IntAF<0.5) ++samples[i].rho; if (is_ti) ++samples[i].eti; if (is_tv) ++samples[i].etv; }
					else if (geno[i]==0) { ++samples[i].nms;                   if (is_exome_rv && IntAF>0.5) ++samples[i].rho; if (is_ti) ++samples[i].eti; if (is_tv) ++samples[i].etv; }
					else				 { ++samples[i].mss; }
				}
				if (is_exome_rv)
				{
					if (het_singleton) { ++samples[last_i].oes; if (!known) ++samples[last_i].nes; }
					if (hom_singleton) { ++samples[last_i].oos; if (!known) ++samples[last_i].nos; }
				}
			}
			else if (is_X)
			{
				if (genepi::expected_ploidy(chr_num,bp,1)==1)
				{
					// previously FdrAF. on 2018-03-06 changed to "(niF>=1000 ? FdrAF : MaxAF)" to be robust to small samples.
					// double UseMAF = (niF>=1000 ? FdrAF : MaxAF);
					
					// But 1000 is a bit too large, making some studies to use MaxAF, which is an overestimate of AF and lead to the overestimation of expected het.
					// FdrAF also have the same problem too. So I add a criteria fdMAC>=5, so that FdrAF is more reliable and I don't need MaxAF at all.
					// It's ok to ignore rare variants as they are not informative for sex inference. I took min(FdrMAF,MaxMAF) so that the var is real and common and AF not inflated.
					double FdrMAF = ( fdMAC>=5 ? std::min(FdrAF,1-FdrAF) : 0);
					double MaxMAF = ( MaxAF>0 ? std::min(MaxAF,1-MaxAF) : 0);
					double UseMAF = std::min(FdrMAF,MaxMAF);
					
					if (!std::isnan(UseMAF) && UseMAF!=0 && UseMAF!=1)
					{
						set<int>::iterator it1 = chrX_bp.lower_bound(bp-par.x_dist);
						set<int>::iterator it2 = chrX_bp.upper_bound(bp+par.x_dist);
						if (it2==it1) // not wthin +- 10k of any previous variants
						{
							chrX_bp.insert(bp);
							double pq2 = 2*UseMAF*(1-UseMAF);
							double p_het_di = pq2*(1-par.x_err1) + (1-pq2)*par.x_err1;
							double p_het_ha = par.x_err1;
							double p_hom_di = 1-p_het_di;
							double p_hom_ha = 1-p_het_ha;
							for (size_t i=0;i<samples.size();++i)
							{
								if		(geno[i]==1) { ++samples[i].nmx; ++samples[i].n1x; samples[i].lrx+=log10(p_het_di/p_het_ha); }
								else if (geno[i]==2) { ++samples[i].nmx; ++samples[i].n2x; samples[i].lrx+=log10(p_hom_di/p_hom_ha); }
								else if (geno[i]==0) { ++samples[i].nmx;				   samples[i].lrx+=log10(p_hom_di/p_hom_ha); }
								else				 { ++samples[i].msx; }
							}
						}
					}
					/*else // removed because these counts are useless. I already missed some due to LD pruning anyway.
					{
						for (size_t i=0;i<samples.size();++i)
						{
							if		(geno[i]==1) { ++samples[i].nmx; ++samples[i].n1x;  }
							else if (geno[i]==2) { ++samples[i].nmx; ++samples[i].n2x;  }
							else if (geno[i]==0) { ++samples[i].nmx;					}
							else				 { ++samples[i].msx;					}
						}
					}*/
				}
			}
			else if (is_Y)
			{
				if (genepi::expected_ploidy(chr_num,bp,1)==1)
				{
					for (size_t i=0;i<samples.size();++i)
					{
						if		(geno[i]==1) { ++samples[i].nmy; ++samples[i].n1y;  }
						else if (geno[i]==2) { ++samples[i].nmy; ++samples[i].n2y;  }
						else if (geno[i]==0) { ++samples[i].nmy;					}
						else				 { ++samples[i].msy;					}
					}
				}
			}
			spl_mt.unlock();
		}
		
		// Below are other filters not related to QC. The excluded variants should be in the summary statistics above (par.spl_qc), so put these parts here.

		// skip by allele frequency.
		if (perch::filXAF && !std::isnan(MaxAF))
		{
			double f = (MaxAF>0.5 ? 1-MaxAF : MaxAF);
			if (!perch::rf_XAF && f> perch::filXAF) { do_log(reasons,this_index,12); if (toFilt) continue; }
			if ( perch::rf_XAF && f<=perch::filXAF) { do_log(reasons,this_index,12); if (toFilt) continue; }
		}
		if (perch::filFAF && !std::isnan(FdrAF))
		{
			double f = (FdrAF>0.5 ? 1-FdrAF : FdrAF);
			if (!perch::rf_FAF && f> perch::filFAF) { do_log(reasons,this_index,12); if (toFilt) continue; }
			if ( perch::rf_FAF && f<=perch::filFAF) { do_log(reasons,this_index,12); if (toFilt) continue; }
		}
		if (perch::filSAF && !std::isnan(SplAF))
		{
			double f = (SplAF>0.5 ? 1-SplAF : SplAF);
			if (!perch::rf_SAF && f> perch::filSAF) { do_log(reasons,this_index,12); if (toFilt) continue; }
			if ( perch::rf_SAF && f<=perch::filSAF) { do_log(reasons,this_index,12); if (toFilt) continue; }
		}
		if (perch::filPAF && !std::isnan(PopAF))
		{
			double f = (PopAF>0.5 ? 1-PopAF : PopAF);
			if (!perch::rf_PAF && f> perch::filPAF) { do_log(reasons,this_index,12); if (toFilt) continue; }
			if ( perch::rf_PAF && f<=perch::filPAF) { do_log(reasons,this_index,12); if (toFilt) continue; }
		}

		// skip non-genic variants to speed up WGS analysis.
		if (par.ftAFNG && is_ii && ( MaxAF>0.5 ? 1-MaxAF<par.ftAFNG : MaxAF<par.ftAFNG))
		{
			if (!(BayesDel>=perch::BDELge)) { do_log(reasons,this_index,40); if (toFilt) continue; }
		}

		// select Ancestry Informative Markers (AIMs)
		if (par.filRos)
		{
			if (std::isnan(Rosenberg))	{ do_log(reasons,this_index,41); if (toFilt) continue; }
			if (Rosenberg<par.filRos)	{ do_log(reasons,this_index,41); if (toFilt) continue; }
		}

		// no more filters. write output
		flag[p.sqd_id][section+p.row_id]=1;
		if (reasons.empty()) do_log(reasons,this_index,0);
		if (par.qc_gtp && !FldInf.no_input() && !pr_qc)
		{
			bool exist_SumQC=false;
			bool exist_MisRt=false;
			bool exist_MisCs=false;
			bool exist_MisCt=false;
			bool exist_PopAF=false;
			bool exist_SplAF=false;
			bool exist_FdrAF=false;
			bool exist_FdrAC=false;
			bool exist_FdrAN=false;
			bool exist_FdrHO=false;
			bool exist_Cs_AF=false;
			bool exist_Ct_AF=false;
			bool exist_spMAC=false;
			bool exist_Rosen=false;
			for (auto &i:INFO)
			{
				if (str_startsw(i,"VICTOR_QC="))	{ i="VICTOR_QC="+reasons;		exist_SumQC=true; }
				if (str_startsw(i,"MissingRate="))	{ i="MissingRate="+ftos(MsAll); exist_MisRt=true; }
				if (str_startsw(i,"MissingCs="))	{ i="MissingCs="+ftos(MsRCS);   exist_MisCs=true; }
				if (str_startsw(i,"MissingCt="))	{ i="MissingCt="+ftos(MsRCT);   exist_MisCt=true; }
				if (str_startsw(i,"PopAF="))		{ i="PopAF="+ftos(PopAF);		exist_PopAF=true; }
				if (str_startsw(i,"SplAF="))		{ i="SplAF="+ftos(SplAF);		exist_SplAF=true; }
				if (str_startsw(i,"AF_Founder="))	{ i="AF_Founder="+ftos(FdrAF);	exist_FdrAF=true; }
				if (str_startsw(i,"AC_Founder="))	{ i="AC_Founder="+ftos(FdrAC);	exist_FdrAC=true; }
				if (str_startsw(i,"AN_Founder="))	{ i="AN_Founder="+itos(FdrAN);	exist_FdrAN=true; }
				if (str_startsw(i,"Hom_Founder="))	{ i="Hom_Founder="+itos(FdrHO);	exist_FdrHO=true; }
				if (str_startsw(i,"Cs_AF="))		{ i="Cs_AF="+ftos(Cs_AF);		exist_Cs_AF=true; }
				if (str_startsw(i,"Ct_AF="))		{ i="Ct_AF="+ftos(Ct_AF);		exist_Ct_AF=true; }
				if (str_startsw(i,"SampleMAC="))	{ i="SampleMAC="+itos(spMAC);	exist_spMAC=true; }
				if (str_startsw(i,"Rosenberg="))	{ i="Rosenberg="+ftos(Rosenberg); exist_Rosen=true; }
			}
			if (!exist_SumQC)							INFO.push_back("VICTOR_QC="+reasons);
			if (!exist_MisRt && !std::isnan(MsAll))		INFO.push_back("MissingRate="+ftos(MsAll));
			if (!exist_MisCs && !std::isnan(MsRCS))		INFO.push_back("MissingCs="+ftos(MsRCS));
			if (!exist_MisCt && !std::isnan(MsRCT))		INFO.push_back("MissingCt="+ftos(MsRCT));
			if (!exist_PopAF && !std::isnan(PopAF))		INFO.push_back("PopAF="+ftos(PopAF));
			if (!exist_SplAF && !std::isnan(SplAF))		INFO.push_back("SplAF="+ftos(SplAF));
			if (!exist_FdrAF && !std::isnan(FdrAF))		INFO.push_back("AF_Founder="+ftos(FdrAF));
			if (!exist_FdrAC 				      )		INFO.push_back("AC_Founder="+itos(FdrAC));
			if (!exist_FdrAN 				      )		INFO.push_back("AN_Founder="+itos(FdrAN));
			if (!exist_FdrHO 				      )		INFO.push_back("Hom_Founder="+itos(FdrHO));
			if (!exist_Cs_AF && !std::isnan(Cs_AF))		INFO.push_back("Cs_AF="+ftos(Cs_AF));
			if (!exist_Ct_AF && !std::isnan(Ct_AF))		INFO.push_back("Ct_AF="+ftos(Ct_AF));
			if (!exist_spMAC && spMAC!=-1)				INFO.push_back("SampleMAC="+itos(spMAC));
			if (!exist_Rosen && !std::isnan(Rosenberg)) INFO.push_back("Rosenberg="+ftos(Rosenberg));
			in[FldInf[0]]=str_of_container(INFO,';');
		}
	}
}

void worker_func(const int num_threads)
{
	for (int done=0, squad=0; done!=kNSq; iome[squad++].unlock())
	{
		if (squad==kNSq) squad=0;
		iome[squad].lock();
		if (rows[squad].empty()) { ++done; continue; }
		flag[squad].assign(rows[squad].size(),0);
		MtJobs<JobData> CompJobs;
		for (int i=0;i<num_threads;++i) CompJobs.push_back(JobData(squad,i,num_threads));
		CompJobs.run(compute, num_threads);
		for (size_t i=0;i<rows[squad].size();++i)
			if (flag[squad][i]) ofld.contents_to_ostream(rows[squad][i],program.outf,DLMTR,true); // print_container(rows[squad][i],program.outf,DLMTR,true);
		rows[squad].clear();
		flag[squad].clear();
		rsrc[squad]=0;
	}
}

int main (int argc, char * const argv[])
{
	// basic parameters
	perch::prelude(argc,argv);
	
	// other parameters
	string			vcf_in;			// genotype file in VCF format
	string			KINGprefix;		// prefix of KING outputs: $KINGprefix.kin and $KINGprefix.kin0. But only read $KINGprefix.kin0.
	double			KINGcutoff=0.1; // kinship cutoff, 0.40612 to remove duplicated samples only, 0.1 to remove related individuals
	double			PedErrPhKn=0.1;	// Phi&kinship cutoff, 0.1 by default
	tfile_format	format;
	format.set_delimiters("\t");
	format.set_option(SKIP_NOTES,false);
	bool			do_nothing=false;
	string			bed_in="minimum.bed";
	string			ex_chr = "1";	// example chromosome string (1)

	// handle program arguments
	perch::clear_flt_af(false);
	program.enable_option("--nt");
	program.enable_option("--prefix");
	program.read_arguments(argc,argv);
	perch::read_arguments();
	vector<string>	qc_file;
	for (size_t argi=1; argi<program.arg().size(); ++argi)
	{
		if (program.arg()[argi]=="--join-sample-qc") {
			while ((argi+1)<program.arg().size() && !program.arg()[argi+1].empty() && program.arg()[argi+1][0]!='-') qc_file.push_back(program.arg()[++argi]);
			if (qc_file.empty()) exit_error("no arguments for --join-sample-qc");
		}
		else if	(str_startsw(program.arg()[argi],"--spl"))				ReadArg(program.arg(),argi,par.spl_in);
		else if	(str_startsw(program.arg()[argi],"--cohort"))			ReadArg(program.arg(),argi,par.coh_in);
		else if	(str_startsw(program.arg()[argi],"--king"))				ReadArg(program.arg(),argi,KINGprefix);
		else if	(str_startsw(program.arg()[argi],"--kin-cutoff"))		ReadArg(program.arg(),argi,KINGcutoff);
		else if	(str_startsw(program.arg()[argi],"--phi-kin"))			ReadArg(program.arg(),argi,PedErrPhKn);
		else if (str_startsw(program.arg()[argi],"--ped"))				ReadArg(program.arg(),argi,par.ped_in);
		else if (str_startsw(program.arg()[argi],"--dup"))				ReadArg(program.arg(),argi,par.dup_in);
		else if (str_startsw(program.arg()[argi],"--rep"))				ReadSet(program.arg(),argi,par.rep_in);
		else if (str_startsw(program.arg()[argi],"--example-chr"))		ReadArg(program.arg(),argi,ex_chr);
		else if (str_startsw(program.arg()[argi],"--rm-var"))		{	ReadArg(program.arg(),argi,par.rm_var); if (!par.rm_var) par.ftAFNG=0; }
		else if (str_startsw(program.arg()[argi],"--do-nothing"))		ReadArg(program.arg(),argi,do_nothing);
		else if (str_startsw(program.arg()[argi],"--filt-no-geno"))		ReadArg(program.arg(),argi,par.ftNoGt);
		else if (str_startsw(program.arg()[argi],"--filt-no-var"))		ReadArg(program.arg(),argi,par.ftNoVr);
		else if (str_startsw(program.arg()[argi],"--filt-mac"))			ReadArg(program.arg(),argi,par.FltMAC);
		else if (str_startsw(program.arg()[argi],"--filt-impute2"))		ReadArg(program.arg(),argi,par.i2info);
		else if (str_startsw(program.arg()[argi],"--filt-miss-pv"))		ReadArg(program.arg(),argi,par.MR_pvl);
		else if (str_startsw(program.arg()[argi],"--rm-ii-MaxAF-lt"))	ReadArg(program.arg(),argi,par.ftAFNG);
		else if (str_startsw(program.arg()[argi],"--filt-Mendelian"))	ReadArg(program.arg(),argi,par.MenErr);
		else if (str_startsw(program.arg()[argi],"--filt-de-novo"))		ReadArg(program.arg(),argi,par.denovo);
		else if (str_startsw(program.arg()[argi],"--filt-uo-dn"))		ReadArg(program.arg(),argi,par.uno_dn);
		else if (str_startsw(program.arg()[argi],"--filt-hh"))			ReadArg(program.arg(),argi,par.hh_prp);
		else if (str_startsw(program.arg()[argi],"--filt-hwe-pv"))		ReadArg(program.arg(),argi,par.HW_pvl);
		else if (str_startsw(program.arg()[argi],"--filt-hwe-info"))	ReadSet(program.arg(),argi,par.InfoHW);
		else if (str_startsw(program.arg()[argi],"--filt-obs-pv"))		ReadArg(program.arg(),argi,par.AF_pvl);
		else if (str_startsw(program.arg()[argi],"--filt-obs-maf"))		ReadArg(program.arg(),argi,par.AF_cut);
		else if (str_startsw(program.arg()[argi],"--filt-af-diff"))		ReadArg(program.arg(),argi,par.AF_dif);
		else if (str_startsw(program.arg()[argi],"--filt-discord"))		ReadArg(program.arg(),argi,par.disDup);
		else if (str_startsw(program.arg()[argi],"--filt-rosenberg"))	ReadArg(program.arg(),argi,par.filRos);
		else if (str_startsw(program.arg()[argi],"--filt-an"))			ReadArg(program.arg(),argi,par.min_AN);
		else if (str_startsw(program.arg()[argi],"--filt-cov-DP"))		ReadArg(program.arg(),argi,par.cov_dp);
		else if (str_startsw(program.arg()[argi],"--filt-cov-pc"))		ReadArg(program.arg(),argi,par.cov_pc);
		else if (str_startsw(program.arg()[argi],"--filt-min-DP"))		ReadArg(program.arg(),argi,par.mnAvDP);
		else if (str_startsw(program.arg()[argi],"--filt-max-DP"))		ReadArg(program.arg(),argi,par.mxAvDP);
		else if (str_startsw(program.arg()[argi],"--hwe-founder"))		ReadArg(program.arg(),argi,par.hwesnf);
		else if (str_startsw(program.arg()[argi],"--hwe-controls"))		ReadArg(program.arg(),argi,par.hweict);
		else if (str_startsw(program.arg()[argi],"--no-filt-known"))	ReadArg(program.arg(),argi,par.noFtKn);
		else if (str_startsw(program.arg()[argi],"--no-filt-bad"))		ReadArg(program.arg(),argi,par.noFtAl);
		else if (str_startsw(program.arg()[argi],"--snv-only"))			ReadArg(program.arg(),argi,par.SNonly);
		else if (str_startsw(program.arg()[argi],"--indel-only"))		ReadArg(program.arg(),argi,par.IDonly);
		else if (str_startsw(program.arg()[argi],"--auto-chr"))			ReadArg(program.arg(),argi,par.AutoCh);
		else if (str_startsw(program.arg()[argi],"--skip-sv"))			ReadArg(program.arg(),argi,par.skipSV);
		else if (str_startsw(program.arg()[argi],"--qc-prv"))			ReadArg(program.arg(),argi,par.qc_prv);
		else if (str_startsw(program.arg()[argi],"--qc-gtp"))			ReadArg(program.arg(),argi,par.qc_gtp);
		else if (str_startsw(program.arg()[argi],"--qc-ploidy"))		ReadArg(program.arg(),argi,par.qc_pld);
		else if (str_startsw(program.arg()[argi],"--qc-mnp-only"))		ReadArg(program.arg(),argi,par.qc_mnp);
		else if (str_startsw(program.arg()[argi],"--mem"))				ReadArg(program.arg(),argi,par.memcap);
		else if (str_startsw(program.arg()[argi],"--tot-spl"))			ReadArg(program.arg(),argi,par.TotSpl);
		else if (str_startsw(program.arg()[argi],"--rv"))				ReadArg(program.arg(),argi,par.sum_RV);
		else if (str_startsw(program.arg()[argi],"--min-bed"))			ReadArg(program.arg(),argi,bed_in);
		else if (str_startsw(program.arg()[argi],"--out-vqs"))			ReadArg(program.arg(),argi,vqs_fn);
		else if (str_startsw(program.arg()[argi],"--qq"))				ReadArg(program.arg(),argi,pqq_fn);
		else if (str_startsw(program.arg()[argi],"--maf-for-qq"))		ReadArg(program.arg(),argi,qq_maf);
		else if (str_startsw(program.arg()[argi],"--excl-for-qq"))		ReadArg(program.arg(),argi,qq_xcl);
		else if (str_startsw(program.arg()[argi],"--log"))				ReadArg(program.arg(),argi,log_fn);
		else if (str_startsw(program.arg()[argi],"--out-dup"))			ReadArg(program.arg(),argi,dup_fn);
		else if (str_startsw(program.arg()[argi],"--out-dn"))			ReadArg(program.arg(),argi,dnv_fn);
		else if (str_startsw(program.arg()[argi],"--insert-dn"))		ReadArg(program.arg(),argi,dnv_in);
		else if (str_startsw(program.arg()[argi],"--sample-qc"))		ReadArg(program.arg(),argi,par.spl_qc);
		else if (str_startsw(program.arg()[argi],"--y-call-numb-m"))	ReadArg(program.arg(),argi,par.y_n_mn);
		else if (str_startsw(program.arg()[argi],"--y-call-rate-m"))	ReadArg(program.arg(),argi,par.y_r_mn);
		else if (str_startsw(program.arg()[argi],"--y-total-var-f"))	ReadArg(program.arg(),argi,par.y_n_wm);
		else if (str_startsw(program.arg()[argi],"--y-call-rate-f"))	ReadArg(program.arg(),argi,par.y_r_wm);
		else if (str_startsw(program.arg()[argi],"--x-err-rate"))		ReadArg(program.arg(),argi,par.x_err1);
		else if (str_startsw(program.arg()[argi],"--x-distance"))		ReadArg(program.arg(),argi,par.x_dist);
		else if (str_startsw(program.arg()[argi],"--info-ac"))			ReadArg(program.arg(),argi,par.InfoAC);
		else if (str_startsw(program.arg()[argi],"--info-an"))			ReadArg(program.arg(),argi,par.InfoAN);
		else if (str_startsw(program.arg()[argi],"--pv"))			{	ReadArg(program.arg(),argi,par.HW_pvl); par.MR_pvl=par.AF_pvl=par.HW_pvl; }
		else if (str_startsw(program.arg()[argi],"-")) exit_error("unknown option "+program.arg()[argi]);
		else if (vcf_in.empty()) vcf_in=program.arg()[argi];
		else { exit_error("excessive parameter "+program.arg()[argi]); }
	}
	
	// show help
	program.help_text_var("_Default_sample_file",par.spl_in);
	program.help_text_var("_Default_cohort",par.coh_in);
	program.help_text_var("_Default_rm_var",str_YesOrNo(par.rm_var));
	program.help_text_var("_Default_filt_mac",itos(par.FltMAC));
	program.help_text_var("_Default_filt_miss_d",ftos(par.MR_pvl));
	program.help_text_var("_Default_filt_i2info",ftos(par.i2info));
	program.help_text_var("_Default_filt_hwe_pv",ftos(par.HW_pvl));
	program.help_text_var("_Default_filt_hwe_info",str_of_container(par.InfoHW,',',false));
	program.help_text_var("_Default_filt_obs_pv",ftos(par.AF_pvl));
	program.help_text_var("_Default_filt_obs_cut",ftos(par.AF_cut));
	program.help_text_var("_Default_filt_af_diff",ftos(par.AF_dif));
	program.help_text_var("_Default_filt_discord",itos(par.disDup));
	program.help_text_var("_Default_filt_an",itos(par.min_AN));
	program.help_text_var("_Default_snv_only",str_YesOrNo(par.SNonly));
	program.help_text_var("_Default_indel_only",str_YesOrNo(par.IDonly));
	program.help_text_var("_Default_auto_chr",str_YesOrNo(par.AutoCh));
	program.help_text_var("_Default_skip_sv",str_YesOrNo(par.skipSV));
	program.help_text_var("_Default_qc_prv",str_YesOrNo(par.qc_prv));
	program.help_text_var("_Default_qc_gtp",str_YesOrNo(par.qc_gtp));
	program.help_text_var("_Default_qc_mnp",str_YesOrNo(par.qc_mnp));
	program.help_text_var("_Default_qc_pld",str_YesOrNo(par.qc_pld));
	program.help_text_var("_Default_noFtKn",str_YesOrNo(par.noFtKn));
	program.help_text_var("_Default_noFtBad",str_YesOrNo(par.noFtAl));
	program.help_text_var("_Default_do_nothing",str_YesOrNo(do_nothing));
	program.help_text_var("_Default_hwe_founder",str_YesOrNo(par.hwesnf));
	program.help_text_var("_Default_hwe_controls",str_YesOrNo(par.hweict));
	program.help_text_var("_Default_filt_no_var",str_YesOrNo(par.ftNoVr));
	program.help_text_var("_Default_filt_AF_NG",ftos(par.ftAFNG));
	program.help_text_var("_Default_spl_qc",par.spl_qc);
	program.help_text_var("_Default_memory",par.memcap);
	program.help_text_var("_Default_dup",par.dup_in);
	program.help_text_var("_Default_rep",str_of_container(par.rep_in,',',false));
	program.help_text_var("_Default_example_chr",ex_chr);
	program.help_text_var("_Default_out_vqs",vqs_fn);
	program.help_text_var("_Default_out_log",log_fn);
	program.help_text_var("_Default_out_dup",dup_fn);
	program.help_text_var("_Default_out_dn",dnv_fn);
	program.help_text_var("_Default_insert_dn",dnv_in);
	program.help_text_var("_Default_y_call_numb_m",itos(par.y_n_mn));
	program.help_text_var("_Default_y_call_rate_m",ftos(par.y_r_mn));
	program.help_text_var("_Default_y_total_var_f",itos(par.y_n_wm));
	program.help_text_var("_Default_y_call_rate_f",ftos(par.y_r_wm));
	program.help_text_var("_Default_x_err_rate",ftos(par.x_err1));
	program.help_text_var("_Default_x_distance",itos(par.x_dist));
	program.help_text_var("_Default_rv",ftos(par.sum_RV));
	program.help_text_var("_Default_min_bed",bed_in);
	program.help_text_var("_Default_tot_spl",itos(par.TotSpl));
	program.help_text_var("_Default_Mendelian",itos(par.MenErr));
	program.help_text_var("_Default_de_novo",itos(par.denovo));
	program.help_text_var("_Default_uno_dn",itos(par.uno_dn));
	program.help_text_var("_Default_ped",par.ped_in);
	program.help_text_var("_Default_hh",ftos(par.hh_prp));
	program.help_text_var("_Default_rosenberg",ftos(par.filRos));
	program.help_text_var("_Default_info_ac",par.InfoAC);
	program.help_text_var("_Default_info_an",par.InfoAN);
	program.help_text_var("_Default_cov_pc",ftos(par.cov_pc));
	program.help_text_var("_Default_cov_dp",itos(par.cov_dp));
	program.help_text_var("_Default_min_dp",itos(par.mnAvDP));
	program.help_text_var("_Default_max_dp",itos(par.mxAvDP));
	program.help_text_var("_Default_qq_file",pqq_fn);
	program.help_text_var("_Default_qq_maf",ftos(qq_maf));
	program.help_text_var("_Default_qq_xcl",qq_xcl);
	program.help_text_var("_Default_king",KINGprefix);
	program.help_text_var("_Default_kin_cutoff",ftos(KINGcutoff));
	program.help_text_var("_Default_phi_kin",ftos(PedErrPhKn));
	program.push_back_help(GTP::help_text());
	perch::check_arguments();

	if (do_nothing)
	{
		for (Lines_in_File(in, vcf_in, &format)) program.outf << in[0] << endl;
		return 0;
	}

	// check errors
	if (perch::preval==0)	exit_error("prevalence is not set");
	int_fast64_t membyt; // par.memcap in bytes
	if (qc_file.empty())
	{
		if (par.SNonly && par.IDonly) exit_error("--snv-only and --indel-only cannot be used at the same time.");
		if (par.sum_RV<0 || par.sum_RV>=0.5) exit_error("--rv must be within [0,0.5)");
		if (!par.rep_in.empty())
		{
			// if (!par.dup_in.empty()) exit_error("--dup and --rep cannot be used together.");
			if (par.rep_in.size()!=2) exit_error("--rep argumnts should be ReplicationSampleFile,ReplicationGenotypeFile");
			if (!str_endsw(par.rep_in[1],".vcf.gz")) exit_error("--rep argument 2 should be something.vcf.gz, which was compressed by bgzip.");
		}
		if (par.InfoHW.size() % 3) exit_error("Number of arguments for --hwe-info must be Nx3 strings, corresponding to AN,Het,Hom.");
		if (perch::FiltQD==0 && \
			perch::FiltMQ==0 && \
			perch::FiltFS==0 && \
			perch::FiltHS==0 && \
			perch::FiltMR==0 && \
			perch::FiltRP==0 && \
			perch::VQSsnv==-INFINITY && \
			perch::VQSidl==-INFINITY)
			lns<<showw<<"no VQSLOD / hard filters. Data may not be clean."<<flush_logger;
		
		if (!pqq_fn.empty())
		{
			if (!str_endsw(pqq_fn,".pdf") || pqq_fn==".pdf") exit_error("--qq filename must be something.pdf");
			pqq_fn.resize(pqq_fn.size()-4);
		}
		
		if		(str_endsw(par.memcap,"G")||str_endsw(par.memcap,"g")) { membyt = extract_double(par.memcap) / kNSq * 1073741824; }
		else if (str_endsw(par.memcap,"M")||str_endsw(par.memcap,"m")) { membyt = extract_double(par.memcap) / kNSq * 1048576; }
		else exit_error("--mem argument must ends with G or M.");
	}
	else
	{
		if (!dup_fn.empty()) exit_error("--out-dup and --sample-qc cannot be used together.");
	}

	// read sample file. columns: SeqID Sex Aff Cov(s).
	// Covariates can be case-insensitive strings or numbers, but cannot be both. Missing values are emptye strings or . or UNKNOWN.
	set< string >					h_csID;	// SeqID of cases
	set< string >					h_ctID;	// SeqID of controls
	set< string >					h_ukID;	// SeqID of unknowns
	map< string, int>				SexMap;	// SeqID => gender (1 for male, 2 for female, 0 for unknown)
	map< string, string> 			PopMap;	// SeqID => population
	map< string, vector<double> >	CovMap;	// SeqID => covariate for analysis
	if (!par.spl_in.empty())
	{
		map< string, double >		DepMap;	// SeqID => dependent variable
		map< string, string >		StrMap;	// SeqID => strata
		vector<string>				CovNID;	//          covariate new ID. [c] is the same as CovMap.
		perch::read_spl(par.spl_in,false,false,false,0.0,h_csID,h_ctID,h_ukID,SexMap,DepMap,PopMap,StrMap,CovMap,CovNID);
	}
	else
	{
		//if (par.TotSpl<=0) exit_error("without --spl, you must specify --tot-spl with a positive integer.");
	}
	
	map< string,string >			CohMap;	// SeqID => cohort
	if (!par.coh_in.empty())
	{
		const int Spl_ID=0;
		const int SplCoh=1;
		for (Rows_in_File(in,par.coh_in,2))
		{
			boost::to_upper(in[SplCoh]);
			if (in[SplCoh]!="" && in[SplCoh]!="." && in[SplCoh]!="UNKNOWN") CohMap[in[Spl_ID]]=in[SplCoh];
		}
	}

	// join_qc
	string sample_qc_header_str="SeqID\tSex\tms\tnon-ms\tHet\tAltHom\tHetRare\tHomRare\tHetPers\tHomPers\tHetSing\tHomSing\tTi\tTv\tmsRate\tHet/Hom\tHetRate\tTi/Tv\tY_ms\tY_call\tY_call%\tY_Sex\tX_call\tX_het\tX_LOD\tX_Sex\tRepGeno\tRepBthG\tRepConc\tCall%\tConc%\tnumGQ\tnumDP\tmeanGQ\tmeanDP";
	if (!qc_file.empty())
	{
		if (program.prefix().empty()) exit_error("--prefix cannot be empty with --join-sample-qc");
		if (!KINGprefix.empty() && par.spl_in.empty()) exit_error("with --join-sample-qc and --king, please also set --spl");
		vector<int> num_samples(2,0);
		vector<int> tot_rarevar(2,0);
		vector<string> sample_qc_header_row;
		boost::split(sample_qc_header_row,sample_qc_header_str,boost::is_any_of("\t"));
		tfile_format spl_qc_format;
		spl_qc_format.set_titlelines(1);
		spl_qc_format.comment_sw()="##";
		spl_qc_format.set_field_nums(sample_qc_header_row.size(),"line(s) do not have sufficent columns",tfile_format::Break);
		int file_num=0;
		for (auto &file:qc_file)
		{
			if (file_num)
				for (Rows_in_File(in,file,&spl_qc_format))
				{
					if (in.is_header())
					{
						if (in.contents()!=sample_qc_header_row) exit_error("The file "+file+" is not a vQC --sample-qc output file, or version wrong");
						continue;
					}
					Sample& spl=samples[in.RowNumber()-1];
					if (spl.ID != in[0]) exit_error("sample QC files did not match");
					double val;
					int nGQ=0, nDP=0;
					if (read_val(in[2] ,val)) spl.mss+=val; // missing
					if (read_val(in[3] ,val)) spl.nms+=val; // non-missing
					if (read_val(in[4] ,val)) spl.het+=val; // heterozygous
					if (read_val(in[5] ,val)) spl.aho+=val; // non-reference homozygous
					if (read_val(in[6] ,val)) spl.rhe+=val; // heterozygous RV
					if (read_val(in[7] ,val)) spl.rho+=val; // homozygous RV
					if (read_val(in[8] ,val)) spl.nes+=val; // New Het Singleton
					if (read_val(in[9] ,val)) spl.nos+=val; // New Hom Singleton
					if (read_val(in[10],val)) spl.oes+=val; // Obs Het Singleton
					if (read_val(in[11],val)) spl.oos+=val; // Obs Hom Singleton
					if (read_val(in[12],val)) spl.eti+=val; // (protein-coding) exon transition
					if (read_val(in[13],val)) spl.etv+=val; // (protein-coding) exon transversion
					if (read_val(in[18],val)) spl.msy+=val; //     missing Y genotypes
					if (read_val(in[19],val)) spl.nmy+=val; // non-missing Y genotypes
					if (read_val(in[22],val)) spl.nmx+=val; // non-missing X genotypes used for sex inference
					if (read_val(in[23],val)) spl.n1x+=val; // het X
					if (read_val(in[24],val)) spl.lrx+=val; // LOD score X
					if (read_val(in[26],val)) spl.rnm+=val; // replication non-missing
					if (read_val(in[27],val)) spl.rbg+=val; // replication both genotyped
					if (read_val(in[28],val)) spl.rcc+=val; // replication concordance
					if (read_val(in[31],nGQ)) spl.nGQ+=nGQ; // number of GQ
					if (read_val(in[32],nDP)) spl.nDP+=nDP; // number of DP
					if (read_val(in[33],val)) spl.tGQ+=val*nGQ; // sum of GQ
					if (read_val(in[34],val)) spl.tDP+=val*nDP; // sum of DP
				}
			else
				for (Rows_in_File(in,file,&spl_qc_format))
				{
					if (in.is_header())
					{
						if (in.contents()!=sample_qc_header_row) exit_error("The file "+file+" is not a vQC --sample-qc output file, or version wrong");
						continue;
					}
					int aff=0; if (exist_element(h_ctID,in[0])) aff=1; else if (exist_element(h_csID,in[0])) aff=2;
					samples.push_back(Sample(0,aff,in[0],in.RowNumber(),true));
					Sample& spl=samples.back();
					read_val(in[1] ,spl.sex);
					read_val(in[2] ,spl.mss); // missing
					read_val(in[3] ,spl.nms); // non-missing
					read_val(in[4] ,spl.het); // heterozygous
					read_val(in[5] ,spl.aho); // non-reference homozygous
					read_val(in[6] ,spl.rhe); // heterozygous RV
					read_val(in[7] ,spl.rho); // homozygous RV
					read_val(in[8] ,spl.nes); // New Het Singleton
					read_val(in[9] ,spl.nos); // New Hom Singleton
					read_val(in[10],spl.oes); // Obs Het Singleton
					read_val(in[11],spl.oos); // Obs Hom Singleton
					read_val(in[12],spl.eti); // (protein-coding) exon transition
					read_val(in[13],spl.etv); // (protein-coding) exon transversion
					read_val(in[18],spl.msy); //     missing Y genotypes
					read_val(in[19],spl.nmy); // non-missing Y genotypes
					read_val(in[22],spl.nmx); // non-missing X genotypes used for sex inference
					read_val(in[23],spl.n1x); // het X
					read_val(in[24],spl.lrx); // LOD score X
					read_val(in[26],spl.rnm); // replication non-missing
					read_val(in[27],spl.rbg); // replication both genotyped
					read_val(in[28],spl.rcc); // replication concordance
					read_val(in[31],spl.nGQ); // number of GQ
					read_val(in[32],spl.nDP); // number of DP
					read_val(in[33],spl.tGQ); spl.tGQ*=spl.nGQ; // sum of GQ
					read_val(in[34],spl.tDP); spl.tDP*=spl.nDP; // sum of DP
				}
			++file_num;
		}
		openOutFile_or_exit(spl_out,program.prefix()+".sample_qc");
		openOutFile_or_exit(sex_out,program.prefix()+".sex_problem");
		spl_out<<sample_qc_header_str<<endl;
		int totXhet=0, totXcall=0, mnY=0, mnYcall=0, wmY=0, wmYcall=0;
		for (auto &spl:samples)
		{
			int Y_tot_var = spl.msy+spl.nmy;
			double Y_call_rate = (Y_tot_var>0 ? (double)spl.nmy/Y_tot_var : std::numeric_limits<double>::signaling_NaN());
			int P_sex = spl.sex;															// provided sex
			int Y_sex = ( (spl.nmy>=par.y_n_mn && Y_call_rate>par.y_r_mn) ? 1 : ( (Y_tot_var>=par.y_n_wm && Y_call_rate<par.y_r_wm) ? 2 : 0) );	// chr Y
			int X_sex = (spl.lrx > par.xSxLOD ? 2 : (spl.lrx<-par.xSxLOD ? 1 : 0) );		// chr X
			int G_sex = ((X_sex && Y_sex && X_sex!=Y_sex) ? 0 : std::max(X_sex,Y_sex));		// genetic sex
			int F_sex = ((G_sex && P_sex && G_sex!=P_sex) ? 0 : std::max(G_sex,P_sex));		// final sex
			spl_out << spl.ID << DLMTR << spl.sex << DLMTR << spl.mss << DLMTR << spl.nms << DLMTR << spl.het << DLMTR << spl.aho
			<< DLMTR << spl.rhe << DLMTR << spl.rho << DLMTR << spl.nes << DLMTR << spl.nos << DLMTR << spl.oes << DLMTR << spl.oos
			<< DLMTR << spl.eti << DLMTR << spl.etv
			<< DLMTR << ((spl.mss+spl.nms)>0 ? (double)spl.mss/(spl.mss+spl.nms) : -1)
			<< DLMTR << (spl.aho>0 ? (double)spl.het/spl.aho : -1)
			<< DLMTR << (spl.nms>0 ? (double)spl.het/spl.nms : -1)
			<< DLMTR << (spl.etv>0 ? (double)spl.eti/spl.etv : -1)
			<< DLMTR << spl.msy << DLMTR << spl.nmy << DLMTR << (std::isnan(Y_call_rate) ? -1 : Y_call_rate) << DLMTR << Y_sex
			<< DLMTR << spl.nmx << DLMTR << spl.n1x << DLMTR << spl.lrx << DLMTR << X_sex
			<< DLMTR << spl.rnm << DLMTR << spl.rbg << DLMTR << spl.rcc << DLMTR << (spl.rnm ? 100.0*spl.rbg/spl.rnm : -1) << DLMTR << (spl.rbg ? 100.0*spl.rcc/spl.rbg : -1)
			<< DLMTR << spl.nGQ << DLMTR << spl.nDP << DLMTR << (spl.nGQ?spl.tGQ/spl.nGQ:0) << DLMTR << (spl.nDP?spl.tDP/spl.nDP:0)
			<< endl;
			if		(X_sex && Y_sex && X_sex!=Y_sex) sex_out<<spl.ID<<"\tX_sex_different_from_Y_sex"<<endl;
			else if (G_sex && P_sex && P_sex!=G_sex) sex_out<<spl.ID<<"\tGenetic_sex_different_from_the_provided_sex"<<endl;
			else if (F_sex==1) { mnY+=Y_tot_var; mnYcall+=spl.nmy; totXhet+=spl.n1x; totXcall+=spl.nmx; }
			else if (F_sex==2) { wmY+=Y_tot_var; wmYcall+=spl.nmy; }
		}
		if (totXcall) spl_out<<"##X-het-hap-rate: "<<(double)totXhet/totXcall<<endl;
		if (mnY) spl_out<<"##Ycall%-mn: "<<(double)mnYcall/mnY<<endl;
		if (wmY) spl_out<<"##Ycall%-wm: "<<(double)wmYcall/wmY<<endl;
		closefile(spl_out);
		closefile(sex_out);
		
		// test the difference in #RVs/sample between cases and controls
		if (CovMap.empty())
		{
			Values<double> rv_in_ct, rv_in_cs;
			for (auto &spl:samples)
			{
				if 		(spl.aff==1) rv_in_ct.push_back(spl.rhe+spl.rho);
				else if	(spl.aff==2) rv_in_cs.push_back(spl.rhe+spl.rho);
			}
			if (rv_in_ct.size()>0) lns<<showl<<"The number of RVs per sample in controls is "<<rv_in_ct.get(STAT::MEAN)<<flush_logger;
			if (rv_in_cs.size()>0) lns<<showl<<"The number of RVs per sample in cases    is "<<rv_in_cs.get(STAT::MEAN)<<flush_logger;
			if (rv_in_ct.size()>1 && rv_in_cs.size()>1)
			{
				double rv_test_pv=two_sample_t_same_sd(rv_in_ct, rv_in_cs);
				if (rv_test_pv<0.05) 	lns<<showw<<"The number of RVs per sample is significantly different between cases and controls by a raw test. p-value="<<rv_test_pv<<flush_logger;
				else 					lns<<showl<<"The number of RVs per sample is not significantly different between cases and controls by a raw test. p-value="<<rv_test_pv<<flush_logger;
			}
		}
		else
		{
			int t_cs=0;
			int t_ct=0;
			for (auto &spl:samples)
			{
				if (!exist_element(CovMap,spl.ID)) continue;
				if 		(spl.aff==1) ++t_ct;
				else if	(spl.aff==2) ++t_cs;
			}
			int				np=2+CovMap.begin()->second.size();	// number of parameters (1 constance, 1 aff status, # covariates)
			int				ni=t_cs+t_ct;						// number of individuals
			Eigen::VectorXd y(ni);								// y
			Eigen::MatrixXd X(ni,np);							// X
			int j=0;
			for (auto &spl:samples)
			{
				if (!exist_element(CovMap,spl.ID)) continue;
				if (spl.aff!=1 && spl.aff!=2) continue;
				y(j)  =spl.rhe+spl.rho;
				X(j,0)=1;
				X(j,1)=spl.aff-1;
				if (perch::_Debug) program.outf<<spl.ID<<DLMTR<<y(j)<<DLMTR<<X(j,1);
				for (size_t c=0;c<CovMap[spl.ID].size();++c) { X(j,c+2)=CovMap[spl.ID][c]; if (perch::_Debug) program.outf<<DLMTR<<X(j,c+2); }
				if (perch::_Debug) program.outf<<endl;
				++j;
			}
			double pv = pv_1st_linear(X,y,false);
			if (pv<0.05) 	lns<<showw<<"The number of RVs per sample is significantly different between cases and controls adjusted for covariates. p-value="<<pv<<flush_logger;
			else 			lns<<showl<<"The number of RVs per sample is not significantly different between cases and controls adjusted for covariates. p-value="<<pv<<flush_logger;
		}
		
		// sample-wise QC based on relatedness
		if (!KINGprefix.empty())
		{
			// read KINGprefix+".kin0" for cryptic correlations and select samples to remove
			if (FileExists(KINGprefix+".kin0"))
			{
				set<string> ped_mem_id;
				if (!par.ped_in.empty())
				{
					for (Rows_in_File(pi,par.ped_in,5))
					{
						// skip header
						if (exist_element(perch::h_pid,boost::to_lower_copy(pi[0]))) continue;
						// skip unwanted
						if (pi[1]=="0" || exist_element(perch::rm_ind,pi[1])) continue;
						ped_mem_id.insert(pi[1]);
					}
				}
				
				vector< set<string> > groups;
				for (Rows_in_File(in,KINGprefix+".kin0",8))
				{
					if (in.RowNumber()==0) continue;
					double Kinship;
					if (!read_val(in[7],Kinship)) exit_error("read "+KINGprefix+".kin0 failed");
					if (Kinship<=KINGcutoff) continue;
					string& id1=in[1];
					string& id2=in[3];
					if (!exist_element(h_csID,id1) && !exist_element(h_ctID,id1)) continue;
					if (!exist_element(h_csID,id2) && !exist_element(h_ctID,id2)) continue;
					bool fnd;
					fnd=false; for (auto &spl:samples) { if (spl.ID==id1) { fnd=true; break; } } if (!fnd) continue;
					fnd=false; for (auto &spl:samples) { if (spl.ID==id2) { fnd=true; break; } } if (!fnd) continue;
					vector< set<string> >::iterator it1 = groups.end();
					vector< set<string> >::iterator it2 = groups.end();
					for (each_element(groups,g)) if (exist_element(*g, id1)) { it1=g; break; }
					for (each_element(groups,g)) if (exist_element(*g, id2)) { it2=g; break; }
					if		(it1==groups.end() && it2!=groups.end()) { it2->insert(id1); }
					else if (it1!=groups.end() && it2==groups.end()) { it1->insert(id2); }
					else if (it1!=groups.end() && it2!=groups.end()) { if (it1!=it2) { for (each_element(*it2,i)) it1->insert(*i); groups.erase(it2); } }
					else { groups.push_back(set<string>()); groups.back().insert(id1); groups.back().insert(id2); }
				}
				
				openOutFile_or_exit(rel_out,program.prefix()+".rel_problem");
				for (auto &g:groups)
				{
					string max_str;
					double max_val=-1;
					for (auto &i:g)
					{
						for (auto &spl:samples)
						{
							if (spl.ID==i)
							{
								double call_rate = ((spl.mss+spl.nms)>0 ? (double)spl.nms/(spl.mss+spl.nms) : -1);
								double in_ped = exist_element(ped_mem_id,spl.ID);
								double score = in_ped + call_rate;
								if (score>max_val) { max_val=score; max_str=spl.ID; }
								break;
							}
						}
					}
					for (auto &i:g)
					{
						if (i!=max_str)
							rel_out<<i<<"\trelated_to_"<<max_str<<endl;
					}
				}
				closefile(rel_out);
			}
			
			// read KINGprefix+".kin" for pedigree structure errors
			if (!par.ped_in.empty() && FileExists(KINGprefix+".kin"))
			{
				// read .kin
				set<string> 		PedErr_ind;
				map<string,double>	PedErr_phi;
				map<string,double>	PedErr_kin;
				for (Rows_in_File(in,KINGprefix+".kin",10))
				{
					if (in[0]=="FID") continue;
					if (in[9]!="1") continue;
					PedErr_ind.insert(in[0]+"\t"+in[1]);
					PedErr_ind.insert(in[0]+"\t"+in[2]);
					double phi; if (!read_val_noNaN(in[5],phi)) exit_error("failed to read Phi from "+KINGprefix+".kin");
					double kin; if (!read_val_noNaN(in[8],kin)) exit_error("failed to read Kinship from "+KINGprefix+".kin");
					if (phi<PedErrPhKn && kin<PedErrPhKn) continue;
					PedErr_phi.insert(std::make_pair(in[0]+"\t"+in[1]+"\t"+in[2],phi));
					PedErr_kin.insert(std::make_pair(in[0]+"\t"+in[1]+"\t"+in[2],kin));
				}
				
				if (!PedErr_ind.empty())
				{
					// check Phi
					vector<string> true_errors;
					if (!par.ped_in.empty())
					{
						string uID=random_string(12);
						openOutFile_or_exit(file1,perch::TMPDIR+"pedpro_"+uID+".ind");
						file1<<"PID\tIID\n";
						print_container(PedErr_ind,file1,'\n',true);
						closefile(file1);
						string result=exec("pedpro --ped \""+par.ped_in+"\" --kin-of \""+perch::TMPDIR+"pedpro_"+uID+".ind\" > \""+perch::TMPDIR+"pedpro_"+uID+".out\"",false);
						map<string,double>	true_kinship; // true kinship coefficient calculated by pedpro
						for (Rows_in_File(in,perch::TMPDIR+"pedpro_"+uID+".out",4))
						{
							if (in[3]=="kinship") continue;
							double kin; if (!read_val_noNaN(in[3],kin)) exit_error("failed to read Kinship from "+perch::TMPDIR+"pedpro_"+uID+".out");
							true_kinship.insert(std::make_pair(in[0]+"\t"+in[1]+"\t"+in[2],kin));
							true_kinship.insert(std::make_pair(in[0]+"\t"+in[2]+"\t"+in[1],kin));
						}
						for (auto &l:PedErr_phi)
						{
							if (exist_element(true_kinship,l.first)&&true_kinship[l.first]==l.second)
								true_errors.push_back(l.first+"\t"+ftos(l.second)+"\t"+ftos(PedErr_kin[l.first]));
						}
					}
					else
					{
						for (auto &l:PedErr_phi)
						{
							true_errors.push_back(l.first+"\t"+ftos(l.second)+"\t"+ftos(PedErr_kin[l.first]));
						}
					}
					
					// write output
					openOutFile_or_exit(PedErr_out,program.prefix()+".ped_errors");
					if (!true_errors.empty()) print_container(true_errors,PedErr_out,'\n',true);
					closefile(PedErr_out);
				}
			}
		}
		return 0;
	}
	
	// read qq_xcl
	if (!qq_xcl.empty()) qq_xbed.setup(qq_xcl,true,0,true);
	
	// read min_bed
	min_bed.setup(perch::find_file(bed_in),true,0,true);
	
	// read VCF and do QC
	if (!log_fn.empty()) { if (!openOutFile(logout, log_fn)) exit_error("cannot open log file."); }
	if (!dup_fn.empty()) { if (!openOutFile(GQ_DP_log, dup_fn)) exit_error("cannot open out file to record GQ,DP."); }
	int num_threads = std::max( program.nt-2, 1); // reduce one for main, one for worker_func
	int	this_squad=0;
	iome[this_squad].lock();
	format.set_storage_to(rows[this_squad]);
	std::thread worker(worker_func,num_threads);
	bool header_not_read = true;
	bool VCFHeader_not_w = true;
	int num_cs=0, num_ct=0, num_uk=0, num_fd=0;
	for (Rows_in_File(in, vcf_in, &format))
	{
		// read header
		if (in[0].size()>1 && in[0][0]=='#' && in[0][1]=='#')
		{
			in.clear_nf();
			perch::read_meta(in[0]);
			if ( par.h_dels.empty() && str_startsw(in[0],"##INFO=<ID=BayesDel")) ColDel=-1;
			if (!par.h_dels.empty() && str_startsw(in[0],"##INFO=<ID="+par.h_dels+",")) ColDel=-1;
			if (str_has(in[0],"##INFO=<ID=VICTOR_QC,")) VCFHeader_not_w=false;
			if (str_startsw(in[0],"##INFO=<ID=CADD,")) CADD_column=-1;
			if (str_startsw(in[0],"##INFO=<ID=fthmMKL_NC,")) FMNC_column=-1;
			if (str_startsw(in[0],"##INFO=<ID="+perch::h_SEGb+",")) ColSEG=-1;
			if (header_not_read && par.qc_gtp) print_container(in.contents(),program.outf,' ',true);
			rows[this_squad].clear();
			continue;
		}
		if (exist_any(perch::h_col1, in.contents()))
		{
			if (VCFHeader_not_w)
			{
				program.outf<<"##INFO=<ID=VICTOR_QC,Number=1,Type=String,Description=\"Result of quality control by VICTOR.\">"<<endl;
				program.outf<<"##INFO=<ID=MissingRate,Number=1,Type=Float,Description=\"Genotype missing rate in sequenced samples.\">"<<endl;
				program.outf<<"##INFO=<ID=MissingCs,Number=1,Type=Float,Description=\"Genotype missing rate in sequenced cases.\">"<<endl;
				program.outf<<"##INFO=<ID=MissingCt,Number=1,Type=Float,Description=\"Genotype missing rate in sequenced controls.\">"<<endl;
				program.outf<<"##INFO=<ID=PopAF,Number=1,Type=Float,Description=\"Population allele frequency calcalted from data.\">"<<endl;
				program.outf<<"##INFO=<ID=SplAF,Number=1,Type=Float,Description=\"Allele frequency in independent cases and controls.\">"<<endl;
				program.outf<<"##INFO=<ID=AF_Founder,Number=1,Type=Float,Description=\"Allele frequency in founders.\">"<<endl;
				program.outf<<"##INFO=<ID=AC_Founder,Number=1,Type=Integer,Description=\"ALT allele count in founders.\">"<<endl;
				program.outf<<"##INFO=<ID=AN_Founder,Number=1,Type=Integer,Description=\"Total allele count in founders.\">"<<endl;
				program.outf<<"##INFO=<ID=Hom_Founder,Number=1,Type=Integer,Description=\"Homozygous ALT genotypes in founders.\">"<<endl;
				program.outf<<"##INFO=<ID=Cs_AF,Number=1,Type=Float,Description=\"Allele frequency in cases.\">"<<endl;
				program.outf<<"##INFO=<ID=Ct_AF,Number=1,Type=Float,Description=\"Allele frequency in controls.\">"<<endl;
				program.outf<<"##INFO=<ID=SampleMAC,Number=1,Type=Integer,Description=\"Minor Allele Count in the pooled case-control samples.\">"<<endl;
				program.outf<<"##INFO=<ID=Rosenberg,Number=1,Type=Float,Description=\"Rosenberg's informativeness of genetic markers for inference of ancestry.\">"<<endl;
				VCFHeader_not_w=false;
			}
			if (!header_not_read) { continue; }
			header_not_read=false;
			format.clear_field_nums();
			FldChr.clear();
			FldPos.clear();
			FldSNP.clear();
			FldRef.clear();
			FldAlt.clear();
			FldFlt.clear();
			FldInf.clear();
			FldFmt.clear();
			FldXAF.clear();
			
			// read dup_in, write DupSpl
			if (!par.dup_in.empty())
			{
				tfile_format FmtDup;
				FmtDup.forbid_nf_rpt();
				size_t MaxCol=0;
				size_t TotCmp=0;
				for (Rows_in_File(dup_in,par.dup_in,&FmtDup))
				{
					// read
					DupGrp DG;
					size_t	NumCmp=0;
					bool	GoldStd=false;
					for (size_t j=0;j<dup_in.contents().size();++j)
					{
						int fn=-1;
						string& id=dup_in[j];
						for (int i=0;i<in.NumFields();++i) if (in[i]==id && !exist_element(perch::rm_ind,in[i])) { fn=i; if (j) ++NumCmp; else GoldStd=true; }
						if (fn>=0)
						{
							DG.IDs.push_back(id);
							DG.FNs.push_back(fn);
							DG.SXs.push_back(SexMap[id]);
							DG.Spl.push_back(-1);
						}
						else
						{
							DG.IDs.push_back("");
							DG.FNs.push_back(-1);
							DG.SXs.push_back(0);
							DG.Spl.push_back(-1);
						}
					}
					if (!GoldStd || NumCmp==0) continue;
					
					// check sex
					set<int> sexes;
					for (auto &s:DG.SXs) if (s) sexes.insert(s);
					if (sexes.size()>1) exit_error("Sex conflict for "+str_of_container(DG.IDs,','));
					// if (sexes.size()<1) { lns<<showw<<"sex unknown for "<<str_of_container(DG.IDs,',')<<". Skipped."<<flush_logger; continue; }
					if (sexes.size()==1) for (auto &s:DG.SXs) s=*sexes.begin();
					
					// keep only one in cs/ct/uk, prefer cs/ct
					int in_cs=0, in_ct=0, in_uk=0;
					vector<string> id_in_uk;
					for (size_t i=0;i<DG.IDs.size();++i)
					{
						if (exist_element(h_csID,DG.IDs[i])) { ++in_cs; }
						if (exist_element(h_ctID,DG.IDs[i])) { ++in_ct; }
						if (exist_element(h_ukID,DG.IDs[i])) { ++in_uk; id_in_uk.push_back(DG.IDs[i]); }
					}
					if		(in_cs+in_ct>1)	exit_error("more than one duplicated sample ("+str_of_container(DG.IDs,',')+") in independent cases and/or controls");
					else if (in_cs+in_ct<1)	for (size_t i=1;i<id_in_uk.size();++i) h_ukID.erase(id_in_uk[i]);
					else					for (size_t i=0;i<id_in_uk.size();++i) h_ukID.erase(id_in_uk[i]);
					
					// add
					TotCmp += NumCmp;
					DupSpl.push_back(DG);
					if (DG.FNs.size()>MaxCol) MaxCol=DG.FNs.size();
				}
				if (MaxCol)
				{
					concord_n.assign(MaxCol,0);
					concord_p.assign(MaxCol,0);
					false_pos.assign(MaxCol,0);
					false_neg.assign(MaxCol,0);
					SplSet1na.assign(MaxCol,0);
					soft_miss.assign(MaxCol,0);
					hard_miss.assign(MaxCol,0);
					both12_na.assign(MaxCol,0);
					het_2_hom.assign(MaxCol,0);
					hom_2_het.assign(MaxCol,0);
					lns<<showl<<"There are "<<DupSpl.size()<<" groups of duplications, "<<TotCmp<<" comparisons."<<flush_logger;
				}
				else
					exit_error("No duplicated samples read from "+par.dup_in);
			}

			// read cases/controls and other header
			for (int i=0;i<in.NumFields();++i)
			{
				if (in[i]==perch::h_SEGb && ColSEG==-2)		ColSEG=i;
				if (in[i]=="CADD" && CADD_column==-2) 		CADD_column=i;
				if (in[i]=="fthmMKL_NC" && FMNC_column==-2) FMNC_column=i;
				if (in[i]=="FILTER")				FldFlt.push_back(i+1);
				if (in[i]=="INFO")					FldInf.push_back(i+1);
				if (in[i]=="FORMAT")				FldFmt.push_back(i+1);
				if (in[i]==perch::h_MxAF)			FldXAF.push_back(i+1);
				if (in[i]=="ID")					FldSNP.push_back(i+1);
				if		(exist_element(h_csID,in[i])) {	samples.push_back(Sample(SexMap[in[i]],2,in[i],i,true)); ++num_cs; }
				else if (exist_element(h_ctID,in[i])) {	samples.push_back(Sample(SexMap[in[i]],1,in[i],i,true)); ++num_ct; }
				else if (exist_element(h_ukID,in[i])) {	samples.push_back(Sample(SexMap[in[i]],0,in[i],i,true)); ++num_uk; }
				if (exist_element(PopMap,in[i])) PopLoc[PopMap[in[i]]].push_back(samples.size()-1);
				if (exist_element(CohMap,in[i])) CohLoc[CohMap[in[i]]].push_back(samples.size()-1);
				if ( par.h_dels.empty() && str_startsw(in[i],"BayesDel"))	{ if (ColDel==-2) ColDel=i; else exit_error("multiple columns for BayesDel"); }
				if (!par.h_dels.empty() && in[i]==par.h_dels)				{ if (ColDel==-2) ColDel=i; else exit_error("multiple columns for "+par.h_dels); }
				if (in[i]=="Chr"   && FldChr.no_input())	FldChr.push_back(i+1); if (in[i]=="#CHROM")		{	FldChr.clear(); FldChr.push_back(i+1); }
				if (in[i]=="Start" && FldPos.no_input())	FldPos.push_back(i+1); if (in[i]=="POS")		{	FldPos.clear(); FldPos.push_back(i+1); }
				if (in[i]=="Ref"   && FldRef.no_input())	FldRef.push_back(i+1); if (in[i]=="REF")		{	FldRef.clear(); FldRef.push_back(i+1); }
				if (in[i]=="Alt"   && FldAlt.no_input())	FldAlt.push_back(i+1); if (in[i]=="ALT")		{	FldAlt.clear(); FldAlt.push_back(i+1); }
			}

			// test the validity of FltMAC.
			if (par.FltMAC)
			{
				if (par.TotSpl)
				{
					if (par.FltMAC>=par.TotSpl)
					{
						lns<<showw<<"sample size is too small; --filt-mac="<<par.FltMAC<<" is ignored. If you still want to use --filt-mac, change it to a smaller threshold."<<flush_logger;
						par.FltMAC=0;
					}
				}
				else
				{
					if (par.FltMAC>=(int)samples.size())
					{
						lns<<showw<<"sample size is too small; --filt-mac="<<par.FltMAC<<" is ignored. If you still want to use --filt-mac, change it to a smaller threshold."<<flush_logger;
						par.FltMAC=0;
					}
				}
			}

			// read pedigrees
			if (!par.ped_in.empty())
			{
				int num_trios=0;
				for (Rows_in_File(pi,par.ped_in,5))
				{
					// skip header
					if (exist_element(perch::h_pid,boost::to_lower_copy(pi[0]))) continue;
					
					// see whether sequenced
					int col_ch=0, col_pa=0, col_ma=0; // 0-based, 0=NotSequenced
					for (int i=1;i<in.NumFields();++i)
					{
						if (in[i]==pi[1]) col_ch=i;
						if (in[i]==pi[2]) col_pa=i;
						if (in[i]==pi[3]) col_ma=i;
					}

					// skip unwanted
					if (pi[1]=="0" || exist_element(perch::rm_ind,pi[1]) || !col_ch) continue;
					
					// read sex
					int	s = perch::read_sex(pi[4]);
					
					// add to Samples with aff=unk if not there. This is important for --filt-no-geno, --filt-no-var, --filt-missing-rate, --filt-uno-dn, etc.
					bool is_founder = (pi[2]=="0"&&pi[3]=="0");
					bool is_new = true;
					for (auto &spl:samples) if (spl.fld==col_ch) { is_new=false; spl.fdr=is_founder; break; }
					if (is_new) { samples.push_back(Sample(s,0,pi[1],col_ch,is_founder)); ++num_uk; }
					
					// skip unwanted
					if (pi[2]=="0" || exist_element(perch::rm_ind,pi[2])) continue;
					if (pi[3]=="0" || exist_element(perch::rm_ind,pi[3])) continue;
					if (!col_pa && !col_ma) continue; // skip both parent unsequenced
					if (col_pa && col_ma) ++num_trios;
					
					// add
					bool new_ped=true;
					for (auto &p:pedigrees)
					{
						if (pi[2]==p.id[0] && pi[3]==p.id[1])
						{
							p.id.push_back(pi[1]);
							p.sex.push_back(s);
							p.col.push_back(col_ch);
							new_ped=false;
							break;
						}
					}
					if (new_ped)
					{
						nuclear_ped p;
						p.id.push_back(pi[2]); p.sex.push_back(1); p.col.push_back(col_pa);
						p.id.push_back(pi[3]); p.sex.push_back(2); p.col.push_back(col_ma);
						p.id.push_back(pi[1]); p.sex.push_back(s); p.col.push_back(col_ch);
						pedigrees.push_back(p);
					}
				}
				lns<<showl<<"There're "<<pedigrees.size()<<" nuclear pedigrees with 1+ parent and 1+ offspring, among which there're "<<num_trios<<" sequenced full trios."<<flush_logger;
			}
			
			for (auto &spl:samples) if (spl.fdr) ++num_fd;
			lns<<showl<<"There're "<<num_cs<<" independent cases, "<<num_ct<<" independent controls, and "<<num_uk<<" others. Total number of founders = "<<num_fd<<flush_logger;
			
			if (!par.rep_in.empty())
			{
				vector<string> id;
				tfile_format rep_geno_format;
				rep_geno_format.set_delimiters("\t");
				rep_geno_format.comment_sw()="##";
				string chr_spe_gtp = boost::algorithm::replace_all_copy(par.rep_in[1],"@",ex_chr);
				if (!FileExists(chr_spe_gtp+".tbi")) exit_error("cannot find "+chr_spe_gtp+".tbi");
				for (Rows_in_File(ri,chr_spe_gtp,&rep_geno_format))
				{
					if (exist_any(perch::h_col1, ri.contents()))
					{
						id=ri.contents();
						break;
					}
					exit_error("faile to read "+chr_spe_gtp+". It should be in VCF format.");
				}
				
				for (Rows_in_File(ri,par.rep_in[0],2))
				{
					DupGrp DG;
					if (exist_element(perch::rm_ind,ri[1])) continue;
					for (size_t i=0;i<id.size();     ++i) if (id[i]==ri[0]) { DG.IDs.push_back(ri[0]); DG.FNs.push_back(i); DG.SXs.push_back(SexMap[ri[0]]); DG.Spl.push_back(-1); break; }
					for (int    i=0;i<in.NumFields();++i) if (in[i]==ri[1]) { DG.IDs.push_back(ri[1]); DG.FNs.push_back(i); DG.SXs.push_back(SexMap[ri[1]]); DG.Spl.push_back(-1); break; }
					if (DG.FNs.size()!=2) continue;
					
					// check sex
					set<int> sexes;
					for (auto &s:DG.SXs) if (s) sexes.insert(s);
					if (sexes.size()>1) exit_error("Sex conflict for "+str_of_container(DG.IDs,','));
					// if (sexes.size()<1) { lns<<showw<<"sex unknown for "<<str_of_container(DG.IDs,',')<<". Skipped."<<flush_logger; continue; }
					if (sexes.size()==1) for (auto &s:DG.SXs) s=*sexes.begin();
					
					// find .spl
					for (size_t i=0; i<DG.FNs.size(); ++i)
						for (size_t j=0; j<samples.size(); ++j)
							if (DG.FNs[i]==samples[j].fld) { DG.Spl[i]=j; break; }
					
					// add
					RepSpl.push_back(DG);
				}
				lns<<showl<<"There're "<<RepSpl.size()<<" samples replicated by another experiment."<<flush_logger;
				concord_n.assign(2,0);
				concord_p.assign(2,0);
				false_pos.assign(2,0);
				false_neg.assign(2,0);
				SplSet1na.assign(2,0);
				soft_miss.assign(2,0);
				hard_miss.assign(2,0);
				both12_na.assign(2,0);
				het_2_hom.assign(2,0);
				hom_2_het.assign(2,0);
			}
			
			if (!dnv_in.empty())
			{
				for (Rows_in_File(dni,dnv_in,5))
				{
					string index = dni[0]+'\t'+dni[1]+'\t'+dni[2]+'\t'+dni[3];
					int col=-1;
					for (int i=0;i<in.NumFields();++i) if (in[i]==dni[4]) { col=i; break; }
					if (col>=0) dnv_db[index].push_back(col);
				}
			}
			
			if (par.filRos)
			{
				if (par.HW_pvl==0) exit_error("For --filt-rosenberg to work, --filt-hwe-pv must not be 0");
				if (PopLoc.empty()) exit_error("For --filt-rosenberg to work, the Sample File must contain population information at the 4th column, and the samples should have an unknown affection status.");
				if (PopLoc.size()<2) exit_error("For --filt-rosenberg to work, there must be more than one populations.");
			}
			if (num_cs==0 && num_ct==0 && num_uk==0)
			{
				if (par.ftNoGt) { par.ftNoGt=false; lns<<showl<<"No case/control/pedigree samples, turnning off --filt-no-geno"<<flush_logger; }
				if (par.ftNoVr) { par.ftNoVr=false; lns<<showl<<"No case/control/pedigree samples, turnning off --filt-no-var"<<flush_logger; }
			}
			if (FldChr.no_input()) exit_error("The #CHROM/Chr column is missing.");
			if (FldPos.no_input()) exit_error("The POS/Start column is missing.");
			// if (FldSNP.no_input()) exit_error("The ID column is missing.");
			if (FldRef.no_input()) exit_error("The REF/Ref column is missing.");
			if (FldAlt.no_input()) exit_error("The ALT/Alt column is missing.");
			// if (FldInf.no_input()) exit_error("The INFO column is missing.");

			for (int i=0;i<in.NumFields();++i)
				if (!exist_element(perch::rm_ind,in[i])) ofld.push_back(i+1);
			
			if (par.qc_gtp)
				in.write_r(program.outf,ofld,true); // print_container(in.contents(),program.outf,DLMTR,true);
			else
				program.outf << "## There are " << num_cs << " cases and " << num_ct << " controls." <<endl
				<<"#CHROM"<<DLMTR<<"POS"<<DLMTR<<"ID"<<DLMTR<<"REF"<<DLMTR<<"ALT"<<DLMTR<<"VQSLOD"<<DLMTR
				<<"A:nm"<<DLMTR<<"A:0/0"<<DLMTR<<"A:0/1"<<DLMTR<<"A:1/1"<<DLMTR
				<<"U:nm"<<DLMTR<<"U:0/0"<<DLMTR<<"U:0/1"<<DLMTR<<"U:1/1"<<endl;
			rows[this_squad].clear();
			continue;
		}
		if (header_not_read) exit_error("Header lines missing.");
		
		for (auto& x:in.contents()) rsrc[this_squad]+=x.size();
		if (rsrc[this_squad]>=membyt && rows[this_squad].size()%num_threads==0)
		{
			int next_squad = ( this_squad+1==kNSq? 0 : this_squad+1 );
			iome[next_squad].lock();
			format.set_storage_to(rows[next_squad]);
			iome[this_squad].unlock();
			this_squad=next_squad;
		}
	}
	iome[this_squad].unlock();
	worker.join();
	if (!pqq_fn.empty())
	{
		qq_plot(pqq_pv, pqq_fn, 8, "pdf", perch::gnuplot, false);
		if (savePV) {
			openOutFile_or_exit(pv_out,pqq_fn+"_pv.txt");
			pv_out << "test_pv\t%cov_ct\t%cov_cs\tGQratio\tDPratio\tMR_pv\tMisRate\tHW_pv\tObsA_pv\tObsU_pv\t#CHROM\tPOS\tID\tREF\tALT\t"<<perch::h_MxAF<<"\tSplAF" << endl;
			pv_out << pqq_log.str();
			closefile(pv_out);
		}
	}
	if (par.sum_RV)
	{
		double tot_spl = std::max(par.TotSpl, num_fd);
		if (tot_spl)
		{
			lns<<showl<<"Within the minimum exome of autosomal chromosomes:";
			lns<<showl<<"  Number of variants per founder: "<<(double)fnd_RV[0]/tot_spl<<" rare; "<<(double)fnd_SV[0]/tot_spl<<" singletons; and "<<(double)fnd_PV[0]/tot_spl<<" personal.";
			lns<<showl<<"  Number of SNVs     per founder: "<<(double)fnd_RS[0]/tot_spl<<" rare; "<<(double)fnd_SS[0]/tot_spl<<" singletons; and "<<(double)fnd_PS[0]/tot_spl<<" personal.";
		}
		if (num_ct)
		{
			lns<<showl<<"  Number of variants per control: "<<(double)fnd_RV[1]/num_ct<<" rare; "<<(double)fnd_SV[1]/num_ct<<" singletons; and "<<(double)fnd_PV[1]/num_ct<<" personal.";
			lns<<showl<<"  Number of SNVs     per control: "<<(double)fnd_RS[1]/num_ct<<" rare; "<<(double)fnd_SS[1]/num_ct<<" singletons; and "<<(double)fnd_PS[1]/num_ct<<" personal.";
		}
		if (num_cs)
		{
			lns<<showl<<"  Number of variants per case:    "<<(double)fnd_RV[2]/num_cs<<" rare; "<<(double)fnd_SV[2]/num_cs<<" singletons; and "<<(double)fnd_PV[2]/num_cs<<" personal.";
			lns<<showl<<"  Number of SNVs     per case:    "<<(double)fnd_RS[2]/num_cs<<" rare; "<<(double)fnd_SS[2]/num_cs<<" singletons; and "<<(double)fnd_PS[2]/num_cs<<" personal.";
		}
		lns<<flush_logger;
	}
	if (!RepSpl.empty() && dup_fn.empty())
	{
		lns<<showl<<"#\tBothMs\tRefMs\tSoftMs\tHardMs\tConcord\t0->1,2\t1,2->0\t1->2\t2->1\tRefNoMs\tMsRate\tBthNoMs\tConc\tSens\tSpec\tPPV\tNPV";
		for (size_t i=1;i<false_pos.size();++i)
		{
			lns<<showl<<"# col"<<i+1<<'\t'<<both12_na[i]<<'\t'<<SplSet1na[i]<<'\t'<<soft_miss[i]<<'\t'<<hard_miss[i]<<'\t'<<concord_n[i]+concord_p[i]<<'\t'<<false_pos[i]<<'\t'<<false_neg[i]<<'\t'<<het_2_hom[i]<<'\t'<<hom_2_het[i];
			// below calculate the missing rate among the genotype-pairs that Ref are non-missing.
			int ref_non_mss = soft_miss[i]+hard_miss[i]+concord_n[i]+concord_p[i]+false_pos[i]+false_neg[i]+het_2_hom[i]+hom_2_het[i];
			lns<<'\t'<< ref_non_mss;
			lns<<'\t'<< ftos_MaxWidth((soft_miss[i]+hard_miss[i])/(double)ref_non_mss);
			// below calculate the condordance rate among the genotype-pairs that both are non-missing. RefAlt and AltAlt are different genotypes.
			int both_non_mss = concord_n[i]+concord_p[i]+false_pos[i]+false_neg[i]+het_2_hom[i]+hom_2_het[i];
			lns<<'\t'<< both_non_mss;
			lns<<'\t'<< ftos_MaxWidth((concord_n[i]+concord_p[i])/(double)both_non_mss);
			// below are summary statistics that also consider missing values, but doesn't differentiate RefAlt and AltAlt.
			double TP = concord_p[i]+het_2_hom[i]+hom_2_het[i];
			double TN = concord_n[i]+soft_miss[i];
			double FP = false_pos[i];
			double FN = false_neg[i]+hard_miss[i];
			lns<<'\t'<< ftos_MaxWidth(TP/(TP+FN));
			lns<<'\t'<< ftos_MaxWidth(TN/(TN+FP));
			lns<<'\t'<< ftos_MaxWidth(TP/(TP+FP));
			lns<<'\t'<< ftos_MaxWidth(TN/(TN+FN));
		}
		lns<<flush_logger;
	}
	if (!vqs_fn.empty())
	{
		openOutFile_or_exit(vqs_out,vqs_fn);
		vqs_out << vqs_log.rdbuf();
		closefile(vqs_out);
	}
	if (!dnv_fn.empty())
	{
		openOutFile_or_exit(dnv_out,dnv_fn);
		dnv_out << dnv_log.rdbuf();
		closefile(dnv_out);
	}
	if (!par.spl_qc.empty())
	{
		openOutFile_or_exit(sqc_out,par.spl_qc);
		sqc_out<<sample_qc_header_str<<endl;
		int totXhet=0, totXcall=0, mnY=0, mnYcall=0, wmY=0, wmYcall=0;
		for (auto &spl:samples)
		{
			int Y_tot_var = spl.msy+spl.nmy;
			double Y_call_rate = (Y_tot_var>0 ? (double)spl.nmy/Y_tot_var : std::numeric_limits<double>::signaling_NaN());
			int P_sex = spl.sex;															// provided sex
			int Y_sex = ( (spl.nmy>=par.y_n_mn && Y_call_rate>par.y_r_mn) ? 1 : ( (Y_tot_var>=par.y_n_wm && Y_call_rate<par.y_r_wm) ? 2 : 0) );	// chr Y
			int X_sex = (spl.lrx > par.xSxLOD ? 2 : (spl.lrx<-par.xSxLOD ? 1 : 0) );		// chr X
			int G_sex = ((X_sex && Y_sex && X_sex!=Y_sex) ? 0 : std::max(X_sex,Y_sex));		// genetic sex
			int F_sex = ((G_sex && P_sex && G_sex!=P_sex) ? 0 : std::max(G_sex,P_sex));		// final sex
			sqc_out << spl.ID << DLMTR << spl.sex << DLMTR << spl.mss << DLMTR << spl.nms << DLMTR << spl.het << DLMTR << spl.aho
			<< DLMTR << spl.rhe << DLMTR << spl.rho << DLMTR << spl.nes << DLMTR << spl.nos << DLMTR << spl.oes << DLMTR << spl.oos
			<< DLMTR << spl.eti << DLMTR << spl.etv
			<< DLMTR << ((spl.mss+spl.nms)>0 ? (double)spl.mss/(spl.mss+spl.nms) : -1)
			<< DLMTR << (spl.aho>0 ? (double)spl.het/spl.aho : -1)
			<< DLMTR << (spl.nms>0 ? (double)spl.het/spl.nms : -1)
			<< DLMTR << (spl.etv>0 ? (double)spl.eti/spl.etv : -1)
			<< DLMTR << spl.msy << DLMTR << spl.nmy << DLMTR << (std::isnan(Y_call_rate) ? -1 : Y_call_rate) << DLMTR << Y_sex
			<< DLMTR << spl.nmx << DLMTR << spl.n1x << DLMTR << spl.lrx << DLMTR << X_sex
			<< DLMTR << spl.rnm << DLMTR << spl.rbg << DLMTR << spl.rcc << DLMTR << (spl.rnm ? 100.0*spl.rbg/spl.rnm : -1) << DLMTR << (spl.rbg ? 100.0*spl.rcc/spl.rbg : -1)
			<< DLMTR << spl.nGQ << DLMTR << spl.nDP << DLMTR << (spl.nGQ?spl.tGQ/spl.nGQ:0) << DLMTR << (spl.nDP?spl.tDP/spl.nDP:0)
			<< endl;
			if		(X_sex && Y_sex && X_sex!=Y_sex) lns<<showl<<spl.ID<<"\tX_sex_different_from_Y_sex"<<flush_logger;
			else if (G_sex && P_sex && P_sex!=G_sex) lns<<showl<<spl.ID<<"\tGenetic_sex_different_from_the_provided_sex"<<flush_logger;
			else if (F_sex==1) { mnY+=Y_tot_var; mnYcall+=spl.nmy; totXhet+=spl.n1x; totXcall+=spl.nmx; }
			else if (F_sex==2) { wmY+=Y_tot_var; wmYcall+=spl.nmy; }
		}
		if (totXcall) sqc_out<<"##X-het-hap-rate: "<<(double)totXhet/totXcall<<endl;
		if (mnY) sqc_out<<"##Ycall%-mn: "<<(double)mnYcall/mnY<<endl;
		if (wmY) sqc_out<<"##Ycall%-wm: "<<(double)wmYcall/wmY<<endl;
		closefile(sqc_out);
	}
	return 0;
}
